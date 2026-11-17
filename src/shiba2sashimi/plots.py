"""Sashimi rendering with aligned transcript structures."""
from contextlib import contextmanager
import logging
import textwrap
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator, ScalarFormatter

from . import bams, tracks
from .transcripts import event_features
from .utils import coord2int

logger = logging.getLogger(__name__)
PALETTE = ["#767676", "#3775A0", "#B56363", "#57917D", "#8A739C", "#B28D4C",
           "#467D87", "#A57785", "#6C7B4D", "#716F9A", "#A16F48", "#4D4D4D"]
EVENT_COLOR = "#B64342"


@contextmanager
def figure_style(font_family=None):
    family = font_family or "Arial"
    # Do not silently replace a requested publication font.
    try:
        font_manager.findfont(font_manager.FontProperties(family=family), fallback_to_default=False)
    except ValueError:
        raise ValueError(f"Font {family!r} is not installed. Install it or select an installed font with --font_family.") from None
    with mpl.rc_context({"font.family": family, "font.size": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "axes.linewidth": .7, "xtick.major.width": .7, "ytick.major.width": .7,
                         "xtick.direction": "out", "ytick.direction": "out",
                         "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42}):
        yield


def bezier_point(t, p0, p1, p2, p3):
    return ((1-t)**3 * np.array(p0) + 3*(1-t)**2*t * np.array(p1)
            + 3*(1-t)*t**2 * np.array(p2) + t**3 * np.array(p3))


def _format_count(count):
    """Fixed-point labels with up to three decimals, never scientific notation."""
    return f"{count:.3f}".rstrip("0").rstrip(".")


def _draw_models(ax, models, features, start, end):
    for i, model in enumerate(models):
        y = len(models) - i - 1
        ax.text(0, y + .28, model["name"], transform=ax.get_yaxis_transform(), fontsize=7, va="bottom")
        for left, right in model["introns"]:
            left, right = max(start, left), min(end, right)
            if left >= right:
                continue
            ax.plot([left, right], [y, y], lw=.7, color="#767676")
            if model["strand"] in ("+", "-"):
                marker = ">" if model["strand"] == "+" else "<"
                positions = np.linspace(left, right, 5)[1:-1]
                ax.plot(positions, [y]*len(positions), linestyle="none", marker=marker,
                        markersize=2.5, color="#767676")
        for left, right in model["exons"]:
            left, right = max(start, left), min(end, right)
            if left >= right:
                continue
            ax.add_patch(Rectangle((left, y-.12), right-left, .24, facecolor="#767676", linewidth=0))
            for a, b in features:
                a, b = max(a, left), min(b, right)
                if a < b:
                    ax.add_patch(Rectangle((a, y-.12), b-a, .24, facecolor=EVENT_COLOR, linewidth=0))
    ax.set_ylim(-.5, len(models) - .15)
    ax.set_xlim(start, end)
    ax.set_axis_off()


def sashimi(coverage_dict, junctions_dict, experiment_dict, samples, groups, colors,
            fig_width, chrom, start, end, output, pos_id=None, coordinate=None,
            strand=None, gene_name=None, junction_direction_dict=None, psi_values_dict=None,
            font_family=None, dpi=300, nolabel=False, nojunc=False, minimum_junc_reads=1,
            group_average=False, smoothing_window_size=1, transcript_models=None, event=None):
    """Render sample tracks or means; return the closed figure for inspection/tests."""
    order = tracks.select_samples(experiment_dict, samples, groups)
    data = tracks.build_tracks(coverage_dict, junctions_dict, experiment_dict, order,
                               psi_values_dict, group_average)
    if not data:
        raise ValueError("No samples selected")
    group_order = list(dict.fromkeys(track["group"] for track in data))
    palette = [c.strip() for c in colors.split(",")] if colors else PALETTE
    if len(palette) < len(group_order):
        raise ValueError("Number of colors is less than number of groups")
    color_map = dict(zip(group_order, palette))
    for color in color_map.values():
        if not mpl.colors.is_color_like(color):
            raise ValueError(f"Invalid color: {color}")
    for track in data:
        track["coverage"] = bams.median_filter(track["coverage"], smoothing_window_size)
        if len(track["coverage"]) != end-start:
            raise ValueError("Coverage length does not match the plotting region")
    models = transcript_models or []
    features = event_features(event)
    model_height = .34 * len(models) + .15 if models else 0
    heights = [1.1]*len(data) + ([model_height] if model_height else []) + [.14]
    with figure_style(font_family):
        fig = plt.figure(figsize=(fig_width, sum(heights) + .8))
        try:
            grid = fig.add_gridspec(len(heights), 1, height_ratios=heights, hspace=.42)
            all_counts = [count for track in data for j, count in track["junctions"].items()
                          if count > 0 and count >= minimum_junc_reads
                          and coord2int(j)[2] - coord2int(j)[1] > 1]
            min_count, max_count = (min(all_counts), max(all_counts)) if all_counts else (0, 0)
            # Shared scales allow direct comparison, including groups of unequal size.
            scale = max(1., max(float(np.max(track["coverage"])) for track in data))
            has_down = not nojunc and any(
                (junction_direction_dict or {}).get(j) == "down" and c > 0 and c >= minimum_junc_reads
                for track in data for j, c in track["junctions"].items())
            boundary_positions = sorted({coord2int(j)[1]
                for track in data for j, c in track["junctions"].items()
                if not nojunc and c > 0 and c >= minimum_junc_reads
                and coord2int(j)[2] - coord2int(j)[1] == 1
                and start <= coord2int(j)[1] <= end})
            # Place labels outward when a short retained intron brings them close.
            boundary_alignment = dict.fromkeys(boundary_positions, "center")
            if len(boundary_positions) == 2 and (boundary_positions[1]-boundary_positions[0])/(end-start) * fig_width < 1.2:
                boundary_alignment = dict(zip(boundary_positions, ["right", "left"]))
            for i, track in enumerate(data):
                ax = fig.add_subplot(grid[i])
                cov, color = track["coverage"], color_map[track["group"]]
                # Base intervals are [x, x+1), matching GTF rectangles and junction boundaries.
                ax.fill_between(np.arange(start, end+1), np.r_[cov, cov[-1]], step="post", color=color, lw=0)
                if not nolabel:
                    label = track["label"]
                    if group_average:
                        label += f" (n = {track['n']})"
                    if track["has_psi"]:
                        value = f"{track['psi']:.2f}%" if track["psi"] is not None else "NA"
                        label += f"; {'mean ' if group_average else ''}PSI = {value}"
                        if group_average:
                            label += f" (n = {track['psi_n']})"
                    ax.text(0, 1.02, label, transform=ax.transAxes, fontsize=8, va="bottom")
                if not nojunc:
                    for junction_id, count in track["junctions"].items():
                        if count <= 0 or count < minimum_junc_reads:
                            continue
                        _, left, right = coord2int(junction_id)
                        if right - left == 1:
                            if left in boundary_alignment:
                                y = -scale * .16
                                align = boundary_alignment[left]
                                ax.plot([left, left], [0, y], color=color, lw=.7, linestyle="--")
                                ax.plot(left, y, marker="_", markersize=5, color=color)
                                ax.annotate(f"EI: {_format_count(count)}", (left, y),
                                            xytext=({"right": -3, "left": 3, "center": 0}[align], -3),
                                            textcoords="offset points", ha=align, va="top", fontsize=7)
                            continue
                        # Shiba flanking exon bases -> zero-based splice boundaries.
                        x1, x2 = left, right-1
                        if not (start <= x1 < x2 <= end):
                            continue
                        direction = (junction_direction_dict or {}).get(junction_id, "up")
                        y1 = cov[max(0, x1-start-1)] if direction == "up" else 0
                        y2 = cov[min(len(cov)-1, x2-start)] if direction == "up" else 0
                        offset = scale*.55 * (1 if direction == "up" else -1)
                        points = [(x1, y1), (x1, y1+offset), (x2, y2+offset), (x2, y2)]
                        width = .6 + (count-min_count)/(max_count-min_count) if max_count > min_count else 1.
                        ax.add_patch(PathPatch(Path(points, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                                               lw=width, edgecolor=color, facecolor="none"))
                        bx, by = bezier_point(.5, *points)
                        label = _format_count(count)
                        ax.text(bx, by, label, ha="center", va="center", fontsize=7,
                                bbox={"facecolor": "white", "edgecolor": "none", "pad": .4})
                ax.set(xlim=(start, end), ylim=(-scale*.65 if boundary_positions else -scale*.55 if has_down else 0, scale*1.65),
                       ylabel="Mean depth" if group_average else "Read depth")
                ax.set_yticks([0, scale])
                ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: _format_count(value)))
                ax.spines["left"].set_bounds(0, scale)
                for side in ("top", "right", "bottom"):
                    ax.spines[side].set_visible(False)
                ax.set_xticks([])
            if model_height:
                _draw_models(fig.add_subplot(grid[-2]), models, features, start, end)
            ax_x = fig.add_subplot(grid[-1])
            ax_x.set_xlim(start, end)
            ax_x.xaxis.set_major_locator(MaxNLocator(nbins=max(2, int(fig_width)), integer=True))
            formatter = ScalarFormatter(useOffset=False)
            formatter.set_scientific(False)
            ax_x.xaxis.set_major_formatter(formatter)
            ax_x.set_xlabel(f"Genomic coordinate ({chrom}; 0-based)")
            ax_x.get_yaxis().set_visible(False)
            for side in ("top", "right", "left"):
                ax_x.spines[side].set_visible(False)
            title = f"{gene_name or ''} ({strand}) | {pos_id}" if pos_id else f"{chrom}:{start}-{end}"
            title = textwrap.fill(title, width=max(30, int(fig_width * 16)), break_long_words=True)
            fig.suptitle(title, fontsize=8, y=.99)
            title_inches = .22 + .14 * len(title.splitlines())
            fig.subplots_adjust(top=1 - title_inches / fig.get_figheight(),
                                bottom=.13 if sum(heights) < 3 else .08)
            fig.savefig(output, dpi=dpi, bbox_inches="tight")
        finally:
            plt.close(fig)
    return fig
