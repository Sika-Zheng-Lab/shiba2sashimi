import sys
import os
import numpy as np
import logging
# Configure logging
logger = logging.getLogger(__name__)
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from matplotlib import font_manager

def bezier_point(t, p0, p1, p2, p3):
	return (
		(1 - t)**3 * np.array(p0) +
		3 * (1 - t)**2 * t * np.array(p1) +
		3 * (1 - t) * t**2 * np.array(p2) +
		t**3 * np.array(p3)
	)

def sashimi(coverage_dict, junctions_dict, experiment_dict, samples, groups, colors, chrom, start, end, output, pos_id = None, coordinate = None, strand = None, gene_name = None, junction_direction_dict = None, psi_values_dict = None, font_family = None, dpi = 300):
	"""
	Create Sashimi plot.
	"""
	# Make sure that fonts can be found in a Docker/Singularity container
	font_dir = '/usr/share/fonts/truetype/msttcorefonts/'
	if os.path.exists(font_dir):
		for font_file in os.listdir(font_dir):
			if font_file.endswith('.ttf'):
				font_manager.fontManager.addfont(os.path.join(font_dir, font_file))
	# Set font family
	if font_family:
		matplotlib.rcParams["font.family"] = font_family
	chrom = f"chr{chrom}" if not chrom.startswith("chr") and (chrom.isdigit() or chrom in ["X", "Y", "M", "MT"]) else chrom
	# Set figure size
	n_samples = len(coverage_dict)
	fig_height = 1 * n_samples
	fig_width = 8
	fig = plt.figure(figsize=(fig_width, fig_height))
	# Subplots for coverage
	gs = fig.add_gridspec(n_samples + 1, 1, hspace=0.8, height_ratios=[1] * n_samples + [0.2])
	# Set sample order
	sample_order = []
	if groups:
		groups_list = groups.split(",")
		for group in groups_list:
			sample_group_order = [sample for sample, info in experiment_dict.items() if info["group"] == group]
			if samples:
				sample_group_order = sorted(sample_group_order, key=samples.split(",").index)
			sample_order += sample_group_order
	else:
		groups_list = []
		seen = set()
		for info in experiment_dict.values():
			group = info["group"]
			if group not in seen:
				groups_list.append(group)
				seen.add(group)
		if samples:
			sample_order = samples.split(",")
		else:
			sample_order = list(experiment_dict.keys())
	# Set colors for each group
	colors_list = colors.split(",") if colors else ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"]
	if len(colors_list) < len(groups_list):
		logger.error("Number of colors is less than number of groups")
		sys.exit(1)
	try:
		color_dict = {group: color for group, color in zip(groups_list, colors_list)}
	except ValueError:
		logger.error("Number of colors does not match number of groups")
		sys.exit(1)
	# Plot coverage for each sample
	junc_reads_max = max([max([junc_reads for junc_ID, junc_reads in region_junctions.items()]) for region_junctions in junctions_dict.values()])
	junc_reads_min = min([min([junc_reads for junc_ID, junc_reads in region_junctions.items()]) for region_junctions in junctions_dict.values()])
	for i, sample_name in enumerate(sample_order):
		ax = fig.add_subplot(gs[i, 0])
		cov = coverage_dict[sample_name]
		cov_max = max(cov)
		x_positions = range(start, end)
		group = experiment_dict[sample_name]["group"]
		color = color_dict[group]
		ax.fill_between(x_positions, cov, step="pre", color=color, alpha=0.8)
		# Add sample name and PSI value
		if psi_values_dict:
			try:
				psi = psi_values_dict[sample_name]
				ax.text(0.01, 0.85, f"{sample_name} (PSI = {psi:.2f})",transform=ax.transAxes, fontsize=10, color="black")
			except:
				psi = "NA"
				ax.text(0.01, 0.85, f"{sample_name} (PSI = {psi})", transform=ax.transAxes, fontsize=10, color="black")
		else:
			ax.text(0.01, 0.85, f"{sample_name}", transform=ax.transAxes, fontsize=10, color="black")
		# Plot junctions
		region_junctions = junctions_dict[sample_name]
		for junc_ID in region_junctions:
			# Get direction of junction
			direction = junction_direction_dict[junc_ID] if junction_direction_dict else "up"
			# Get junction coordinates
			junc_start = int(junc_ID.split(":")[1].split("-")[0]) - 1 # 0-based
			junc_end = int(junc_ID.split(":")[1].split("-")[1])
			# Get number of reads
			junc_reads = region_junctions[junc_ID]
			# Ignore if junction is out of range
			if not (start < junc_start < end and start < junc_end < end):
				continue
			# Draw arc
			(x1, y1) = (junc_start, cov[junc_start - start]) if direction == "up" else (junc_start, 0)
			(x2, y2) = (junc_end, cov[junc_end - start]) if direction == "up" else (junc_end, 0)
			# Calculate control point for quadratic Bezier curve
			arc_height = cov_max * 0.5
			ctrl1 = (x1, y1 + arc_height) if direction == "up" else (x1, y1 - arc_height)
			ctrl2 = (x2, y2 + arc_height) if direction == "up" else (x2, y2 - arc_height)
			verts = [ (x1, y1), ctrl1, ctrl2, (x2, y2) ]
			codes = [ Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4 ]
			# Set linewidth according to the number of reads
			linewidth_factor = (1.5 - 0.5) / (junc_reads_max - junc_reads_min) if junc_reads_max != junc_reads_min else 1  # Scale linewidth from 0.5 to 1.5
			arc_linewidth = 0.5 + (junc_reads - junc_reads_min) * linewidth_factor
			if junc_reads == 0:
				arc_linewidth = 0.25
			# Create a Bezier curve patch
			path = Path(verts, codes)
			bezier = PathPatch(path, linewidth=arc_linewidth, edgecolor=color, facecolor='none', clip_on=False)
			ax.add_patch(bezier)
			# Calculate midpoint (to use as the center of the arc)
			bx, by = bezier_point(0.5, (x1, y1), ctrl1, ctrl2, (x2, y2))
			# Add junc_reads as text on the arc
			ax.text(
				bx, by, str(junc_reads),
				fontsize=8, ha='center', va='center', color='black',
				backgroundcolor='white', bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0'),
				clip_on=False  # Allow text to be drawn outside the axes
			)
		ax.set_xlim(start, end)
		ax.set_ylim(bottom = 0, top = max(cov) * 1.4)
		ax.set_ylabel("Coverage", fontsize=8)
		# Despine top, right, and bottom
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		# Remove xticks for all samples
		ax.set_xticks([])
	# Create a separate x-axis at the bottom
	ax_x = fig.add_subplot(gs[-1, 0])
	ax_x.set_xlim(start, end)
	ax_x.xaxis.set_major_locator(MaxNLocator(nbins=10, integer=True)) # Set number of ticks
	# Disable scientific notation
	formatter = ScalarFormatter(useOffset=False, useMathText=False)
	formatter.set_scientific(False)
	ax_x.xaxis.set_major_formatter(formatter)
	ax_x.tick_params(axis='x', labelrotation=45, labelsize=6)
	for label in ax_x.get_xticklabels():
		label.set_ha('right')  # Set horizontal alignment to 'right'
	ax_x.set_xlabel(f"Genomic coordinate ({chrom})", fontsize=10)
	ax_x.spines['top'].set_visible(False)
	ax_x.spines['right'].set_visible(False)
	ax_x.spines['left'].set_visible(False)
	ax_x.get_yaxis().set_visible(False)
	# Put title on the top subplot
	title = ""
	if coordinate:
		title += f"{chrom}:{start}-{end}"
	if pos_id:
		title += f"\n{pos_id}, {gene_name} ({strand})"
	if title:
		ax_top = fig.axes[0]
		ax_top.annotate(
			title,
			xy=(0.5, 1.5),
			xycoords='axes fraction',
			ha='center',
			va='bottom',
			fontsize=12
		)
	# Save plot
	plt.savefig(output, dpi=dpi, bbox_inches="tight")
