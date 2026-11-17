"""Command-line entry point."""
import argparse
import logging
from . import tables, bams, plots, junc, utils, tracks, transcripts

logger = logging.getLogger(__name__)
VERSION = "v0.2.0"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=f"shiba2sashimi {VERSION} - Create Sashimi plots from Shiba output")
    parser.add_argument("--version", action="version", version=VERSION)
    parser.add_argument("-e", "--experiment", required=True, help="Shiba bulk experiment table (sample, bam, group; extra columns allowed)")
    parser.add_argument("-s", "--shiba", required=True, help="Shiba working directory")
    parser.add_argument("-o", "--output", required=True, help="Output file (PDF, SVG, PNG, etc.)")
    parser.add_argument("--id", help="Positional ID (pos_id) of the event to plot")
    parser.add_argument("-c", "--coordinate", help="Plot region: chr:start-end (0-based, end-exclusive); overrides the event window")
    parser.add_argument("--samples", help="Comma-separated samples in display order; default: all")
    parser.add_argument("--groups", help="Comma-separated groups in display order; overrides --samples; supports more than two groups")
    parser.add_argument("--group-average", action="store_true", help="One mean track per group; PSI is the mean of nonmissing sample PSI")
    parser.add_argument("--gtf", help="GTF or GTF.gz for transcript models; default: annotation/assembled_annotation.gtf, otherwise event-local structures")
    parser.add_argument("--no-transcripts", action="store_true", help="Hide transcript structures")
    parser.add_argument("--max-transcripts", type=int, default=20, help="Maximum GTF transcripts to show, prioritizing event junction matches (default: %(default)s)")
    parser.add_argument("--colors", help="Comma-separated colors in selected group order")
    parser.add_argument("--width", default=8, type=float, help="Figure width in inches (default: %(default)s)")
    parser.add_argument("--extend_up", default=500, type=int, help="Extend event window upstream in transcript orientation (default: %(default)s)")
    parser.add_argument("--extend_down", default=500, type=int, help="Extend event window downstream in transcript orientation (default: %(default)s)")
    parser.add_argument("--smoothing_window_size", default=1, type=int, help="Median filter window, applied after averaging; 1 disables smoothing, even values round up (default: %(default)s)")
    parser.add_argument("--font_family", default="Arial", help="Installed font family (default: %(default)s)")
    parser.add_argument("--nolabel", action="store_true", help="Hide sample/group labels and PSI")
    parser.add_argument("--nojunc", action="store_true", help="Hide junction arcs and counts, including RI boundaries")
    parser.add_argument("--minimum_junc_reads", default=1, type=float, help="Minimum raw or group-mean count to show a junction or RI boundary (default: %(default)s)")
    parser.add_argument("--dpi", default=300, type=int, help="Raster resolution (default: %(default)s)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase verbosity")
    args = parser.parse_args(argv)
    if not (args.id or args.coordinate):
        parser.error("provide --id or --coordinate")
    for name in ("width", "dpi", "max_transcripts", "smoothing_window_size"):
        value = getattr(args, name)
        if value <= 0 or value == float("inf") or value != value:
            parser.error(f"--{name.replace('_', '-')} must be positive and finite")
    for name in ("extend_up", "extend_down", "minimum_junc_reads"):
        value = getattr(args, name)
        if value < 0 or value == float("inf") or value != value:
            parser.error(f"--{name} must be nonnegative and finite")
    return args


def run(args):
    experiment = tables.load_experiment_table(args.experiment)
    order = tracks.select_samples(experiment, args.samples, args.groups)
    if args.groups and args.samples:
        logger.info("Both --samples and --groups are provided. Ignoring --samples")
    event = utils.read_event(args.id, args.shiba) if args.id else None
    psi = tables.get_psi_values(args.id, args.shiba) if event else None
    if psi is not None:
        missing = set(order) - set(psi)
        if missing:
            raise ValueError(f"Samples missing from PSI matrix: {', '.join(sorted(missing))}")
    if args.coordinate:
        chrom, start, end = utils.coord2int(args.coordinate)
        if event:
            event_chrom, event_start, event_end = utils.event_region(event, 0, 0)
            if not utils.same_chrom(chrom, event_chrom) or end <= event_start or start >= event_end:
                raise ValueError("--coordinate must overlap the selected event on the same chromosome")
    else:
        chrom, start, end = utils.event_region(event, args.extend_up, args.extend_down)
    directions = utils.event_junctions(event) if event else None
    count_ids = list(directions) + utils.event_boundaries(event) if event else None
    junctions = {} if args.nojunc else junc.extract_junctions_in_region(
        args.shiba, chrom, start, end, count_ids)
    if not args.nojunc:
        missing = set(order) - set(junctions)
        if missing:
            raise ValueError(f"Samples missing from junction table: {', '.join(sorted(missing))}")
    models = [] if args.no_transcripts else transcripts.load_models(
        args.shiba, chrom, start, end, event, args.gtf, args.max_transcripts)
    coverage = {}
    for sample in order:
        logger.info("Calculating coverage: %s", sample)
        coverage[sample] = bams.get_coverage(experiment[sample]["bam"], chrom, start, end)
    window = args.smoothing_window_size
    if window % 2 == 0:
        window += 1
        logger.info("Using odd smoothing window: %d", window)
    if args.group_average:
        logger.info("Group means: raw depth/counts over all selected replicates; PSI over valid sample values")
    return plots.sashimi(
        coverage, junctions, experiment, ",".join(order), args.groups, args.colors,
        args.width, chrom, start, end, args.output, pos_id=args.id, coordinate=args.coordinate,
        strand=event["strand"] if event else None,
        gene_name=(event.get("gene_name") or event.get("gene_id")) if event else None,
        junction_direction_dict=directions, psi_values_dict=psi, font_family=args.font_family,
        dpi=args.dpi, nolabel=args.nolabel, nojunc=args.nojunc,
        minimum_junc_reads=args.minimum_junc_reads, group_average=args.group_average,
        smoothing_window_size=window, transcript_models=models,
        event=event if not args.no_transcripts else None)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(format="[%(asctime)s] %(levelname)7s %(message)s",
                        level=logging.DEBUG if args.verbose else logging.INFO)
    logger.info("Running shiba2sashimi (%s)", VERSION)
    try:
        run(args)
    except (OSError, ValueError, KeyError) as error:
        logger.error("%s", error)
        return 1
    logger.info("Output file: %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
