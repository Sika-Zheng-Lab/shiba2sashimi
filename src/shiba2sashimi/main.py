import argparse
import sys
import logging
import time
from . import tables, bams, plots, junc, utils
# Configure logger
logger = logging.getLogger(__name__)
# Set version
VERSION = "v0.1.4"

def parse_args():
	parser = argparse.ArgumentParser(
		description=f"shiba2sashimi {VERSION} - Create Sashimi plot from Shiba output"
	)

	parser.add_argument("-e", "--experiment", required = True, help = "Experiment table used for Shiba")
	parser.add_argument("-s", "--shiba", required = True, help = "Shiba output directory")
	parser.add_argument("-o", "--output", required = True, help = "Output file")
	parser.add_argument("--id", required = False, help = "Positional ID (pos_id) of the event to plot")
	parser.add_argument("-c", "--coordinate", required = False, help = "Coordinates of the region to plot")
	parser.add_argument("--samples", required = False, help = "Samples to plot. e.g. sample1,sample2,sample3 Default: all samples in the experiment table")
	parser.add_argument("--groups", required = False, help = "Groups to plot. e.g. group1,group2,group3 Default: all groups in the experiment table. Overrides --samples")
	parser.add_argument("--colors", required = False, help = "Colors for each group. e.g. red,orange,blue")
	parser.add_argument("--extend_up", default = 500, type = int, help = "Extend the plot upstream. Only used when not providing coordinates. Default: %(default)s")
	parser.add_argument("--extend_down", default = 500, type = int, help = "Extend the plot downstream. Only used when not providing coordinates. Default: %(default)s")
	parser.add_argument("--smoothing_window_size", default = 21, type = int, help = "Window size for median filter to smooth coverage plot. Greater value gives smoother plot. Default: %(default)s")
	parser.add_argument("--font_family", help = "Font family for labels")
	parser.add_argument("--dpi", default = 300, type = int, help = "DPI of the output figure. Default: %(default)s")
	parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase verbosity")
	args = parser.parse_args()
	return args

def main():

	# Get arguments
	args = parse_args()

	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)

	# Validate input and config
	logger.info(f"Running shiba2sashimi ({VERSION})")
	time.sleep(1)
	logger.debug(f"Arguments: {args}")

	# Load experiment table
	logger.info(f"Loading experiment table from {args.experiment}")
	experiment_dict = tables.load_experiment_table(args.experiment)
	logger.debug(f"Experiment table: {experiment_dict}")

	# Check if provided samples exist in the experiment table
	if args.samples:
		for sample in args.samples.split(","):
			if sample not in experiment_dict:
				logger.error(f"Sample not found in the experiment table: {sample}")
				logger.error("Please double check and provide a valid sample")
				sys.exit(1)
	# Check if provided groups exist in the experiment table
	if args.groups:
		for group in args.groups.split(","):
			if group not in set([info["group"] for info in experiment_dict.values()]):
				logger.error(f"Group not found in the experiment table: {group}")
				logger.error("Please double check and provide a valid group")
				sys.exit(1)

	# Get PSI values from Shiba output
	if args.id:
		logger.debug(f"Get PSI values for positional ID: {args.id}")
		psi_values_dict = tables.get_psi_values(args.id, args.shiba)

	# Get coordinates of the target region from positional ID or coordinate
	if args.id:
		logger.debug(f"Extracting coordinates from positional ID: {args.id}")
		# Get coordinates from positional ID
		chrom, start, end, strand, gene_name, junction_list, junction_direction_dict = utils.posid2int(args.id, args.shiba, args.extend_up, args.extend_down)
		logger.debug(f"junction_list: {junction_list}")
		if args.coordinate:
			logger.debug(f"Using provided coordinate: {args.coordinate}")
			# Get coordinates from provided coordinate
			chrom, start, end = utils.coord2int(args.coordinate)
	elif args.coordinate:
		logger.debug(f"Using provided coordinate: {args.coordinate}")
		# Get coordinates from provided coordinate
		chrom, start, end = utils.coord2int(args.coordinate)
		strand = None
		gene_name = None
		junction_list = None
		junction_direction_dict = None
	else:
		logger.error("Please provide either positional ID or coordinate to define the target region")
		sys.exit(1)

	# Get coverage of the target region for each sample
	logger.info("Calculating coverage for each sample")
	coverage_dict = {}
	if args.groups and args.samples:
		logger.info("Both --samples and --groups are provided. Ignoring --samples")
	for sample, info in experiment_dict.items():
		if args.groups:
			if info["group"] not in args.groups.split(","):
				continue
		elif args.samples:
			if sample not in args.samples.split(","):
				continue
		logger.info(f"Calculating coverage for {sample}")
		window_size = args.smoothing_window_size if args.smoothing_window_size % 2 == 1 else args.smoothing_window_size + 1
		coverage = bams.get_coverage(info["bam"], chrom, start, end, window_size)
		coverage_dict[sample] = coverage

	# Get information of target junctions
	logger.info("Extracting junctions in the target region")
	logger.debug(f"Target region: {chrom}:{start}-{end}")
	junctions_dict = junc.extract_junctions_in_region(args.shiba, chrom, start, end, junction_list)
	logger.debug(f"Junctions in the target region: {junctions_dict}")

	# Create Sashimi plot
	logger.info("Creating Sashimi plot")
	plots.sashimi(
		coverage_dict = coverage_dict,
		junctions_dict = junctions_dict,
		experiment_dict = experiment_dict,
		samples = args.samples if not args.groups else ",".join(coverage_dict.keys()),
		groups = args.groups,
		colors = args.colors,
		chrom = chrom,
		start = start,
		end = end,
		output = args.output,
		pos_id = args.id if args.id else None,
		coordinate = args.coordinate if args.coordinate else None,
		strand = strand,
		gene_name = gene_name,
		junction_direction_dict = junction_direction_dict if args.id else None,
		psi_values_dict = psi_values_dict if args.id else None,
		font_family = args.font_family if args.font_family else None,
		dpi = args.dpi
	)

	# Finish
	logger.info("shiba2sashimi finished successfully")
	logger.info(f"Output file: {args.output}")

	# Return
	return 0

if __name__ == "__main__":
	main()

