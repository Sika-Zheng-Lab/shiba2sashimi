import sys
import os
import logging
# Configure logging
logger = logging.getLogger(__name__)

def extract_junctions_in_region(shiba_path, chrom, start, end, junction_list = None) -> dict:
	"""
	Get read number for each junction in the specified region from Shiba output.
	"""
	junctions_bed = os.path.join(shiba_path, "junctions", "junctions.bed")
	# Check if junctions.bed file exists
	if not os.path.exists(junctions_bed):
		logger.error(f"Junctions file not found: {junctions_bed}")
		logger.error("Please double check and provide a valid Shiba output path")
		raise FileNotFoundError(f"Junctions file not found: {junctions_bed}")
		sys.exit(1)
	# Initialize junction dictionary
	junctions_dict = {}
	# Read junctions.bed file
	with open(junctions_bed, "r") as junctions:
		for line in junctions:
			line = line.strip()
			# Skip header
			if line.startswith("chr\tstart"):
				samples = line.split("\t")[4:]
				samples_col_dict = {sample: i for i, sample in enumerate(samples)}
				continue
			# Parse junction information
			junc_cols = line.split("\t")[:4]
			junc_chrom = junc_cols[0]
			junc_start = junc_cols[1]
			junc_end = junc_cols[2]
			junc_ID = junc_cols[3]
			junction_values = line.split("\t")[4:]
			if junction_list:
				if junc_ID not in junction_list:
					continue
				for sample, col in samples_col_dict.items():
					if sample not in junctions_dict:
						junctions_dict[sample] = {}
					junctions_dict[sample][junc_ID] = int(junction_values[col])
			else:
				# Check if junction is within the specified region
				if (junc_chrom == f"chr{chrom}" or junc_chrom == chrom) and (start < int(junc_start) < end) and (start < int(junc_end) < end):
					# Get read number
					for sample, col in samples_col_dict.items():
						if sample not in junctions_dict:
							junctions_dict[sample] = {}
						junctions_dict[sample][junc_ID] = int(junction_values[col])
	return junctions_dict
