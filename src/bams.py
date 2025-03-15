import sys
import os
import logging
# Configure logging
logger = logging.getLogger(__name__)
import pysam
import numpy as np
from scipy.ndimage import median_filter

def coord2int(coordinate):
	"""
	Convert string coordinate to chr, start, and end.
	"""
	try:
		# Split by colon and dash
		chrom, pos = coordinate.split(":")
		start, end = pos.split("-")
		# Convert to integer
		start = int(start)
		end = int(end)
	except:
		logger.error(f"Invalid coordinate format: {coordinate}")
		logger.error("Please provide a valid coordinate format: chr:start-end")
		raise ValueError(f"Invalid coordinate format: {coordinate}")
		sys.exit(1)
	return chrom, start, end

def posid2int(positional_id, shiba_path) -> tuple:
	"""
	Get chromosome, start, and end of the target region and junction coordinates from positional ID.
	"""
	# Parse event type from positional ID (i.e. SE@chr10@20327164-20327391@20326099-20328102)
	event_type = positional_id.split("@")[0] # (SE, FIVE, THREE, MXE, RI, AFE, ALE, MSE)
	# Get PSI file path
	psi_file_path = os.path.join(shiba_path, "results", "splicing", f"PSI_{event_type}.txt")
	# Check if PSI file exists
	if not os.path.exists(psi_file_path):
		logger.error(f"PSI file not found: {psi_file_path}")
		logger.error("Please double check and provide a valid positional ID")
		raise FileNotFoundError(f"PSI file not found: {psi_file_path}")
		sys.exit(1)
	# Read PSI file
	psi_file_col_dict = {}
	header_dict = {}
	with open(psi_file_path, "r") as psi_file:
		for line in psi_file:
			line = line.strip()
			# Skip header
			if line.startswith("event_id"):
				header = line.split("\t")
				for i, col in enumerate(header):
					header_dict[col] = i
				continue
			# Store values in dictionary if the 2nd column matches positional ID
			pos_id_col = line.split("\t")[1]
			if pos_id_col == positional_id:
				# Get all columns and store in dictionary
				columns = line.split("\t")
				for col, i in header_dict.items():
					psi_file_col_dict[col] = columns[i]
				break
	# Check if positional ID exists in PSI file
	if not psi_file_col_dict:
		logger.error(f"Positional ID not found in PSI file: {positional_id}")
		logger.error("Please double check and provide a valid positional ID")
		raise ValueError(f"Positional ID not found in PSI file: {positional_id}")
		sys.exit(1)
	# Get chromosome, start, and end from positional ID
	chrom = positional_id.split("@")[1].replace("chr", "")
	# Extend 1kb upstream and downstream
	intron_key_map = {
		"SE": "intron_c",
		"FIVE": "intron_b",
		"THREE": "intron_b",
		"MXE": ["intron_a1", "intron_a2"],
		"RI": "intron_a",
		"AFE": "intron_a",
		"ALE": "intron_a",
		"MSE": "intron"
	}
	if event_type in intron_key_map:
		intron_keys = intron_key_map[event_type]
		if isinstance(intron_keys, list):
			start = int(psi_file_col_dict[intron_keys[0]].split(":")[1].split("-")[0]) - 1000
			end = int(psi_file_col_dict[intron_keys[1]].split(":")[1].split("-")[1]) + 1000
		else:
			intron = psi_file_col_dict[intron_keys]
			if event_type == "MSE":
				intron = intron.split(";")[-1]
			start = int(intron.split(":")[1].split("-")[0]) - 1000
			end = int(intron.split(":")[1].split("-")[1]) + 1000
		if start < 0:
			start = 0
	else:
		logger.error(f"Unsupported event type: {event_type}")
		logger.error("Please provide a valid positional ID with supported event type")
		raise ValueError(f"Unsupported event type: {event_type}")
		sys.exit(1)

	# Get junction coordinates
	junction_list = []
	junction_key_map = {
		"SE": ["intron_a", "intron_b", "intron_c"],
		"FIVE": ["intron_a", "intron_b"],
		"THREE": ["intron_a", "intron_b"],
		"MXE": ["intron_a1", "intron_a2", "intron_b1", "intron_b2"],
		"RI": ["intron_a"],
		"AFE": ["intron_a", "intron_b"],
		"ALE": ["intron_a", "intron_b"],
		"MSE": "intron"
	}
	junction_keys = junction_key_map[event_type]
	if event_type == "MSE":
		junction_list += psi_file_col_dict[junction_keys].split(";")
	else:
		for key in junction_keys:
			junction_list += [psi_file_col_dict[key]]
	return chrom, start, end, junction_list

def get_coverage(bam_path, chrom, start, end):
	"""
	Return coverage array of the specified region (chrom, start, end) using pysam.
	"""
	# Check if the BAM file and its index exist
	if not os.path.exists(bam_path):
		logger.error(f"BAM file not found: {bam_path}")
		logger.error("Please double check and provide a valid BAM file")
		raise FileNotFoundError(f"BAM file not found: {bam_path}")
		sys.exit(1)
	if not os.path.exists(f"{bam_path}.bai"):
		logger.error(f"BAM index not found: {bam_path}.bai")
		logger.error("Please create index using samtools index")
		raise FileNotFoundError(f"BAM index not found: {bam_path}.bai")
		sys.exit(1)
	# Initialize coverage array
	arr_len = end - start
	coverage = np.zeros(arr_len, dtype=int)
	# Open BAM file with pysam.AlignmentFile
	with pysam.AlignmentFile(bam_path, "rb") as bam:
		# Get coverage using pileup
		count = bam.count_coverage(chrom, start, end)
		for i in range(arr_len):
			total = 0
			# Sum up coverage for each base
			for j in range(4):
				total += count[j][i]
			coverage[i] = total
	# Apply median filter for smooth coverage plot
	window_size = 21
	coverage = median_filter(coverage, size=window_size)
	return coverage

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

