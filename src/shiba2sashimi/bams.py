import sys
import os
import logging
# Configure logging
logger = logging.getLogger(__name__)
import pysam
import numpy as np

def median_filter(data, window_size):
	"""
	Apply a median filter to the data with the specified window size.
	"""
	filtered_data = np.zeros_like(data)
	half_window = window_size // 2
	for i in range(len(data)):
		start = max(0, i - half_window)
		end = min(len(data), i + half_window + 1)
		filtered_data[i] = np.median(data[start:end])
	return filtered_data

def get_coverage(bam_path, chrom, start, end, window_size=21):
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
		try:
			count = bam.count_coverage(chrom, start, end)
		except KeyError:
			logger.debug(f"Chromosome {chrom} not found in BAM file")
			# Remove "chr" prefix and try again
			logger.debug(f"Trying without 'chr' prefix")
			chrom = chrom.replace("chr", "")
			count = bam.count_coverage(chrom, start, end)
		except Exception as e:
			logger.error(f"Failed to get coverage for {chrom}:{start}-{end}")
			logger.error(f"Error: {e}")
			sys.exit(1)
		for i in range(arr_len):
			total = 0
			# Sum up coverage for each base
			for j in range(4):
				total += count[j][i]
			coverage[i] = total
	# Apply median filter for smooth coverage plot
	coverage = median_filter(coverage, window_size)
	return coverage
