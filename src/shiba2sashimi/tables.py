import os
import sys
import logging
# Configure logging
logger = logging.getLogger(__name__)

def load_experiment_table(experiment_table) -> dict:
	"""
	Load experiment table from file

	Parameters
	- experiment_table: str
		Path to experiment table file

	Returns
	- experiment_dict: dict
		Dictionary of experiment table
		{sample: {bam, group, bam_index}}
	"""
	experiment_dict = {}
	# Check if the experiment table file exists
	if not os.path.exists(experiment_table):
		logger.error(f"Experiment table not found: {experiment_table}")
		logger.error("Please double check and provide a valid experiment table")
		raise FileNotFoundError(f"Experiment table not found: {experiment_table}")
		sys.exit(1)
	# Load experiment table
	with open(experiment_table, "r") as experiment:
		for line in experiment:
			line = line.strip()
			if not line or line.startswith("sample"):
				continue
			sample, bam, group = line.split(maxsplit=2)
			experiment_dict[sample] = {
				"bam": bam,
				"group": group,
			}
	return experiment_dict

def get_psi_values(positional_id, shiba_path) -> dict:
	"""
	Load PSI values from PSI file

	Parameters
	- positional_id: str
		Positional ID (e.g. SE@chr1:1000-2000)
	- shiba_path: str
		Path to Shiba output directory

	Returns
	- psi_values_dict: dict
		Dictionary of PSI values
		{sample: psi}
	"""
	psi_matrix_path = os.path.join(shiba_path, "results", "splicing", "PSI_matrix_sample.txt")
	# Check if PSI matrix file exists
	if not os.path.exists(psi_matrix_path):
		logger.error(f"PSI matrix file not found: {psi_matrix_path}")
		logger.error("Please double check and provide a valid Shiba output directory")
		raise FileNotFoundError(f"PSI matrix file not found: {psi_matrix_path}")
		sys.exit(1)
	# Load PSI matrix file
	psi_values_dict = {}
	with open(psi_matrix_path, "r") as psi_matrix:
		for line in psi_matrix:
			line = line.strip()
			if not line or line.startswith("event_id"):
				samples = line.split("\t")[2:]
				continue
			pos_id_in_file = line.split("\t")[1]
			if pos_id_in_file == positional_id:
				psi_values = line.split("\t")[2:]
				for sample, psi in zip(samples, psi_values):
					if psi.replace('.', '', 1).isdigit():
						psi_values_dict[sample] = float(psi) * 100
					else:
						psi_values_dict[sample] = "NA"
				break
	if not psi_values_dict:
		# PSI values are all NA
		psi_values_dict = {sample: "NA" for sample in samples}
	return psi_values_dict
