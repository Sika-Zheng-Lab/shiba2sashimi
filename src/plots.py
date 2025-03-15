import sys
import os
import numpy as np
import logging
# Configure logging
logger = logging.getLogger(__name__)
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

def sashimi(coverage_dict, junctions_dict, experiment_dict, samples, groups, colors, chrom, start, end, strand, output, exon_s = 1, intron_s = 1, pos_id = None, psi_values_dict = None, font_family = None):
	"""
	Create Sashimi plot.
	"""
	if font_family:
		matplotlib.rcParams["font.family"] = font_family
	chrom = f"chr{chrom}" if not chrom.startswith("chr") and (chrom.isdigit() or chrom in ["X", "Y", "M", "MT"]) else chrom
	# Set figure size
	n_samples = len(coverage_dict)
	fig_height = 1 * n_samples
	fig_width = 8
	fig = plt.figure(figsize=(fig_width, fig_height))
	# Subplots for coverage
	gs = fig.add_gridspec(n_samples, 1, hspace=0.3)
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
		x_positions = range(start, end)
		group = experiment_dict[sample_name]["group"]
		color = color_dict[group]
		ax.fill_between(x_positions, cov, step="pre", color=color, alpha=0.8)
		# Add sample name and PSI value
		if psi_values_dict:
			psi = psi_values_dict[sample_name]
			ax.text(0.01, 0.85, f"{sample_name} (PSI = {psi:.2f})",transform=ax.transAxes, fontsize=11, color="black")
		else:
			ax.text(0.01, 0.85, f"{sample_name}", transform=ax.transAxes, fontsize=11, color="black")
		# Plot junctions
		region_junctions = junctions_dict[sample_name]
		for junc_ID in region_junctions:
			junc_start = int(junc_ID.split(":")[1].split("-")[0]) - 1 # 0-based
			junc_end = int(junc_ID.split(":")[1].split("-")[1])
			junc_reads = region_junctions[junc_ID]
			# Ignore if junction is out of range
			if not (start < junc_start < end and start < junc_end < end):
				continue
			# Draw arc
			x1, y1 = junc_start, cov[junc_start - start]
			x2, y2 = junc_end, cov[junc_end - start]
			# Calculate midpoint (to use as the center of the arc)
			mx = (x1 + x2) / 2
			my = (y1 + y2) / 2
			# Calculate distance between the two points
			dx = x2 - x1
			dy = y2 - y1
			dist = np.hypot(dx, dy)
			# Angle (in degrees) of the line between the two points
			angle_deg = np.degrees(np.arctan2(dy, dx))
			# Reduce the height of the arc
			def set_arc_height_factor(dist): # to control arc curvature
				if dist < 1000:
					return 0.05
				elif dist < 2500:
					return 0.01
				elif dist < 5000:
					return 0.005
				elif dist < 10000:
					return 0.001
				else:
					return 0.00001
			arc_height = dist * set_arc_height_factor(dist)  # Reduce height
			# Set linewidth according to the number of reads
			linewidth_factor = (2 - 1) / (junc_reads_max - junc_reads_min) if junc_reads_max != junc_reads_min else 1 # Scale linewidth from 1 to 2
			arc_linewidth = 1 + (junc_reads - junc_reads_min) * linewidth_factor
			# Create an Arc patch
			arc = Arc(
				(mx, my),
				dist,
				arc_height,
				angle=angle_deg,
				theta1=0,
				theta2=180,
				linewidth=arc_linewidth,
				edgecolor=color
			)
			ax.add_patch(arc)
			# Add junc_reads as text on the arc
			ax.text(
				mx, my + arc_height / 2, str(junc_reads),
				fontsize=8, ha='center', va='center', color='black',
				backgroundcolor='white', bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0')
			)
		ax.set_xlim(start, end)
		ax.set_ylim(bottom = 0, top = max(cov) * 1.4)
		ax.set_ylabel("Coverage", fontsize=10)
		if i < n_samples - 1:
			ax.set_xticklabels([])
		else:
			ax.set_xlabel(f"Genomic coordinate ({chrom})", fontsize=10)
		# Despine top, right, and bottom
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
	# Put title on the top subplot
	if pos_id:
		fig.suptitle(f"{chrom}:{start}-{end} ({strand})\n{pos_id}", fontsize=12, y=1.4 - 0.1 * n_samples)
	else:
		fig.suptitle(f"{chrom}:{start}-{end} ({strand})", fontsize=12, y=1.4 - 0.1 * n_samples)
	# Save plot
	plt.savefig(output, dpi=800, bbox_inches="tight")
