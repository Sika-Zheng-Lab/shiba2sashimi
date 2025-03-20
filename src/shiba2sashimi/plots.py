import sys
import os
import numpy as np
import logging
# Configure logging
logger = logging.getLogger(__name__)
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

def sashimi(coverage_dict, junctions_dict, experiment_dict, samples, groups, colors, chrom, start, end, output, pos_id = None, coordinate = None, strand = None, gene_name = None, junction_direction_dict = None, psi_values_dict = None, font_family = None, dpi = 300):
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
	gs = fig.add_gridspec(n_samples, 1, hspace=0.8)
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
			# Calculate midpoint (to use as the center of the arc)
			mx = (x1 + x2) / 2
			my = (y1 + y2) / 2
			# Calculate distance between the two points
			dx = x2 - x1
			dy = y2 - y1
			dist = np.hypot(dx, dy)
			# Angle (in degrees) of the line between the two points
			angle_deg = np.degrees(np.arctan2(dy, dx))
			# Set arc height according to the coverage
			arc_height = cov_max/3 if direction == "up" else -cov_max/2
			# Set linewidth according to the number of reads
			linewidth_factor = (2 - 1) / (junc_reads_max - junc_reads_min) if junc_reads_max != junc_reads_min else 1 # Scale linewidth from 1 to 2
			arc_linewidth = 1 + (junc_reads - junc_reads_min) * linewidth_factor
			if junc_reads == 0:
				arc_linewidth = 0.5
			# Create an Arc patch
			arc = Arc(
				(mx, my),
				dist,
				abs(arc_height),
				angle=angle_deg,
				theta1=0 if direction == "up" else 180,
				theta2=180 if direction == "up" else 360,
				linewidth=arc_linewidth,
				edgecolor=color,
				clip_on=False  # Allow drawing outside the axes
			)
			ax.add_patch(arc)
			text_y_offset = arc_height / 2
			# Add junc_reads as text on the arc
			ax.text(
				mx, my + text_y_offset, str(junc_reads),
				fontsize=8, ha='center', va='center', color='black',
				backgroundcolor='white', bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0'),
				clip_on=False  # Allow text to be drawn outside the axes
			)
		ax.set_xlim(start, end)
		ax.set_ylim(bottom = 0, top = max(cov) * 1.4)
		ax.set_ylabel("Coverage", fontsize=8)
		if i < n_samples - 1:
			ax.set_xticklabels([])
		else:
			ax.set_xlabel(f"Genomic coordinate ({chrom})", fontsize=10)
		# Despine top, right, and bottom
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		# Remove xticks for all samples
		ax.set_xticks([])
	# Add xticks to the bottom subplot
	ax.set_xticks(np.arange(start, end, step=(end - start) // 10))
	ax.set_xticklabels(np.arange(start, end, step=(end - start) // 10), rotation=45, ha='right', fontsize=6)
	# Put title on the top subplot
	title = ""
	if coordinate:
		title += f"{chrom}:{start}-{end}"
	if pos_id:
		title += f"\n{pos_id}, {gene_name} ({strand})"
	if title:
		ax_top = fig.axes[0]  # Get the top subplot
		bbox = ax_top
		bbox = ax_top.get_position()
		x_center = (bbox.x0 + bbox.x1) / 2
		y_top = bbox.y1
		plt.gcf().text(
			x_center, y_top + 0.1,
			title, fontsize=12, ha='center', va='bottom'
		)
	# Save plot
	plt.savefig(output, dpi=dpi, bbox_inches="tight")
