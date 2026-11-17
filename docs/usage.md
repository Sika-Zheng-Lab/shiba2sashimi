# shiba2sashimi usage

## Inputs and compatibility

Use a **bulk RNA-seq** Shiba/SnakeShiba working directory and its experiment table.
Both legacy outputs and Shiba v1.0.0 Fisher/beta-binomial outputs are supported;
PSI labels always come from the observed `PSI_matrix_sample.txt`, not regression
predictions. scShiba's barcode experiment tables and cell-filtered BAM coverage
are not supported.

The experiment table has named `sample`, `bam`, and `group` columns. Extra columns
(such as Shiba v1.0.0 `technology`) are accepted, and column order is flexible:

```tsv
sample	bam	group	technology
control_1	/path/to/control_1.bam	control	short
control_2	/path/to/control_2.bam	control	short
treatment_1	/path/to/treatment_1.bam	treatment	short
```

BAM files must be coordinate-sorted and indexed with BAI or CSI. Paths are used
as supplied, relative to the current working directory when not absolute. The
working directory supplies:

```text
workdir/
├── annotation/assembled_annotation.gtf  # optional
├── events/EVENT_<TYPE>.txt              # required with --id
├── junctions/junctions.bed              # required unless --nojunc
└── results/splicing/PSI_matrix_sample.txt  # required with --id
```

Use the merged `junctions.bed` written by Shiba, including its Tosa-based v1.0.0
output. Raw Tosa output is not a substitute. For an RI event selected with `--id`,
exon–intron boundary counts appear below the coverage baseline as dashed ticks
at the two boundary coordinates, labeled `EI: count`. These are unspliced
boundary counts, distinct from the spliced junction arc. Group plots use their
arithmetic means, including zeros; `--minimum_junc_reads` and `--nojunc` apply
to both junctions and boundaries. Absent boundary rows are not drawn, and
coordinate-only plots continue to show splice junctions only. Count labels use
ordinary decimal notation with up to three decimal places, never scientific
notation. An empty junction region is valid; missing
selected sample columns are an error rather than an assumed zero.

Shiba v1.0.0 defines MXE exon A/B in transcript orientation, including minus-strand
events. shiba2sashimi follows the event table's exon/junction assignments and uses
the extrema of all junction coordinates to obtain the plotting window. Legacy
MXE tables remain readable; their PSI meaning follows that original Shiba run.
See [Shiba's release notes](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/CHANGELOG.md)
and [event definitions](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/docs/output/shiba.md).

## Sample and group tracks

```bash
# Default: one track per sample
shiba2sashimi -e experiment.tsv -s workdir -o samples.pdf \
  --id 'SE@chr2@157561213-157561293@157560260-157561542'

# One averaged track per group, including more than two groups
shiba2sashimi -e experiment.tsv -s workdir -o groups.svg \
  --id 'SE@chr2@157561213-157561293@157560260-157561542' \
  --groups control,treatment1,treatment2 --group-average

# Subset replicates before grouping; without --groups, --samples determines membership
shiba2sashimi -e experiment.tsv -s workdir -o subset.pdf \
  -c chr2:157560000-157562000 --samples control_1,treatment_1 --group-average
```

`--groups` selects/orders groups and takes precedence over `--samples`. It does
not itself average samples. Without either selector, all experiment-table samples
are included in table order; group averages follow first occurrence of each group.
Selected groups are not limited to Shiba's two-group differential comparison.

For `--group-average`:

| Quantity | Calculation |
|---|---|
| Depth | Arithmetic mean of raw per-base BAM depth across all selected replicates in the group, including zeros |
| Junction count | Arithmetic mean of counts across all selected replicates; absent junction IDs count as zero |
| PSI | Arithmetic mean of finite, nonmissing sample PSI, expressed as a percentage |
| Group `n` | Number of selected replicates contributing to depth/count means |
| PSI `n` | Number of valid sample PSI values; an all-missing group shows `NA`, `n = 0` |

These are unweighted sample means. They are **not** Shiba's pooled-count group PSI
or beta-binomial model predictions. No library-size normalization is applied.
Coverage uses `pysam.count_coverage` defaults (A/C/G/T bases, base quality ≥15,
skipping unmapped, secondary, QC-failed, and duplicate reads). Coverage and Shiba
junction counts can therefore use different read filters; depth is not recomputed
from junction counts.

All tracks share the same depth scale. Junction widths scale with displayed
positive counts. `--minimum_junc_reads` is applied **after averaging**; fractional
thresholds are accepted, and zero-count arcs are never drawn. Fractional labels
use up to three significant digits; integer counts are shown in full.

## Transcript structures and event annotation

```bash
shiba2sashimi -e experiment.tsv -s workdir -o event.pdf \
  --id 'SE@chr2@157561213-157561293@157560260-157561542' \
  --group-average --gtf reference.gtf.gz
```

Transcript models appear below the read tracks, with exon blocks and directional
intron arrows. Red blocks and a labeled target-event marker identify the event's
exons; for RI, the retained intron is highlighted instead. All tracks use the same
genomic x-axis. Features outside the requested window are clipped.

By default, `annotation/assembled_annotation.gtf` is used when present, with
`annotation/assembled.gtf` accepted for older Shiba runs. For
reference-only Shiba runs, supply the reference annotation with `--gtf`. Plain
GTF and gzipped GTF are supported; `gene_id` and `transcript_id` attributes identify
models. With `--id`, models are restricted to the event's gene ID. StringTie
`ref_gene_id` attributes are supported. Novel transcripts without this attribute
inherit the most common reference-gene assignment within their assembled gene,
matching Shiba's mapping; explicit reference assignments are preserved.

Up to 20 transcripts are displayed by default, prioritizing those with the most
matching event junctions and then transcript ID. A warning reports truncation;
use `--max-transcripts` to increase the limit. Models show annotation structure,
not evidence that every transcript is expressed.

When there are no matching GTF models, `--id` still produces two **local event
structures** from the event table. Only recorded exons and introns are shown;
unknown flanking exon extents are not invented. Local paths are labeled `(local)`
and are not full-length transcript assignments. All eight Shiba event types are
supported. Coordinate-only plots can show GTF models but have no target-event
highlight. Use `--no-transcripts` to hide the lower annotation panel.

## Coordinates, smoothing, and export

`--coordinate chr:start-end` uses **zero-based, end-exclusive** coordinates,
matching BAM queries. An explicit coordinate overrides the default event window;
when combined with `--id`, it must overlap that event on the same chromosome.
Event windows include all event junctions and extend 500 bases upstream/downstream
by default. `--extend_up` and `--extend_down` follow transcript orientation, and
negative starts are clipped to zero. A region beyond the BAM chromosome length is
reported as an error; use `--coordinate` to select a valid window.

Shiba event exons and GTF exons are one-based inclusive. Shiba junction IDs refer
to the flanking exon bases. These are converted internally to align coverage,
exons, and splice boundaries. Event IDs shown in the title retain their original
Shiba notation. Exact chromosome names are preferred in BAMs, with an unambiguous
`chr` alias accepted when necessary.

Coverage is **unsmoothed by default**. To request median smoothing:

```bash
shiba2sashimi -e experiment.tsv -s workdir -o smoothed.pdf \
  -c chr2:157560000-157562000 --group-average --smoothing_window_size 21
```

The filter is applied after averaging; even window sizes round up to the next odd
number. It does not alter junction counts or PSI. Use `1` to disable it.

Choose `.svg` or `.pdf` for editable vector output, or `.png` for a raster preview
(default 300 dpi). Arial is the default and must be installed; the program does
not silently substitute a font. To explicitly choose another available font,
use `--font_family 'DejaVu Sans'`. `--width` accepts inches, including fractional
values; `--colors` sets colors in selected group order. `--nolabel` hides sample
labels/PSI and `--nojunc` hides junction arcs/counts.

Run `shiba2sashimi --help` for all options and `shiba2sashimi --version` for the
installed version.

## Tests and examples

```bash
python -m pip install -e '.[test]'
python -m pytest --cov=shiba2sashimi --cov-report=term-missing
```

The [Nov24 example gallery](https://github.com/Sika-Zheng-Lab/shiba2sashimi/tree/develop/img/examples)
contains all eight event types from Control/Upf2_cKO bulk RNA-seq results, with
sample and group plots. Example-generation scripts and source BAMs remain local
and are not distributed. Automated tests generate their own synthetic data and
do not require samtools, Shiba/Tosa installations, or external sequencing datasets.
