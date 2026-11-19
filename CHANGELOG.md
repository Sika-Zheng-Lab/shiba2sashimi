# Change log

All notable changes to this shiba2sashimi project will be documented in this file.

## [v0.2.0] - 2026-11-19

### Added

- RI exon–intron boundary counts shown as labeled dashed ticks at each boundary,
  with arithmetic group averaging and the same count threshold as junctions.
- `--group-average`: one mean coverage/junction track per selected group; mean
  nonmissing sample PSI and valid sample counts are labeled explicitly.
- Transcript structures beneath sashimi tracks, with strand arrows and target
  event highlighting; automatic assembled-GTF discovery, `--gtf`/GTF.gz,
  `--max-transcripts`, and `--no-transcripts` options. Event-local structures are
  used when no matching GTF models are available, without inventing exon boundaries.
- Automated unit and CLI integration tests using generated indexed BAM files;
  GitHub Actions for Python 3.9/3.12; real-data examples for all eight
  event types, each with individual-sample and group-mean plots.

### Fixed

- Recognize legacy `assembled.gtf` and StringTie reference-gene assignments,
  including novel transcripts without an explicit `ref_gene_id`.
- Parse experiment tables by column name, accepting Shiba v1.0.0 `technology`
  metadata and sample names beginning with `sample`.
- Handle empty/trailing and nonfinite PSI values without dropping sample columns.
- Derive event windows from all junction coordinates, including Shiba v1.0.0
  minus-strand MXE and semicolon-separated AFE/ALE/MSE paths.
- Align Shiba's flanking-exon junction IDs with zero-based BAM/GTF plotting
  coordinates; exclude RI boundary counts from splice-junction arcs.
- Handle empty junction regions and zero coverage; support standard BAI/CSI names,
  symmetric chromosome aliases, and useful input validation errors.

### Changed

- Space transcript strand arrows uniformly on a shared coordinate grid (24 pt
  on the page), with clearance from exon boundaries.
- Remove the separate target-event label and bar above transcript models; retain
  event highlighting within the transcript structures.
- Use fixed-point count/depth labels, retaining up to three decimal places for
  fractional group means without scientific notation.
- Plot raw coverage by default (`--smoothing_window_size 1`); explicit smoothing
  is applied after group averaging. Use `21` to request the previous window.
- Share depth scales across displayed tracks; retain fractional group counts and
  apply `--minimum_junc_reads` to the group mean when averaging.
- Use installed Arial by default, a restrained palette, editable SVG text and
  embedded TrueType PDF fonts. An unavailable font produces an actionable error.
- `--coordinate` is explicitly zero-based/end-exclusive. It must overlap the
  selected event when combined with `--id`.

## [v0.1.7] - 2025-08-06

### Added

- Added `--minimum_junc_reads` option to allow users to specify the minimum number of reads required to draw a junction arc.

### Changed

- Changed internal logic to handle coordinates more flexibly to support future releases of Shiba.

## [v0.1.6] - 2025-06-20

### Fixed

- Fixed a bug that caused an error when the *csi* index file for the BAM file exists instead of the *bai* index file.

### Added

- Added `--width` option to allow users to specify the width of the plot.
- Added `--nolabel` and `--nojunc` options to allow users to disable sample labels and junction arcs, respectively.

### Changed

- Adjusted font size of labels to improve readability.

## [v0.1.5] - 2025-06-02

### Fixed

- Corrected chromosome prefix handling in read coverage calculation to avoid KeyError when chromosome names do not match the expected format.

## [v0.1.4] - 2025-04-13

### Changed

- Improved x-axis labeling to avoid overlapping with the plot area.

## [v0.1.3] - 2025-03-25

### Changed

- Changed a way to draw junction arcs to avoid drawing arcs outside the plot area; use bezier curves instead of circular arcs.

## [v0.1.2] - 2025-03-19

### Fixed

- Fixed a bug that caused bioconda build to fail.

## [v0.1.1] - 2025-03-19

### Added

- Dockerfile: Added Arial font to the image.

### Fixed

- Fixed title position in the plot.
- Fixed a logic to calculate arc height.
- Fixed other minor issues.

### Changed

- Dependencies: Updated `python` to `>3.9`.
- Dockerfile: Changed the base image to `python:3.9-bullseye` from `python:3.7-buster`.

## [v0.1.0] - 2025-03-17

Initial release of 🐕 shiba2sashimi 🍣
