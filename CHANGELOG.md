# Change log

All notable changes to this shiba2sashimi project will be documented in this file.

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
