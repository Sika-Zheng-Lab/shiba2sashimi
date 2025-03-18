# 🐕 shiba2sashimi 🍣 (v0.1.0)

[![GitHub License](https://img.shields.io/github/license/Sika-Zheng-Lab/shiba2sashimi)](https://github.com/Sika-Zheng-Lab/shiba2sashimi/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/947608002.svg)](https://doi.org/10.5281/zenodo.15042265)
[![GitHub Release](https://img.shields.io/github/v/release/Sika-Zheng-Lab/shiba2sashimi?style=flat)](https://github.com/Sika-Zheng-Lab/shiba2sashimi/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/Sika-Zheng-Lab/shiba2sashimi)](https://github.com/Sika-Zheng-Lab/shiba2sashimi/releases)
[![Create Release](https://github.com/Sika-Zheng-Lab/shiba2sashimi/actions/workflows/release.yaml/badge.svg)](https://github.com/Sika-Zheng-Lab/shiba2sashimi/actions/workflows/release.yaml)
[![Publish PyPI](https://github.com/Sika-Zheng-Lab/shiba2sashimi/actions/workflows/publish.yaml/badge.svg)](https://github.com/Sika-Zheng-Lab/shiba2sashimi/actions/workflows/publish.yaml)
[![PyPI](https://img.shields.io/pypi/v/shiba2sashimi)](https://pypi.org/project/shiba2sashimi/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/shiba2sashimi)](https://pypi.org/project/shiba2sashimi/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/shiba2sashimi)](https://pypi.org/project/shiba2sashimi/)
[![Docker Pulls](https://img.shields.io/docker/pulls/naotokubota/shiba2sashimi)](https://hub.docker.com/r/naotokubota/shiba2sashimi)
[![Docker Image Size](https://img.shields.io/docker/image-size/naotokubota/shiba2sashimi)](https://hub.docker.com/r/naotokubota/shiba2sashimi)

A utility to create Sashimi plots, a publication-quality visualization of RNA-seq data, from [Shiba](https://github.com/Sika-Zheng-Lab/Shiba) output. Greatly inspired by [rmats2sashimiplot](https://github.com/Xinglab/rmats2sashimiplot) and [MISO](https://miso.readthedocs.io/en/fastmiso/sashimi.html)'s original implementation.

## Quick start

```bash
shiba2sashimi -e /path/to/Shiba/experiment_table.tsv \
-s /path/to/Shiba/output/ -o img/sashimi_example.png \
--id "SE@chr11@5091460-5091477@5088146-5091859"
```

![Sashimi plot example](img/sashimi_example.png)

## How to install

```bash
pip install shiba2sashimi
```

or

```bash
git clone https://github.com/Sika-Zheng-Lab/shiba2sashimi.git
cd shiba2sashimi
pip install .
```

You can also run shiba2sashimi without pip installing by running the script directly if you have installed all dependencies.

```bash
git clone https://github.com/Sika-Zheng-Lab/shiba2sashimi.git
cd shiba2sashimi
python -m shiba2sashimi.main
```

## Dependencies

- numpy (>=1.18.0,<2.0.0)
- matplotlib (>=3.1.0)
- pysam (>=0.22.0)

## Usage

```bash
usage: shiba2sashimi [-h] -e EXPERIMENT -s SHIBA -o OUTPUT [--id ID] [-c COORDINATE] [--samples SAMPLES] [--groups GROUPS] [--colors COLORS] [--extend_up EXTEND_UP] [--extend_down EXTEND_DOWN]
                     [--smoothing_window_size SMOOTHING_WINDOW_SIZE] [--font_family FONT_FAMILY] [--dpi DPI] [-v]

shiba2sashimi v0.1.0 - Create Sashimi plot from Shiba output

optional arguments:
  -h, --help            show this help message and exit
  -e EXPERIMENT, --experiment EXPERIMENT
                        Experiment table used for Shiba
  -s SHIBA, --shiba SHIBA
                        Shiba output directory
  -o OUTPUT, --output OUTPUT
                        Output file
  --id ID               Positional ID (pos_id) of the event to plot
  -c COORDINATE, --coordinate COORDINATE
                        Coordinates of the region to plot
  --samples SAMPLES     Samples to plot. e.g. sample1,sample2,sample3 Default: all samples in the experiment table
  --groups GROUPS       Groups to plot. e.g. group1,group2,group3 Default: all groups in the experiment table. Overrides --samples
  --colors COLORS       Colors for each group. e.g. red,orange,blue
  --extend_up EXTEND_UP
                        Extend the plot upstream. Only used when not providing coordinates. Default: 500
  --extend_down EXTEND_DOWN
                        Extend the plot downstream. Only used when not providing coordinates. Default: 500
  --smoothing_window_size SMOOTHING_WINDOW_SIZE
                        Window size for median filter to smooth coverage plot. Greater value gives smoother plot. Default: 21
  --font_family FONT_FAMILY
                        Font family for labels
  --dpi DPI             DPI of the output figure. Default: 300
  -v, --verbose         Increase verbosity
```

## Authors

- Naoto Kubota ([0000-0003-0612-2300](https://orcid.org/0000-0003-0612-2300))
- Liang Chen ([0000-0001-6164-4553](https://orcid.org/0000-0001-6164-4553))
- Sika Zheng ([0000-0002-0573-4981](https://orcid.org/0000-0002-0573-4981))
