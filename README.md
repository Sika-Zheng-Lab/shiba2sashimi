# 🐕 shiba2sashimi 🍣 (v0.1.4)

[![GitHub License](https://img.shields.io/github/license/Sika-Zheng-Lab/shiba2sashimi)](https://github.com/Sika-Zheng-Lab/shiba2sashimi/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/947608002.svg)](https://doi.org/10.5281/zenodo.15042265)
[![GitHub Release](https://img.shields.io/github/v/release/Sika-Zheng-Lab/shiba2sashimi?style=flat)](https://github.com/Sika-Zheng-Lab/shiba2sashimi/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/Sika-Zheng-Lab/shiba2sashimi)](https://github.com/Sika-Zheng-Lab/shiba2sashimi/releases)
[![Create Release](https://github.com/Sika-Zheng-Lab/shiba2sashimi/actions/workflows/release.yaml/badge.svg)](https://github.com/Sika-Zheng-Lab/shiba2sashimi/actions/workflows/release.yaml)
[![Publish PyPI](https://github.com/Sika-Zheng-Lab/shiba2sashimi/actions/workflows/publish.yaml/badge.svg)](https://github.com/Sika-Zheng-Lab/shiba2sashimi/actions/workflows/publish.yaml)
[![Python](https://img.shields.io/pypi/pyversions/shiba2sashimi.svg?label=Python&color=blue)](https://pypi.org/project/shiba2sashimi/)
[![PyPI](https://img.shields.io/pypi/v/shiba2sashimi.svg?label=PyPI&color=orange)](https://pypi.org/project/shiba2sashimi/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/shiba2sashimi.svg?label=PyPI%20-%20Downloads&color=orange)](https://pypi.org/project/shiba2sashimi/)
[![Conda](https://img.shields.io/conda/v/bioconda/shiba2sashimi?color=3EB049)](https://anaconda.org/bioconda/shiba2sashimi)
[![Conda - Downloads](https://img.shields.io/conda/dn/bioconda/shiba2sashimi?label=Conda%20-%20Downloads&color=3EB049)](https://anaconda.org/bioconda/shiba2sashimi)
[![Docker](https://img.shields.io/docker/v/naotokubota/shiba2sashimi?color=blue&label=Docker)](https://hub.docker.com/r/naotokubota/shiba2sashimi)
[![Docker Pulls](https://img.shields.io/docker/pulls/naotokubota/shiba2sashimi)](https://hub.docker.com/r/naotokubota/shiba2sashimi)
[![Docker Image Size](https://img.shields.io/docker/image-size/naotokubota/shiba2sashimi)](https://hub.docker.com/r/naotokubota/shiba2sashimi)

A utility to create Sashimi plots, a publication-quality visualization of RNA-seq data, from [Shiba](https://github.com/Sika-Zheng-Lab/Shiba) output. Greatly inspired by [rmats2sashimiplot](https://github.com/Xinglab/rmats2sashimiplot) and [MISO](https://miso.readthedocs.io/en/fastmiso/sashimi.html)'s original implementation.

## Quick start

```bash
shiba2sashimi -e /path/to/Shiba/experiment_table.tsv \
-s /path/to/Shiba/output/ -o img/sashimi_example.png \
--id "SE@chr2@157561213-157561293@157560260-157561542"
```

![Sashimi plot example](https://raw.githubusercontent.com/Sika-Zheng-Lab/shiba2sashimi/main/img/sashimi_example.png)

## How to install

### pip

```bash
pip install shiba2sashimi
```

or

```bash
git clone https://github.com/Sika-Zheng-Lab/shiba2sashimi.git
cd shiba2sashimi
pip install .
```

### conda

```bash
conda create -n shiba2sashimi python=3.9
conda activate shiba2sashimi
conda install -c bioconda shiba2sashimi
```

### Docker

```bash
docker pull naotokubota/shiba2sashimi
```

## Dependencies

- python (>=3.9)
- numpy (>=1.18.0,<2.0.0)
- matplotlib (>=3.1.0)
- pysam (>=0.22.0)

## Usage

```bash
usage: shiba2sashimi [-h] -e EXPERIMENT -s SHIBA -o OUTPUT [--id ID] [-c COORDINATE] [--samples SAMPLES] [--groups GROUPS] [--colors COLORS] [--extend_up EXTEND_UP] [--extend_down EXTEND_DOWN]
                     [--smoothing_window_size SMOOTHING_WINDOW_SIZE] [--font_family FONT_FAMILY] [--dpi DPI] [-v]

shiba2sashimi v0.1.4 - Create Sashimi plot from Shiba output

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

## Contributing

Thank you for wanting to improve shiba2sashimi! If you have any bugs or questions, feel free to [open an issue](https://github.com/Sika-Zheng-Lab/shiba2sashimi/issues) or pull request.

## Citation

If you use shiba2sashimi in your research, please cite the original Shiba paper:

Kubota N, Chen L, Zheng S. [Shiba: a versatile computational method for systematic identification of differential RNA splicing across platforms](https://academic.oup.com/nar/article/53/4/gkaf098/8042001). *Nucleic Acids Research*  53(4), 2025, gkaf098.

## Authors

- Naoto Kubota ([0000-0003-0612-2300](https://orcid.org/0000-0003-0612-2300))
- Liang Chen ([0000-0001-6164-4553](https://orcid.org/0000-0001-6164-4553))
- Sika Zheng ([0000-0002-0573-4981](https://orcid.org/0000-0002-0573-4981))
