# shiba2sashimi

## Usage

```bash
usage: shiba2sashimi.py [-h] -e EXPERIMENT -s SHIBA -o OUTPUT [--id ID] [-c COORDINATE] [--samples SAMPLES] [--groups GROUPS] [--colors COLORS] [--extend_up EXTEND_UP] [--extend_down EXTEND_DOWN] [--font_family FONT_FAMILY]
                        [--dpi DPI] [-v]

shiba2sashimi v0.1.0 - Create Sashimi plot from Shiba output

options:
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
                        Extend the plot upstream. Only used when not providing coordinates. Default: 1000
  --extend_down EXTEND_DOWN
                        Extend the plot downstream. Only used when not providing coordinates. Default: 1000
  --font_family FONT_FAMILY
                        Font family for labels
  --dpi DPI             DPI of the output figure. Default: 300
  -v, --verbose         Increase verbosity
```
