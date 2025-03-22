# shiba2sashimi usage

All you need is installing `shiba2sashimi` and running it with the appropriate inputs. Here is a quick guide to help you get started.

## Prepare inputs

You guys should have done the Shiba analysis, right? Then you don't need to worry about the inputs for `shiba2sashimi`. You are ready to go! But just in case, make sure you have the following inputs ready:

- **Experiment table**: A tab-separated file containing the sample information used for Shiba. The first row should be the header. The first column should be the sample name, the second column should be path to the BAM file, and the third column should be the group name. For example:

``` tsv
sample  bam group
sample1 /path/to/sample1.bam  Ref
sample2 /path/to/sample2.bam  Ref
sample3 /path/to/sample3.bam  Alt
sample4 /path/to/sample4.bam  Alt
```

!!! Important

    Please make sure you can access the BAM files and index files (.bai) from the path provided in the experiment table as `shiba2sashimi` will read the BAM files to create the Sashimi plot.

- **Shiba output directory**: The directory containing the output files from Shiba. The typical structure of the Shiba output directory is as follows:

```bash
path/to/Shiba/output/
├── annotation
├── events
├── junctions
├── plots
├── report.txt
└── results
```

You are ready to create a Sashimi plot! Please refer to the following sections for more details.

- [Basic usage](examples/basic_usage.md)
- [Customizing the plot](examples/customizing_the_plot.md)


## All options of `shiba2sashimi`

``` bash
usage: shiba2sashimi [-h] -e EXPERIMENT -s SHIBA -o OUTPUT [--id ID] [-c COORDINATE] [--samples SAMPLES] [--groups GROUPS] [--colors COLORS] [--extend_up EXTEND_UP] [--extend_down EXTEND_DOWN]
                     [--smoothing_window_size SMOOTHING_WINDOW_SIZE] [--font_family FONT_FAMILY] [--dpi DPI] [-v]

shiba2sashimi v0.1.2 - Create Sashimi plot from Shiba output

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