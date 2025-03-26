# 🐕 shiba2sashimi 🍣 (v0.1.2)

Here is the documentation page of [shiba2sashimi](https://github.com/Sika-Zheng-Lab/shiba2sashimi), a utility to create Sashimi plots from [Shiba](https://github.com/Sika-Zheng-Lab/Shiba) output. Greatly inspired by [rmats2sashimiplot](https://github.com/Xinglab/rmats2sashimiplot) and [MISO](https://miso.readthedocs.io/en/fastmiso/sashimi.html)'s original implementation.

## Quick start

```bash
shiba2sashimi -e /path/to/Shiba/experiment_table.tsv \
-s /path/to/Shiba/output/ -o img/sashimi_example.png \
--id "SE@chr2@157561213-157561293@157560260-157561542"
```

![Sashimi plot example](https://raw.githubusercontent.com/Sika-Zheng-Lab/shiba2sashimi/main/img/sashimi_example.png){ width=80% }

## Contents

- [Installation](installation.md)
- [Usage](usage.md)
- Examples
    - [Basic usage](examples/basic_usage.md)
    - [Customizing the plot](examples/customizing_the_plot.md)