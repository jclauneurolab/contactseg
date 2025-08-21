# Command line interface


## contactseg Command-line interface

 The following can also be seen by entering ``contactseg -h`` into your terminal. 

These are all the required and optional arguments contactseg accepts in order to run flexibly on many different input data types and with many options, but in most cases only the required arguments are needed. 

<!-- the below code needs to be debugged -->
```{argparse}
---
filename: ../contactseg/run.py
func: get_parser
prog: contactseg
---
```


## Snakemake command-line interface

In addition to the above command-line arguments, Snakemake arguments are also be passed at the `contactseg` command-line. 

The most critical of these is the `--cores` or `-c` argument, which is a **required** argument for contactseg. 

The complete list of [Snakemake](https://snakemake.readthedocs.io/en/stable/) arguments are below, and mostly act to determine your environment and App behaviours. They will likely only need to be used for running in cloud environments or troubleshooting. These can be listed from the command-line with `contactseg --help-snakemake`.  

```{argparse}
---
module: snakemake.cli
func: get_argument_parser
prog: snakemake
---
```