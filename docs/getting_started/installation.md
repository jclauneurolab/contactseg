# Installation

**Automatic Contact Segmentation** is a BIDS App for localizing stereoelectroencephalography (SEEG) contacts from post-operative CT scans. It uses the [nnUNetv2 framework](https://github.com/MIC-DKFZ/nnUNet) to train a 3D U-Net model for automatic SEEG contact segmentation. The app integrates [Snakemake](https://snakemake.readthedocs.io/) and [SnakeBIDS](https://github.com/akhanf/snakebids) for workflow management and reproducibility. The project is managed using [uv](https://github.com/astral-sh/uv), which provides fast dependency resolution and environment management.

To install the app, first install `uv`. On macOS and Linux, run:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

On Windows, run:
```
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Alternatively, you can install it from PyPI using:

```
pip install uv
```

or with pipx:

```
pipx install uv
```

Once uv is installed, create and activate a virtual environment:

```
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

Then install the package in editable mode:

```
uv pip install -e .
```

To verify that the installation was successful, run:

```
contactseg -h
```

This should display the help message for the CLI. If it doesn’t, ensure that your virtual environment is activated and that the installation completed without errors.

To use the BIDS App, run:

```
contactseg /path/to/bids/dataset /path/to/output/derivatives participant --cores all
```

Replace /path/to/bids/dataset with the path to your BIDS-compliant input dataset and /path/to/output/derivatives with the desired output directory.

If you're developing the app and want to install development dependencies such as linters and formatters, run:

```
uv pip install -e .[dev]
```


# Running the app

Do a dry-run first (`-n`) and simply print (`-p`) what would be run:

```bash
contactseg /path/to/bids/dir /path/to/output/dir participant -np
```

Run the app, using all cores::

```bash
contactseg /path/to/bids/dir /path/to/output/dir participant --cores all
```

If any workflow rules require containers, then run with the `--use-singularity` option.

# Generating a report

After your processing is complete, you can use snakemake's `--report` feature to generate
an HTML report. This report will include a graph of all the jobs run, with clickable nodes
to inspect the shell command or python code used in each job, along with the config files and
run times for each job. Workflows may also contain append images for quality assurance or to
summarize outputs, by using the `report(...)` function on any snakemake output.

To generate a report, run:

```bash
contactseg /path/to/bids/dir /path/to/output/dir participant --report
```

# Compute Canada Instructions

## Setting up a dev environment

Here are some instructions to get your python environment set-up on graham to run contactseg:

### Create a virtualenv and activate it:

```bash
cd $SCRATCH
module load python/3
virtualenv venv_contactseg
source venv_contactseg/bin/activate
```

### Follow the steps above to install from github repository

## Install job submission helpers

Snakemake can submit jobs with SLURM, but you need to set-up a Snakemake profile to enable this. The Khan lab has a
snakemake profile that is configured for graham but is customizable upon install, please see `cc-slurm <https://github.com/khanlab/cc-slurm>` for more info.

If you don't need Snakemake to parallelize jobs across different nodes, you can make use of the simple job submission wrappers in `neuroglia-helpers <https://github.com/khanlab/neuroglia-helpers>`, e.g. `regularSubmit` or `regularInteractive` wrappers.

These are used in the instructions below.

## Running jobs on Compute Canada

In an interactive job (for testing):

```bash
regularInteractive -n 8 contactseg bids_dir out_dir participant --participant_label 001 -j 8
```

Submitting a job (for larger cores, more subjects), still single job, but snakemake will parallelize over the 32 cores:

```bash
regularSubmit -j Fat contactseg bids_dir out_dir participant  -j 32
```

Scaling up to ~hundred subjects (needs cc-slurm snakemake profile installed), submits 1 16core job per subject:

```bash
contactseg bids_dir out_dir participant  --profile cc-slurm
```

Scaling up to even more subjects (uses group-components to bundle multiple subjects in each job), 1 32core job for N subjects (e.g. 10):

```bash
contactseg bids_dir out_dir participant  --profile cc-slurm --group-components subj=10
```
