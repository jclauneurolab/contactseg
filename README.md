# Automatic Contact Segmentation

**Automatic Contact Segmentation** is a BIDS App for localizing stereoelectroencephalography (SEEG) contacts from post-operative CT scans. It uses the [nnUNetv2 framework](https://github.com/MIC-DKFZ/nnUNet) to train a 3D U-Net model for automatic SEEG contact segmentation. The app integrates [Snakemake](https://snakemake.readthedocs.io/) and [SnakeBIDS](https://github.com/akhanf/snakebids) for workflow management and reproducibility. The project is managed using [pixi](https://github.com/prefix-dev/pixi), which provides fast dependency resolution and environment management.

To install the app, first install `pixi`. On macOS and Linux, run:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

On Windows, run:
```
powershell -ExecutionPolicy Bypass -c "irm -useb https://pixi.sh/install.ps1 | iex"
```

Alternatively, you can install it from PyPI using:

```
pip install pixi
```

or with pipx:

```
pipx install pixi
```

Once pixi is installed, setup the environment:

```
pixi install
```

To verify that the installation was successful, run:

```
pixi run contactseg --help
```

This should display the help message for the CLI. 

To use the BIDS App, run:

```
pixi run contactseg /path/to/bids/dataset /path/to/output/derivatives participant --cores all
```

Replace /path/to/bids/dataset with the path to your BIDS-compliant input dataset and /path/to/output/derivatives with the desired output directory.

If you're developing the app and want to install development dependencies such as linters and formatters, run:

```
pixi install -e dev
```
You can then run quality checks and formatting with:

```
pixi run -e dev quality_check
pixi run -e dev quality_fix
```