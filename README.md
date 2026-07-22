# Automatic Contact Segmentation

**Automatic Contact Segmentation** is a BIDS App for localizing stereoelectroencephalography (SEEG) contacts from post-operative CT scans. It uses the [nnUNetv2 framework](https://github.com/MIC-DKFZ/nnUNet) to train a 3D U-Net model for automatic SEEG contact segmentation. The app integrates [Snakemake](https://snakemake.readthedocs.io/) and [SnakeBIDS](https://github.com/akhanf/snakebids) for workflow management and reproducibility. The project is managed using [pixi](https://github.com/prefix-dev/pixi), which provides fast dependency resolution and environment management.

---

## Installation & Setup

### 1. Install Pixi
If you do not have `pixi` installed, run:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

> **Note:** If `pixi` is already installed but yields a lock file version mismatch error, update it by running `pixi self-update`.

### 2. Set Up Environment & Git LFS
This repository uses **Git LFS** to manage sample datasets and large binary assets. After cloning the repository, initialize LFS via `pixi` to fetch the complete NIfTI datasets rather than 134-byte pointer files:

```bash
# Install environment dependencies
pixi install

# Initialize Git LFS and pull sample datasets
pixi run git-lfs install --local
pixi run git-lfs pull
```

---

## Verification

To verify that the installation was successful, run:

```bash
pixi run contactseg -h
```

If you see a help message listing all available command-line options, you’re ready to use `contactseg`!

---

## Usage

To execute the BIDS App on a dataset:

```bash
pixi run contactseg /path/to/bids/dataset /path/to/output/derivatives participant --cores all
```

### Running with GPU Acceleration (Recommended)
When executing on a CUDA-enabled GPU node:

```bash
pixi run contactseg /path/to/bids/dataset /path/to/output/derivatives participant --transform --label --cores all --use_gpu
```

---

## Development

If you are developing the app and need development dependencies (linters and formatters):

```bash
# Install development environment
pixi install -e dev

# Run quality checks and formatting
pixi run -e dev quality_check
pixi run -e dev quality_fix
```
