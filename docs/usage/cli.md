# Command line interface


## contactseg Command-line interface

 The following can also be seen by entering ``contactseg -h`` into your terminal. 

These are all the required and optional arguments contactseg accepts in order to run flexibly on many different input data types and with many options, but in most cases only the required arguments are needed. 

<!-- the below code needs to be debugged -->
```{argparse}
---
module: contactseg.run
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


# Common Usage Scenarios

## Using pre-computed sMRIPrep derivatives

`contactseg` needs to register your native images to standard MNI space for certain features. If you have already processed your dataset using sMRIPrep, you can pass the output directory to contactseg. The pipeline will automatically skip computing its own nonlinear registrations and reuse the standard `.h5` transforms. Using sMRIPrep derivatives also adds a column to the spreadsheet that states whether the contacts are in white matter, grey matter or cerebral spinal fluid:

```bash
contactseg /path/to/bids/dataset /path/to/output/derivatives participant --cores all --SMRIPREP-DIR /path/to/smriprep/derivatives
```

## Automatic Atlas Labeling

You can automatically assign anatomical regions (e.g., Amygdala, Hippocampus) to the segmented SEEG contacts using standard atlases (like the CerebrA atlas).

**Labeling in MNI Space (Default):** By default, the pipeline transforms the native contact coordinates into MNI space to extract standard atlas labels.

```bash
contactseg /path/to/bids/dataset /path/to/output/derivatives participant --cores all --atlas_labels
```

**Labeling strictly in Native Space:** If you prefer to warp the standard atlas back into the patient's native anatomy to extract labels, add the `--use_native_space` flag:

```bash
contactseg /path/to/bids/dataset /path/to/output/derivatives participant --cores all --atlas_labels --use_native_space
```

`contactseg` uses ANTsPy SyN to non-linearly register the native T1w image to the template space. The resulting transformations, such as forward warp, inverse warp, and affine matrices, are used to either transform the native points to the MNI space or the MNI atlas to native space. The `lookup_atlas_labels` function uses a fuzzy search to locate the nearest voxel within a specified radius that has a non-zero voxel label. That numerical label is then matched to the atlas TSV file to retrieve the anatomical label name along with whether it is in the left or right hemisphere. The `--atlas_labels` flag outputs the `atlas_labelled_contactseg.fcsv` file, which provides the coordinates in native space with their respective anatomical labels, as well as the `mni_transformed_contactseg.fcsv` file and a spreadsheet containing the native coordinates, MNI coordinates, anatomical label, numbered label, and tissue type.


# BIDS Compliant Inputs

Your dataset must follow the Brain Imaging Data Structure (BIDS). Below is a specific example of a compliant directory structure for a subject with pre- and post-operative sessions.

**Example: `sub-P167`**

```text
sub-P167/
├── ses-post
│   └── ct
│       ├── sub-P167_ses-post_acq-Electrode_run-01_ct.json
│       └── sub-P167_ses-post_acq-Electrode_run-01_ct.nii.gz
├── ses-pre
│   ├── anat
│   │   ├── sub-P167_ses-pre_run-01_FLAIR.json
│   │   ├── sub-P167_ses-pre_run-01_FLAIR.nii.gz
│   │   ├── sub-P167_ses-pre_run-01_T1w.json
│   │   ├── sub-P167_ses-pre_run-01_T1w.nii.gz
│   │   ├── sub-P167_ses-pre_run-02_T1w.json
│   │   └── sub-P167_ses-pre_run-02_T1w.nii.gz
│   ├── ct
│   │   └── sub-P167_ses-pre_acq-Electrode_run-01_ct.nii.gz
│   ├── pet
│   │   ├── sub-P167_ses-pre_task-rest_run-01_pet.json
│   │   └── sub-P167_ses-pre_task-rest_run-01_pet.nii.gz
│   └── sub-P167_ses-pre_scans.tsv
├── sub-P167_run-01_planned.fcsv
└── sub-P167_run-01_xfm.txt
```
