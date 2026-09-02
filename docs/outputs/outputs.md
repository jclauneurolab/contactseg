# Output Directory Structure

In the specified `/path/to/output/dir`, the pipeline generates a standardized root structure containing metadata and subject derivatives:

```text
/path/to/output/dir/
├── config/        # Records of code versions, parameters, and input paths
├── logs/          # Detailed execution logs for each processing step
├── sub-{subject}/ # Primary subject-level derivatives
├── .snakebids/    # Hidden metadata for Snakebids tracking
└── .snakemake/    # Hidden metadata for Snakemake state and checkpoints
```

## The `sub-{subject}` Directory

The main outputs are organized according to BIDS spatial and datatype conventions. The contents of this folder expand based on the optional flags provided during execution.

### 1. Default Output

By default (without any optional flags), the pipeline executes the `model_inference` rule to generate a deep-learning segmentation of electrode contacts in native CT space.

```text
sub-{subject}/
 └── sub-{subject}_run-01_desc-contacts_nnUNet_dseg.nii.gz
```

### 2. Spatial Registration (`--transform`)

The `--transform` flag triggers the registration of the post-operative CT to the pre-operative T1w MRI. This executes the `get_registration_matrix` and `register_contacts` steps.

```text
sub-{subject}/
├── registration/
│   ├── sub-{subject}_run-01_desc-from_ct-to_T1w_slicer.mat
│   └── sub-{subject}_run-01_desc-from_ct-to_T1w_xfm.txt
├── ses-pre/
│   └── anat/
│       └── sub-{subject}_ses-pre_desc-n4_T1w.nii.gz
├── ses-post/
│   └── anat/
│       └── sub-{subject}_ses-post_run-01_space-T1w_ct.nii.gz
├── slicer_fcsv/
│   ├── sub-{subject}_run-01_contactseg.fcsv
│   └── sub-{subject}_run-01_transformed_contactseg.fcsv
├── sub-{subject}_run-01_desc-contacts_nnUNet_dseg.nii.gz
└── sub-{subject}_run-01_space-T1w_desc-contacts_nnUNet_dseg.nii.gz
```

* **`registration/`**: Contains Slicer-compatible transformation matrices for CT-to-T1w alignment.
* **`ses-pre/anat/`**: Contains the N4 bias-corrected T1w image used as the anatomical reference.
* **`slicer_fcsv/`**: Includes raw and MRI-aligned coordinate files (`.fcsv`) for 3D Slicer visualization.

### 3. Trajectory-Based Numerical Labeling (`--label`)

To produce numerically labeled electrodes in T1w native space, both the `--transform` and `--label` flags are required.

This process utilizes the `sub-{subject}_run-01_planned.fcsv` file as a geometric prior. The pipeline identifies the Entry and Target coordinates for each lead and numbers all segmented contacts found along that trajectory.

```text
sub-{subject}/
├── ses-post/
│   └── ieeg/
│       ├── sub-{subject}_ses-post_run-01_space-T1w_coordsystem.json
│       └── sub-{subject}_ses-post_run-01_space-T1w_electrodes.tsv
├── slicer_fcsv/
│   └── sub-{subject}_run-01_labelled_contactseg.fcsv
└── ... (includes registration outputs)
```

### 4. Anatomical Atlas Labeling (`--atlas_labels`)

The `--atlas_labels` flag is used to automatically assign anatomical regions (e.g., Amygdala, Hippocampus) to the segmented SEEG contacts using standard atlases (like the CerebrA atlas). `contactseg` uses ANTsPy SyN to non-linearly register the native T1w image to the template space.

By default, the pipeline transforms the native contact coordinates into MNI space to extract standard atlas labels. If you prefer to warp the standard atlas back into the patient's native anatomy to extract labels, you can add the `--use_native_space` flag.

To account for minor registration misalignments or contacts located precisely on structural boundaries, the pipeline heavily relies on a `lookup_atlas_labels` function. This performs a fuzzy search to locate the nearest voxel within a specified radius that has a non-zero voxel label. That numerical label is then matched to the atlas TSV file to retrieve the anatomical label name along with whether it is in the left or right hemisphere.

```text
sub-{subject}/
├── atlas/
│   ├── **sub-{subject}_desc-from_T1w-to-MNI_slicer.mat**
│   ├── **sub-{subject}_desc-mni_to_t1_InverseWarp.nii.gz**
│   ├── **sub-{subject}_desc-t1_to_mni_Warp.nii.gz**
│   ├── sub-{subject}_run-01_atlas_in_t1_space.nii.gz
│   ├── sub-{subject}_run-01_atlas_labelled_contactseg.csv
│   ├── sub-{subject}_run-01_atlas_labelled_contactseg.fcsv
│   └── sub-{subject}_run-01_mni_transformed_contactseg.fcsv
└── ... (includes registration and trajectory outputs)
```

**Key Outputs in `atlas/`:**

* **`*_atlas_labelled_contactseg.fcsv`**: Provides the coordinates in native space along with their respective anatomical labels.
* **`*_mni_transformed_contactseg.fcsv`**: Contains the contact coordinates mapped into standard MNI space.
* **`*_atlas_labelled_contactseg.csv`**: A spreadsheet containing the native coordinates, MNI coordinates, anatomical label, numbered label, and tissue type.
* **Transformation Files**: If `contactseg` computes its own registrations, this `atlas/` directory will also contain the non-linear forward (`*Warp.nii.gz`) and inverse warps (`*InverseWarp.nii.gz`) generated by ANTsPy SyN, as well as the Slicer-compatible `.mat` transform (**bolded above**). Let it be noted that if pre-computed sMRIPrep derivatives are provided using the `--SMRIPREP-DIR` flag, `contactseg` will **not** create or use these transformation files; instead, it reuses the standard sMRIPrep `.h5` transforms. Additionally, using sMRIPrep derivatives populates the tissue type column in the spreadsheet, explicitly stating whether the contacts are in white matter, grey matter, or cerebral spinal fluid (CSF).