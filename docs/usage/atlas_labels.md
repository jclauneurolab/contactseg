# Atlas labelling

`--atlas_labels` gives every contact an anatomical label. The labelling itself
is always the same operation — look the contact up in a segmentation that
shares its coordinate frame — but the segmentation can be reached three
different ways, chosen with `--atlas_source`.

```bash
contactseg bids_dir out_dir participant --label --atlas_labels
```

`--atlas_labels` implies `--transform`, since the lookup runs on the contact
coordinates in T1w space.

## Sources

### `--atlas_source antspy` (default)

The workflow computes the subject-to-template registration itself, with a
single ANTsPy SyN registration against the packaged template
(`rule reg_t1w_to_template`). Tissue probabilities come from a three-class
Atropos segmentation of the same T1w (`rule tissue_seg`). Nothing outside the
BIDS dataset is needed.

### `--atlas_source smriprep`

```bash
contactseg bids_dir out_dir participant --label --atlas_labels \
    --atlas_source smriprep --smriprep_dir /path/to/derivatives/smriprep
```

No registration is run. The composite transforms sMRIPrep already wrote are
used instead:

```
sub-<label>/anat/sub-<label>_from-T1w_to-<template>_mode-image_xfm.h5
sub-<label>/anat/sub-<label>_from-<template>_to-T1w_mode-image_xfm.h5
sub-<label>/anat/sub-<label>_label-{GM,WM,CSF}_probseg.nii.gz
```

sMRIPrep has to have been run with the same template as `--template`; by
default that is `MNI152NLin2009cSym`, which is the space the CerebrA atlas is
defined in. A `ses-<label>/anat/` layout is searched as well as the
session-less one, with the session taken from `derivatives_session` in the
config.

### `--atlas_source freesurfer`

```bash
contactseg bids_dir out_dir participant --label --atlas_labels \
    --atlas_source freesurfer --freesurfer_dir /path/to/derivatives/freesurfer \
    --atlas aparcaseg DKTatlas
```

Labels come from the subject's own FreeSurfer (or FastSurfer) run, so no
template registration happens at all and no template coordinates are reported.
`aparc+aseg.mgz` is read on its own grid — the conformed volume shares its
scanner RAS frame with the T1w, so the labels are never resampled and never
interpolated. Tissue classes are derived from the label indices themselves
rather than from probability maps.

## Atlases

`--atlas` takes one or more entries from the `atlases` config block. Each
declares where it lives and how it is stored, and the workflow picks the route
that turns it into a segmentation in subject space:

| Atlas | Space | Type | Route into subject space |
| --- | --- | --- | --- |
| `CerebrA` | template | volume | warped with the transforms from `--atlas_source` |
| `Yale` | fsaverage | surface | resampled to the subject sphere, then painted into the volume |
| `DKTatlas` | native | surface | annotation converted to GIFTI, then painted into the volume |
| `aparcaseg` | native | volume | read where it is |

Adding an atlas means adding an entry to that block, not writing a rule.

## Volume and surface representations

Surface atlases are carried into the volume with
`wb_command -label-to-volume-mapping`, ribbon-constrained between the white and
pial surfaces, one hemisphere at a time. `key_offset` shifts each hemisphere's
label values before the two are merged, so parcellations that number their
parcels per hemisphere stay distinguishable — with 1000 and 2000 the result
carries FreeSurfer's own `ctx-lh-*` / `ctx-rh-*` indices.

`--map_atlas_surfaces` runs the mapping the other way as well, so volumetric
atlases get a surface representation via
`wb_command -volume-to-surface-mapping` and `-metric-label-import`. Both
directions need `--freesurfer_dir`, since they are built on the subject's
surfaces.

FreeSurfer surfaces are converted to GIFTI with the c_ras offset of the
conformed volume applied, which puts them in the same scanner RAS frame as the
contacts.

## Uncertainty

The lookup does not stop at the voxel a contact rounds into. The neighbourhood
is sampled with a Gaussian weight and reported as a distribution over
structures, so a contact on a boundary says so:

```
sigma_total = sqrt(sigma_contact**2 + sigma_reg**2)
```

`--sigma_contact` covers the physical extent of the contact plus its
localisation error; `--sigma_reg` covers the residual error of whichever
registration brought the atlas and the contacts together. Raising `--sigma_reg`
for a registration you trust less widens the distribution rather than letting
the pipeline report a crisp wrong answer.

Each lookup prints a concordance line — how many contacts landed on a label,
and how many agree in hemisphere with the L/R prefix of their name. Both
numbers collapsing is the signature of a segmentation that is not in the frame
the workflow thinks it is.

## Outputs

```
sub-<label>/ses-post/ieeg/
  sub-<label>_..._space-T1w_desc-atlas_electrodes.tsv   merged table
  sub-<label>_..._space-T1w_desc-atlas_electrodes.fcsv  labelled contacts for Slicer
  sub-<label>_..._space-T1w_atlas-<atlas>_desc-atlas_labels.tsv
  sub-<label>_..._space-T1w_desc-tissue_labels.tsv
  sub-<label>_..._space-<template>_desc-contactseg_coords.fcsv
sub-<label>/ses-pre/atlasreg/
  sub-<label>_..._from-T1w_to-<template>_mode-image_desc-{affine,warp}_xfm.{mat,nii.gz}
  sub-<label>_..._space-T1w_atlas-<atlas>_dseg.nii.gz
sub-<label>/ses-pre/surf/
  sub-<label>_..._hemi-<L|R>_space-T1w_{white,pial,midthickness}.surf.gii
  sub-<label>_..._hemi-<L|R>_space-T1w_atlas-<atlas>_dseg.label.gii
```

Per-contact columns in the merged table, prefixed by atlas name:
`structure`, `tissue`, `top_structure`, `p_top`, `second_structure`,
`p_second`, `margin`, `entropy`, `confidence`, `p_GM`, `p_WM`, `p_CSF`,
`n_structures`, `dist_to_boundary_mm`, `nearest_structure`,
`nearest_dist_mm`, `hemi_match`.

`--export_nrrd` additionally writes each segmentation as a 3D Slicer
`.seg.nrrd`, with one named, coloured segment per label, so the atlas loads
into the Segmentations module next to the contacts.
