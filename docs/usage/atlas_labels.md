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
scanner RAS frame with the image FreeSurfer was run on, so the labels are never
resampled and never interpolated (see the next section when that image is not
the one contactseg localises against). Tissue classes are derived from the label indices themselves
rather than from probability maps.

## When the derivatives were run on a different image

sMRIPrep and FreeSurfer are often run on a different acquisition than the one
contactseg localises contacts against — a non-contrast T1w for the same
patient, while `pre_t1w` points at the contrast-enhanced run. Everything those
pipelines produce (transforms, segmentations, surfaces) then lives in *that*
image's scanner RAS frame, offset from `pre_t1w` by however much the patient
moved between the two scans. Reading them as if the frames agreed puts every
contact off by that offset, silently.

The image each dataset was computed on is named in the config. A bare string is
resolved inside the derivatives dataset:

```yaml
derivatives_anat:
  smriprep: "anat/sub-{subject}_desc-preproc_T1w.nii.gz"
  freesurfer: "mri/orig.mgz"
```

If the image the pipeline was run on is the raw acquisition in your BIDS
dataset — the non-contrast T1w sitting next to the contrast-enhanced run that
`pre_t1w` points at — name it with a `root` instead:

```yaml
derivatives_anat:
  freesurfer:
    root: bids
    path: "anat/sub-{subject}_ses-{session}_run-01_T1w.nii.gz"
```

Paths are relative to the subject (or subject/session) directory of whichever
root, so both forms read the same way, and `{subject}` and `{session}` are
filled in per subject (`{session}` comes from `derivatives_session`).

Which to prefer: the derivatives dataset's own reference image, when it has
one. sMRIPrep's `desc-preproc_T1w` and FreeSurfer's `orig.mgz` *are*, by
definition, the frames their transforms and segmentations are expressed in. A
raw acquisition matches that frame only when the pipeline consumed exactly that
one image — feed sMRIPrep two T1w runs and it averages them, and the averaged
reference is not identical to either input. Point at the BIDS image when the
derivatives directory does not ship a reference image, or when you want the
bridge pinned to a specific acquisition and know only that one was used.

`--derivatives_reg rigid` (the default) registers that image to the contactseg
T1w once per dataset and writes a 4x4 RAS matrix, in exactly the convention the
workflow already uses for the CT-to-T1w matrix:

```
sub-<label>/ses-pre/atlasreg/
  sub-<label>_..._from-<smriprep|freesurfer>_to-T1w_mode-image_xfm.txt
  sub-<label>_..._space-T1w_desc-<smriprep|freesurfer>_T1w.nii.gz   (QC)
```

Every consumer then applies it the cheapest exact way there is:

| Consumer | How the bridge is applied |
| --- | --- |
| FreeSurfer segmentations | composed into the image affine — no voxel is resampled |
| FreeSurfer surfaces | applied to the vertices, after the c_ras shift |
| Contacts going to template | applied in RAS before the ANTs chain, as `transform_coords` does |
| Template atlas coming back | ANTs resamples onto the derivatives grid, then the affine is composed |
| sMRIPrep tissue maps | its inverse is applied to the contacts before sampling |

Nothing is interpolated twice, and the nonlinear ANTs chain stays the
two-transform form that the tests cover.

A run that takes transforms from sMRIPrep and surfaces from FreeSurfer is
reading two frames, so it gets two bridges — one per dataset.

Set `--derivatives_reg identity` to skip it, which is correct only when the
derivatives were computed on the very same image as `pre_t1w`. If you are not
sure, leave it on: registering an image to itself costs about a minute and
returns an identity, whereas guessing wrong is the failure the concordance
check exists to catch.

Point `pre_t1w` at whichever acquisition you want the contacts localised in by
editing `pybids_inputs.pre_t1w` in the config — the current filter pins
`run: "02"`.

## Atlases

`--atlas` takes one or more entries from the `atlases` config block. Each
declares where it lives and how it is stored, and the workflow picks the route
that turns it into a segmentation in subject space:

| Atlas | Space | Type | Route into subject space |
| --- | --- | --- | --- |
| `CerebrA` | template | volume | warped with the transforms from `--atlas_source` |
| `Yale` | fsaverage | surface | annotation converted, resampled to the subject sphere, then painted into the volume |
| `DKTatlas` | native | surface | annotation converted to GIFTI, then painted into the volume |
| `aparcaseg` | native | volume | read where it is |

Adding an atlas means adding an entry to that block, not writing a rule. A
surface atlas may ship as a freesurfer `.annot` or as a `.label.gii`; an
annotation is converted once and shared across subjects. The Yale Brain Atlas
is the worked example — two fsaverage annotations dropped into
`contactseg/resources/atlases/`, with no lookup table to fetch, since the
parcel names travel inside the annotation. See the README there.

Resampling an fsaverage atlas needs fsaverage's registration sphere, taken
from `<freesurfer_dir>/fsaverage/` or `$FREESURFER_HOME/subjects/fsaverage/`
rather than being redistributed.

## Volume and surface representations

An fsaverage atlas is resampled onto the subject with
`wb_command -label-resample ... ADAP_BARY_AREA -area-surfs`, which weights each
parcel by the cortical area it covers. Plain `BARYCENTRIC` ignores that, which
matters for a fine parcellation: small parcels can otherwise be lost where
fsaverage and the subject differ in surface area.

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
`structure`, `tissue`, `top_structure`, `probability`, `entropy`,
`confidence`, `second_structure`, `p_second`, `p_GM`, `p_WM`, `p_CSF`,
`n_structures`, `dist_to_boundary_mm`, `nearest_structure`,
`nearest_dist_mm`.

Hemisphere agreement — whether the L/R prefix of a contact name matches the
hemisphere of the structure it landed in — is reported in the concordance line
rather than carried as a column. It is a check on the run, not a property of a
contact.

Every atlas that carries a lookup table also gets viewer colour tables, so a
plain labelmap opens as anatomy rather than as numbers:

```
sub-<label>_..._space-T1w_atlas-<atlas>_desc-colors_dseg.ctbl   3D Slicer
sub-<label>_..._space-T1w_atlas-<atlas>_desc-colors_dseg.txt    ITK-SNAP
```

In Slicer, load the `.ctbl` in the Colors module, then set it as the colour
node of the loaded `dseg.nii.gz`. Colours come from the atlas itself where it
has them — a GIFTI label table, or r/g/b columns in the lookup table — and are
otherwise generated from the label index, so the same atlas always looks the
same. `aparcaseg` is the exception: it uses freesurfer's indices, and both
Slicer and freesurfer already ship `FreeSurferColorLUT.txt` for it.

`--export_nrrd` additionally writes each segmentation as a 3D Slicer
`.seg.nrrd`, with one named, coloured segment per label, so the atlas loads
into the Segmentations module next to the contacts.
