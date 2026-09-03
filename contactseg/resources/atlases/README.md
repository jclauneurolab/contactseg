# Atlas resources

Files referenced by the `templates` and `atlases` blocks of
`contactseg/config/snakebids.yml`. Paths in those blocks are resolved relative
to this directory (`atlas_dir`), except for atlases whose `space` is `native`,
which are read from the FreeSurfer subject directory instead.

## Tracked here

Everything the packaged atlases need is committed, so a fresh clone can run
`--atlas_labels` without fetching anything by hand:

| File | Used by |
| --- | --- |
| `tpl-MNI152NLin2009cSym_res-1_T1w.nii.gz` | fixed image of `reg_t1w_to_template` |
| `tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii.gz` | `--atlas CerebrA` |
| `tpl-MNI152NLin2009cSym_atlas-CerebrA_dseg.tsv` | `--atlas CerebrA` label names |
| `YBA_696_LH_fsaverage.annot`, `YBA_696_RH_fsaverage.annot` | `--atlas Yale` |

The two `.nii.gz` go through git-lfs, as every `.nii.gz` in this repository
does, so **clone with git-lfs installed** — `git lfs install` once per machine,
then a plain `git clone`. Without it they arrive as 130-byte pointer files and
the workflow fails on the first read. `git lfs pull` fixes an existing clone.

The `.annot` files are ordinary git objects, about 1.3 MB each.

A patch or a bundle carries git history but not lfs objects, so the two volumes
have to come from a real clone or fetch of the remote, not from a `git am` of an
emailed patch.

The CerebrA segmentation is stored gzipped rather than as the raw `.nii` the
atlas ships as — 17 MB down to under 300 kB, and the config names it `.nii.gz`.
No voxel values or labels were altered in any of the files here.

## Attribution

The files in this directory are third-party data redistributed under their own
licences, which are **not** the MIT licence covering contactseg's code. Anyone
using an atlas through this workflow inherits its terms and should cite its
paper, not only this tool.

### CerebrA, and the MNI-ICBM152 2009c template it is defined on

`tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii.gz`,
`tpl-MNI152NLin2009cSym_atlas-CerebrA_dseg.tsv`,
`tpl-MNI152NLin2009cSym_res-1_T1w.nii.gz`

Released under a [Creative Commons Attribution 4.0 International licence
(CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/), which permits
redistribution — including commercially — provided the source is credited, the
licence named, and any changes indicated. Redistributed here unmodified except
for gzip compression.

> Manera, A.L., Dadar, M., Fonov, V. & Collins, D.L. CerebrA, registration and
> manual label correction of Mindboggle-101 atlas for MNI-ICBM152 template.
> *Scientific Data* **7**, 237 (2020).
> https://doi.org/10.1038/s41597-020-0557-9

CerebrA is a manually corrected registration of the Mindboggle-101 labels onto
the MNI-ICBM152 2009c symmetric template. Both underlying works carry their own
citations, and the template ships here as its own file:

> Fonov, V., Evans, A.C., Botteron, K., Almli, C.R., McKinstry, R.C. &
> Collins, D.L. Unbiased average age-appropriate atlases for pediatric studies.
> *NeuroImage* **54**(1), 313–327 (2011).
> https://doi.org/10.1016/j.neuroimage.2010.07.033

> Klein, A. & Tourville, J. 101 labeled brain images and a consistent human
> cortical labeling protocol. *Frontiers in Neuroscience* **6**, 171 (2012).
> https://doi.org/10.3389/fnins.2012.00171

### Yale Brain Atlas

`YBA_696_LH_fsaverage.annot`, `YBA_696_RH_fsaverage.annot`

> McGrath, H., Zaveri, H.P., Collins, E., Jafar, T., Chishti, O., Obaid, S.,
> Ksendzovsky, A., Wu, K., Papademetris, X. & Spencer, D.D. High-resolution
> cortical parcellation based on conserved brain landmarks for localization of
> multimodal data to the nearest centimeter. *Scientific Reports* **12**, 18778
> (2022). https://doi.org/10.1038/s41598-022-21543-3

Project site: https://yalebrainatlas.github.io/YaleBrainAtlas/ — the authors ask
that it be linked wherever the atlas is used.

Note the licence is stated two different ways by the two places the atlas is
distributed from: the [NITRC entry](https://www.nitrc.org/projects/yale_atlas_2021/)
records **Attribution Non-Commercial**, while the
[GitHub repository](https://github.com/YaleBrainAtlas/YaleBrainAtlas) declares
**MIT** for the repository as a whole. This directory treats the stricter of the
two as binding: the annotations are redistributed for academic and research use
with attribution, and anyone intending commercial use should settle the question
with the atlas authors first rather than reading it off either page.

## Surface atlases

### Yale Brain Atlas

`--atlas Yale` reads the 696-parcel atlas as the two fsaverage annotations the
project distributes, both committed here. They were fetched with:

```bash
cd contactseg/resources/atlases
base=https://raw.githubusercontent.com/YaleBrainAtlas/YaleBrainAtlas/master/data/YBA_696parcels
curl -LO $base/YBA_696_LH_fsaverage.annot
curl -LO $base/YBA_696_RH_fsaverage.annot
```

The parcel names and colours travel inside the annotation, so there is no
lookup table to fetch: the workflow writes one from the annotation when it maps
the atlas into the volume.

If the right-hemisphere file is named differently in the repository, correct
the pattern in the `atlases` block rather than renaming the download:

```yaml
  Yale:
    label: "YBA_696_{hemi_up}_fsaverage.annot"
```

`{hemi_up}` expands to LH/RH, `{hemi_fs}` to lh/rh and `{hemi}` to L/R.

The two hemispheres number their parcels from the same range, so the config
shifts the right by `key_offset: {L: 0, R: 1000}` before the hemispheres are
merged into one volume. If the atlas is ever redistributed with hemisphere-
unique numbering, set both offsets to zero.

### fsaverage

Resampling an fsaverage atlas onto a subject needs fsaverage's own
registration sphere. It is not redistributed here either — it ships with
freesurfer, and is looked for in `<freesurfer_dir>/fsaverage/` first and
`$FREESURFER_HOME/subjects/fsaverage/` second. If neither has it, copy the
fsaverage subject into your freesurfer derivatives directory.

### Adding another surface atlas

Give it an entry under `atlases` with `space: fsaverage`, `type: surface` and
its own `label` file. Both `.annot` and `.label.gii` are accepted; an
annotation is converted once and shared across subjects.
