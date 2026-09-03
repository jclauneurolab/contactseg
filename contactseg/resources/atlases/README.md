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

### Provenance

| File | Source |
| --- | --- |
| MNI152NLin2009cSym T1w, CerebrA | the CerebrA atlas release, `nist.mni.mcgill.ca/cerebra` |
| YBA 696 parcels | `github.com/YaleBrainAtlas/YaleBrainAtlas`, `data/YBA_696parcels` |

The CerebrA segmentation is stored gzipped rather than as the raw `.nii` the
atlas ships as — 17 MB down to under 300 kB, and the config names it `.nii.gz`.

Check each atlas's own licence before redistributing this directory outside the
lab; they are not all under the same terms as the code.

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
