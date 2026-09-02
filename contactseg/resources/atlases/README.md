# Atlas resources

Files referenced by the `templates` and `atlases` blocks of
`contactseg/config/snakebids.yml`. Paths in those blocks are resolved relative
to this directory (`atlas_dir`), except for atlases whose `space` is `native`,
which are read from the FreeSurfer subject directory instead.

## Tracked here

| File | Used by |
| --- | --- |
| `tpl-MNI152NLin2009cSym_atlas-CerebrA_dseg.tsv` | `--atlas CerebrA` label names |

## Volumes to add

Nothing in this directory is tracked except the table above and this file.
`.gitignore` excludes the volumes and the annotations, so a patch or a bundle
never carries them and never overwrites what you downloaded. Populating them is
a setup step, done once per checkout.

The two volumes below already exist on the `jthrower/atlas_labels` branch, so
the quickest way to get them is:

```bash
git checkout origin/jthrower/atlas_labels -- resources/atlases/
mv resources/atlases/tpl-MNI152NLin2009cSym_res-1_T1w.nii.gz \
   contactseg/resources/atlases/
gzip -9 resources/atlases/tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii
mv resources/atlases/tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii.gz \
   contactseg/resources/atlases/
rm -r resources/atlases
```

| File | Used by |
| --- | --- |
| `tpl-MNI152NLin2009cSym_res-1_T1w.nii.gz` | fixed image of `reg_t1w_to_template` |
| `tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii.gz` | `--atlas CerebrA` |

Store the CerebrA segmentation gzipped rather than as the raw `.nii` the atlas
ships as — it goes from 17 MB to under 300 kB, and the config names it
`.nii.gz`.

`.nii.gz` is tracked with git-lfs elsewhere in this repository, which is the
reason these two are ignored rather than committed: a bundle or a patch carries
git history but not LFS objects, so a committed volume arrives as a 130-byte
pointer whose content cannot be fetched from anywhere.

## Surface atlases

### Yale Brain Atlas

`--atlas Yale` reads the 696-parcel atlas as the two fsaverage annotations the
project distributes. Put both here:

```bash
cd contactseg/resources/atlases
base=https://raw.githubusercontent.com/YaleBrainAtlas/YaleBrainAtlas/master/data/YBA_696parcels
curl -LO $base/YBA_696_LH_fsaverage.annot
curl -LO $base/YBA_696_RH_fsaverage.annot
```

Roughly 1.3 MB each. The parcel names and colours travel inside the
annotation, so there is no lookup table to fetch: the workflow writes one from
the annotation when it maps the atlas into the volume.

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
