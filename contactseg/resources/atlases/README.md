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

The two volumes below are large and are not carried in this patch. They already
exist on the `jthrower/atlas_labels` branch, so the quickest way to populate
them is:

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

`.nii.gz` is tracked with git-lfs in this repository, and the CerebrA
segmentation compresses from 17 MB to under 300 kB, so store it gzipped rather
than as the raw `.nii` the atlas ships as.

## Surface atlases

`--atlas Yale` expects the Yale Brain Atlas as per-hemisphere GIFTI label files
on the fsaverage surface, together with the fsaverage sphere they are defined
on:

```
tpl-fsaverage_hemi-L_atlas-Yale_dseg.label.gii
tpl-fsaverage_hemi-R_atlas-Yale_dseg.label.gii
tpl-fsaverage_hemi-L_sphere.surf.gii
tpl-fsaverage_hemi-R_sphere.surf.gii
tpl-fsaverage_atlas-Yale_dseg.tsv
```

These are not redistributed here. If the two hemispheres number their parcels
from the same range, set `key_offset` for the atlas in the config so that the
merged volume keeps them apart — the offsets are added to the label values
before the hemispheres are combined, and the lookup table has to agree with the
result.

Any other surface atlas can be added the same way: give it an entry under
`atlases` with `space: fsaverage`, `type: surface`, and its own `label`,
`sphere` and `lut` files.
