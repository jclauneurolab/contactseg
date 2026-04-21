from templateflow import api as tflow

template_name = snakemake.params["template"]
atlas_name = snakemake.params.get("atlas")

# 1. Download T1w template
if template_name in ["MNI305", "MNIColin27"]:
    tflow.get(str(template_name), desc=None, suffix="T1w", extension="nii.gz")
else:
    tflow.get(str(template_name), desc=None, resolution=1, suffix="T1w", extension="nii.gz")

# 2. Download Atlas files (required so Snakemake finds them!)
if atlas_name:
    tflow.get(str(template_name), atlas=str(atlas_name), suffix="dseg", extension="nii.gz")
    tflow.get(str(template_name), atlas=str(atlas_name), extension="tsv")