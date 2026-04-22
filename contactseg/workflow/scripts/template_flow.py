from templateflow import api as tflow

# Force TemplateFlow to download directly into your output directory
# instead of the hidden ~/.cache folder
template_name = snakemake.params["template"]
atlas_name = snakemake.params.get("atlas")
atlas_desc = snakemake.params.get("atlas_desc")

# 1. Download T1w template
if template_name in ["SUIT"]:
    tflow.get(str(template_name), desc=None, suffix="T1w", extension="nii.gz")
elif template_name in ["Fischer344"]:
    tflow.get(str(template_name), desc=None, suffix="T2w", extension="nii.gz")
elif template_name in ["VALiDATe29"]:
    tflow.get(str(template_name), desc="brain", suffix="T2w", extension="nii.gz")
else:
    tflow.get(str(template_name), desc=None, resolution=1, suffix="T1w", extension="nii.gz")

# 2. Download Atlas files (required so Snakemake finds them!)
if atlas_name in ["VALiDATe29", "Diedrichsen2009", "v4"]:
    tflow.get(str(template_name), atlas=str(atlas_name), suffix="dseg", extension="nii.gz")
    tflow.get(str(template_name), atlas=str(atlas_name), suffix="dseg",extension="tsv")
elif atlas_name in ["schaefer2018", "HOCPAL", "HOCPA", "HOSPA"]:
    tflow.get(str(template_name), resolution=1, atlas=str(atlas_name), desc=str(atlas_desc), suffix="dseg", extension="nii.gz")
    tflow.get(str(template_name), atlas=str(atlas_name), suffix="dseg",extension="tsv")

else:
    tflow.get(str(template_name), resolution=1, atlas=str(atlas_name), suffix="dseg", extension="nii.gz")
    tflow.get(str(template_name), atlas=str(atlas_name), suffix="dseg", extension="tsv")
    