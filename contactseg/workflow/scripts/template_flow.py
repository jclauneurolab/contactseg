from templateflow import api as tflow
import shutil

# Force TemplateFlow to download directly into your output directory
# instead of the hidden ~/.cache folder
template_name = snakemake.params["template"]
atlas_name = snakemake.params.get("atlas")
atlas_desc = snakemake.params.get("atlas_desc")

# 1. Download T1w template
if template_name in ["SUIT"]:
    template_path = tflow.get(str(template_name), desc=None, suffix="T1w", extension="nii.gz")
elif template_name in ["Fischer344"]:
    template_path = tflow.get(str(template_name), desc=None, suffix="T2w", extension="nii.gz")
elif template_name in ["VALiDATe29"]:
    template_path = tflow.get(str(template_name), desc="brain", suffix="T2w", extension="nii.gz")
else:
    template_path = tflow.get(str(template_name), desc=None, resolution=1, suffix="T1w", extension="nii.gz")

# 2. Download Atlas files (required so Snakemake finds them!)
if atlas_name in ["VALiDATe29", "Diedrichsen2009", "v4"]:
    atlas_path = tflow.get(str(template_name), atlas=str(atlas_name), suffix="dseg", extension="nii.gz")
    tsv_path =tflow.get(str(template_name), atlas=str(atlas_name), suffix="dseg",extension="tsv")
elif atlas_name in ["schaefer2018", "HOCPAL", "HOCPA", "HOSPA"]:
    atlas_path = tflow.get(str(template_name), resolution=1, atlas=str(atlas_name), desc=str(atlas_desc), suffix="dseg", extension="nii.gz")
    tsv_path = tflow.get(str(template_name), atlas=str(atlas_name), suffix="dseg",extension="tsv")

else:
    atlas_path = tflow.get(str(template_name), resolution=1, atlas=str(atlas_name), suffix="dseg", extension="nii.gz")
    tsv_path = tflow.get(str(template_name), atlas=str(atlas_name), suffix="dseg", extension="tsv")

# Instead of letting tflow.get decide where to save, 
# you MUST move the result to the Snakemake output paths:
shutil.copy2(str(template_path), snakemake.output.template_t1w)
shutil.copy2(str(atlas_path), snakemake.output.atlas_dseg)
shutil.copy2(str(tsv_path), snakemake.output.atlas_tsv)