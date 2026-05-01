from templateflow import api as tflow

template_name=snakemake.params.template
atlas_name=snakemake.params.atlas
atlas_desc=snakemake.params.atlas_desc

default_template_path= snakemake.input.default_template_path
default_atlas_path= snakemake.input.default_atlas_path
default_tsv_path= snakemake.input.default_tsv_path

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

with open(snakemake.output.template_txt, "w") as f:
            if not template_path:
                template_path = default_template_path
            if not atlas_path:
                atlas_path = default_atlas_path
            if not tsv_path:
                tsv_path = default_tsv_path
            f.write(f"Template Path: {str(template_path)}\n")
            f.write(f"Atlas Path: {str(atlas_path)}\n")
            f.write(f"Atlas TSV Path: {str(tsv_path)}\n")