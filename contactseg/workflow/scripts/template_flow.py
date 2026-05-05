from templateflow import api as tflow
import sys

template_name = snakemake.params.template
atlas_name = snakemake.params.atlas
atlas_desc = snakemake.params.atlas_desc

default_template_path = snakemake.input.default_template_path
default_atlas_path = snakemake.input.default_atlas_path
default_tsv_path = snakemake.input.default_tsv_path

# If atlas requires a description, fail fast with a helpful message
atlas_requires_desc = {"Schaefer2018", "HOCPAL", "HOSPA", "HOCPA"}
if atlas_name:
    if atlas_name in atlas_requires_desc and not atlas_desc:
        sys.exit(
            f"Error: atlas '{atlas_name}' requires an --atlas-desc (e.g. '400Parcels7Networks' or 'th0'). "
            "Please provide atlas_desc through params or CLI."
        )

# 1. Download T1w template
if template_name:
    template_path = tflow.get(
        str(template_name), desc=None, resolution=1, suffix="T1w", extension="nii.gz"
    )

# 2. Download Atlas files (required so Snakemake finds them!)
if atlas_name in ["HOCPAL", "HOCPA", "HOSPA"]:
    atlas_path = tflow.get(
        str(template_name),
        resolution=1,
        atlas=str(atlas_name),
        desc=str(atlas_desc),
        suffix="dseg",
        extension="nii.gz",
    )
    tsv_path = tflow.get(
        str(template_name),
        atlas=str(atlas_name),
        suffix="dseg",
        extension="tsv",
    )
elif atlas_name == "Schaefer2018":
    atlas_path = tflow.get(
        str(template_name),
        resolution=1,
        atlas=str(atlas_name),
        desc=str(atlas_desc),
        suffix="dseg",
        extension="nii.gz",
    )
    tsv_path = tflow.get(
        str(template_name),
        atlas=str(atlas_name),
        desc=str(atlas_desc),
        suffix="dseg",
        extension="tsv",
    )
elif atlas_name == "MASSP20":
    atlas_path = tflow.get(
        str(template_name),
        resolution=1,
        atlas=str(atlas_name),
        suffix="dseg",
        extension="nii.gz",
    )
    tsv_path = tflow.get(
        str(template_name),
        resolution=1,
        atlas=str(atlas_name),
        suffix="dseg",
        extension="tsv",
    )
else:
    atlas_path = tflow.get(
        str(template_name),
        resolution=1,
        atlas=str(atlas_name),
        suffix="dseg",
        extension="nii.gz",
    )
    tsv_path = tflow.get(
        str(template_name), atlas=str(atlas_name), suffix="dseg", extension="tsv"
    )

print(
    f"Template Path: {str(template_path)}",
    f"Atlas Path: {str(atlas_path)}",
    f"Atlas TSV Path: {str(tsv_path)}",
)
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
