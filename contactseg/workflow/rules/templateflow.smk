# --- 1. Define the logic as a single dictionary or variables ---
# We use .get() with defaults to avoid the 'None' naming issue during initial parsing

tf_temp = config.get("template-flow", "MNI152NLin2009cAsym")
tf_atl  = config.get("template-atlas", "HOCPA")
tf_d    = config.get("atlas-desc", "th0")
out     = config.get("output_dir", "results")

# Pre-calculate the concrete strings
res_tag = "" if tf_atl in ["VALiDATe29", "Diedrichsen2009", "v4"] else "_res-01"
d_tag   = f"_desc-{tf_d}" if tf_d else ""

target_t1w   = os.path.join(out, "templateflow", f"tpl-{tf_temp}_res-01_T1w.nii.gz")
target_atlas_nii = os.path.join(out, "templateflow", f"tpl-{tf_temp}{res_tag}_atlas-{tf_atl}{d_tag}_dseg.nii.gz")
target_atlas_tsv   = os.path.join(out, "templateflow", f"tpl-{tf_temp}_atlas-{tf_atl}_dseg.tsv")

# --- 2. Use these CONCRETE STRINGS in the output ---
rule download_templateflow:
    params:
        template = tf_temp,
        atlas = tf_atl,
        atlas_desc = tf_d
    output:
        template_t1w = target_t1w,
        atlas_dseg   = target_atlas_nii,
        atlas_tsv    = target_atlas_tsv
    conda:
        "../envs/templateflow.yaml"
    script:
        "../scripts/template_flow.py"