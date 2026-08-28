if needs_template_reg() and config["atlas_source"] != "smriprep":

    rule reg_t1w_to_template:
        input:
            t1w=rules.n4biascorr.output.corrected_t1w,
            template=get_template_t1w(),
        output:
            xfm_affine=get_affine_xfm(),
            warp=get_warp_xfm(),
            invwarp=get_warp_xfm(inverse=True),
            xfm_ras=bids(
                root=config["output_dir"],
                datatype="atlasreg",
                session="pre",
                desc="affine",
                suffix="xfm",
                extension=".txt",
                **{"from": "T1w", "to": config["template"]},
                **inputs["pre_t1w"].wildcards,
            ),
            warped_t1w=bids(
                root=config["output_dir"],
                datatype="anat",
                session="pre",
                space=config["template"],
                desc="n4",
                suffix="T1w",
                extension=".nii.gz",
                **inputs["pre_t1w"].wildcards,
            ),
        group:
            "subj"
        conda:
            "../envs/image_processing.yaml"
        threads: 8
        resources:
            mem_mb=16000,
            time=60,
        script:
            "../scripts/template_registration.py"


if get_atlases_by_kind("template", "volume"):

    rule warp_atlas_to_native:
        input:
            atlas_dseg=get_template_atlas_dseg,
            t1w=rules.n4biascorr.output.corrected_t1w,
            transforms=get_transforms("template_to_T1w"),
        output:
            atlas_dseg=get_atlas_dseg_in_native(),
        wildcard_constraints:
            atlas="|".join(get_atlases_by_kind("template", "volume")),
        group:
            "subj"
        conda:
            "../envs/image_processing.yaml"
        params:
            invert_flags=get_invert_flags("template_to_T1w"),
        script:
            "../scripts/warp_atlas_to_native.py"


if get_atlases_by_kind("native", "volume"):

    rule import_fs_atlas:
        input:
            atlas_dseg=get_fs_atlas_dseg,
        output:
            atlas_dseg=get_atlas_dseg_in_native(),
        wildcard_constraints:
            atlas="|".join(get_atlases_by_kind("native", "volume")),
        group:
            "subj"
        conda:
            "../envs/analysis.yaml"
        script:
            "../scripts/import_fs_atlas.py"


if config["atlas_source"] == "antspy":

    rule tissue_seg:
        input:
            t1w=rules.n4biascorr.output.corrected_t1w,
        output:
            probseg=expand(
                bids(
                    root=config["output_dir"],
                    datatype="atlasreg",
                    session="pre",
                    label="{label}",
                    desc="atropos3seg",
                    suffix="probseg",
                    extension=".nii.gz",
                    **inputs["pre_t1w"].wildcards,
                ),
                label=config["tissue_labels"],
                allow_missing=True,
            ),
        group:
            "subj"
        conda:
            "../envs/image_processing.yaml"
        threads: 4
        resources:
            mem_mb=16000,
        params:
            tissue_labels=config["tissue_labels"],
        script:
            "../scripts/tissue_seg.py"
