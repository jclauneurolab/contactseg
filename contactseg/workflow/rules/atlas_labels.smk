if needs_template_reg():

    rule warp_contacts_to_template:
        input:
            coords=get_native_coords(),
            transforms=get_transforms("template_to_T1w"),
            xfm_ras=get_bridge_xfm(),
        output:
            warped_coords=get_template_coords(),
        group:
            "subj"
        conda:
            "../envs/image_processing.yaml"
        params:
            invert_flags=get_invert_flags("template_to_T1w"),
        script:
            "../scripts/warp_coords.py"


rule lookup_atlas_labels:
    input:
        coords=get_native_coords(),
        atlas_dseg=get_atlas_dseg_in_native(),
        lut=get_atlas_lut,
    output:
        labels=bids(
            root=config["output_dir"],
            datatype="ieeg",
            session="post",
            space="T1w",
            atlas="{atlas}",
            desc="atlas",
            suffix="labels",
            extension=".tsv",
            **inputs["post_ct"].wildcards,
        ),
    group:
        "subj"
    conda:
        "../envs/analysis.yaml"
    params:
        sigma_contact=config["sigma_contact"],
        sigma_reg=config["sigma_reg"],
        n_sigma=config["n_sigma"],
        max_search=config["max_search_dist"],
    script:
        "../scripts/lookup_atlas_labels.py"


if config["atlas_source"] != "freesurfer":

    rule lookup_tissue_labels:
        input:
            coords=get_native_coords(),
            probseg=get_tissue_probseg(),
            xfm_ras=(
                get_bridge_xfm()
                if config["atlas_source"] == "smriprep"
                else []
            ),
        output:
            tissue=bids(
                root=config["output_dir"],
                datatype="ieeg",
                session="post",
                space="T1w",
                desc="tissue",
                suffix="labels",
                extension=".tsv",
                **inputs["post_ct"].wildcards,
            ),
        group:
            "subj"
        conda:
            "../envs/analysis.yaml"
        params:
            tissue_labels=config["tissue_labels"],
        script:
            "../scripts/lookup_tissue_labels.py"


rule gen_atlas_electrodes:
    input:
        coords=get_native_coords(),
        template_coords=(get_template_coords() if needs_template_reg() else []),
        atlas_labels=expand(
            rules.lookup_atlas_labels.output.labels,
            atlas=config["atlas"],
            allow_missing=True,
        ),
        tissue_labels=(
            rules.lookup_tissue_labels.output.tissue
            if config["atlas_source"] != "freesurfer"
            else []
        ),
    output:
        electrodes_tsv=bids(
            root=config["output_dir"],
            datatype="ieeg",
            session="post",
            space="T1w",
            desc="atlas",
            suffix="electrodes",
            extension=".tsv",
            **inputs["post_ct"].wildcards,
        ),
        electrodes_fcsv=bids(
            root=config["output_dir"],
            datatype="ieeg",
            session="post",
            space="T1w",
            desc="atlas",
            suffix="electrodes",
            extension=".fcsv",
            **inputs["post_ct"].wildcards,
        ),
    group:
        "subj"
    conda:
        "../envs/analysis.yaml"
    params:
        atlases=config["atlas"],
        template=config["template"],
    script:
        "../scripts/gen_atlas_electrodes.py"


if config["export_nrrd"]:

    rule atlas_dseg_to_nrrd:
        input:
            atlas_dseg=get_atlas_dseg_in_native(),
            lut=get_atlas_lut,
        output:
            seg_nrrd=bids(
                root=config["output_dir"],
                datatype="atlasreg",
                session="pre",
                space="T1w",
                atlas="{atlas}",
                suffix="dseg.seg.nrrd",
                **inputs["pre_t1w"].wildcards,
            ),
        group:
            "subj"
        conda:
            "../envs/analysis.yaml"
        params:
            orientation=config["nrrd_orientation"],
        script:
            "../scripts/atlas_dseg_to_nrrd.py"
