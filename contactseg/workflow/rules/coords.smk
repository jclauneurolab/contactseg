rule get_coords:
    input:
        model_seg=rules.model_inference.output.contact_seg,
    output:
        model_coords=bids(
            root=config["output_dir"],
            suffix="contactseg.fcsv",
            datatype="slicer_fcsv",
            **inputs["post_ct"].wildcards,
        ),
    group:
        "subj"
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/nnUNet_coords.py"


if config["transform"]:

    rule transform_coords:
        input:
            coords=rules.get_coords.output.model_coords,
            transformation_matrix=get_reg_matrix(),
        output:
            transformed_coords=bids(
                root=config["output_dir"],
                suffix="transformed_contactseg.fcsv",
                datatype="slicer_fcsv",
                **inputs["post_ct"].wildcards,
            ),
        group:
            "subj"
        conda:
            "../envs/analysis.yaml"
        script:
            "../scripts/transform_coords.py"


if config["label"]:

    rule label_coords:
        input:
            coords=rules.transform_coords.output.transformed_coords,
            planned_fcsv=bids(
                root=config["bids_dir"],
                suffix="planned",
                extension=".fcsv",
                **inputs["post_ct"].wildcards,
            ),
        output:
            labelled_coords=bids(
                root=config["output_dir"],
                suffix="labelled_contactseg.fcsv",
                datatype="slicer_fcsv",
                **inputs["post_ct"].wildcards,
            ),
        params:
            electrode_type=str(
                Path(workflow.basedir).parent / config["electrode_type"]
            ),
        group:
            "subj"
        conda:
            "../envs/analysis.yaml"
        script:
            "../scripts/label.py"

    rule gen_labelled_ieeg_electrodes:
        input:
            fcsv=rules.label_coords.output.labelled_coords,
            ref_ct=get_registered_ct_image(),
        output:
            electrodes_tsv=bids(
                root=config["output_dir"],
                datatype="ieeg",
                space="T1w",
                suffix="electrodes",
                session="post",
                extension=".tsv",
                **inputs["post_ct"].wildcards,
            ),
            coordsystem_json=bids(
                root=config["output_dir"],
                datatype="ieeg",
                space="T1w",
                suffix="coordsystem",
                session="post",
                extension=".json",
                **inputs["post_ct"].wildcards,
            ),
        conda:
            "../envs/analysis.yaml"
        script:
            "../scripts/generate_tsv.py"