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
            electrode_type=str(Path(workflow.basedir).parent / config["electrode_type"]),
        group:
            "subj"
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
        script:
            "../scripts/generate_tsv.py"

if config["atlas_labels"]:

    rule lookup_atlas_labels:
        input:
            mni_coords=bids(  
                root=config["output_dir"],
                suffix="mni_transformed_contactseg.fcsv",
                datatype="atlas",
                **inputs["post_ct"].wildcards,
            ),
            native_coords=bids(
                root=config["output_dir"],
                suffix="labelled_contactseg.fcsv",
                datatype="slicer_fcsv",
                **inputs["post_ct"].wildcards,
            ),
            atlas_segmentation_in_native=bids(
                root=config["output_dir"],
                suffix="atlas_in_t1_space.nii.gz",
                datatype="atlas",
                session="post",
                **inputs["post_ct"].wildcards,
            ),
            atlas_segmentation_in_mni=str(Path(workflow.basedir).parent.parent / "resources/atlases/tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii"),
            atlas_labels=str(Path(workflow.basedir).parent.parent / "resources/atlases/tpl-MNI152NLin2009cSym_atlas-CerebA_dseg.tsv")
        output:
            atlas_labelled_t1w_contactseg=bids(
                root=config["output_dir"],
                suffix="atlas_labelled_contactseg.fcsv",
                datatype="atlas",
                **inputs["post_ct"].wildcards,
            ),
             csv_file=bids(
                root=config["output_dir"],
                suffix="atlas_labelled_contactseg.csv",
                datatype="atlas",
                **inputs["post_ct"].wildcards,
            ),
        params:
            fuzzy_dist=2,
            native_space = config["use_native_space"]
        script:
            "../scripts/lookup_atlas_labels.py"

            