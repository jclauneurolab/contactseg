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
    template_name = config.get("template_flow") or config.get("template-flow", False)
    atlas_name = config.get("template_atlas") or config.get("template-atlas", False)

    if template_name and atlas_name:
        rule download_templateflow:
            output:
                atlas_dseg=Path(config["output_dir"]) / f"templateflow/tpl-{template_name}_atlas-{atlas_name}_dseg.nii.gz",
                atlas_tsv=Path(config["output_dir"]) / f"templateflow/tpl-{template_name}_atlas-{atlas_name}_dseg.tsv",
                template_t1w=Path(config["output_dir"]) / f"templateflow/tpl-{template_name}_T1w.nii.gz"
            params:
                template=template_name,
                atlas=atlas_name,
                out_dir=str(Path(config["output_dir"]) / "templateflow")
            conda:
                "../envs/templateflow.yaml"
            script:
                "../scripts/template_flow.py"


    def get_smriprep_dseg(wildcards):
        smriprep_dir = config.get("SMRIPREP_DIR") or config.get("SMRIPREP-DIR")
        if smriprep_dir:
            session = getattr(wildcards, "session", "pre") 
            return f"{smriprep_dir}/sub-{wildcards.subject}/ses-{session}/anat/sub-{wildcards.subject}_ses-{session}_dseg.nii.gz"
        return [] 
    def get_smriprep_probseg(label):
        def get_probseg(wildcards):
            smriprep_dir = config.get("SMRIPREP_DIR") or config.get("SMRIPREP-DIR")
            if smriprep_dir:
                session = getattr(wildcards, "session", "pre") 
                return f"{smriprep_dir}/sub-{wildcards.subject}/ses-{session}/anat/sub-{wildcards.subject}_ses-{session}_label-{label}_probseg.nii.gz"
            return []
        return get_probseg


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
                **inputs["post_ct"].wildcards,
            ),
            atlas_segmentation_in_mni=get_template_atlas(),
            atlas_labels=get_template_atlas_tsv(),
            native_dseg=get_smriprep_dseg,
            native_prob_seg_GM =get_smriprep_probseg("GM"),
            native_prob_seg_WM =get_smriprep_probseg("WM"),
            native_prob_seg_CSF =get_smriprep_probseg("CSF"),

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
            native_space = config["use_native_space"],
            GWmatter_labels = config["SMRIPREP_DIR"]
        script:
            "../scripts/lookup_atlas_labels.py"

            