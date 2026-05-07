rule get_coords:
    input:
        model_seg=rules.model_inference.output.contact_seg,
    output:
        model_coords=bids(
            root=deriv_root,
            datatype="ieeg",
            suffix="coords",
            desc="model",
            session="post",
            space="ct",
            extension=".fcsv",
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
                root=deriv_root,
                datatype="ieeg",
                suffix="coords",
                desc="transformed",
                session="post",
                space="T1w",
                extension=".fcsv",
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
                root=deriv_root,
                datatype="ieeg",
                suffix="coords",
                desc="labeled",
                space="T1w",
                session="post",
                extension=".fcsv",
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
                root=deriv_root,
                datatype="ieeg",
                space="T1w",
                suffix="electrodes",
                session="post",
                extension=".tsv",
                **inputs["post_ct"].wildcards,
            ),
            coordsystem_json=bids(
                root=deriv_root,
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
                root=deriv_root,
                datatype="ieeg",
                suffix="coords",
                extension=".fcsv",
                space="TEMPLATE",
                session="post",
                desc="labeled",
                **inputs["post_ct"].wildcards,
            ),
            native_coords=rules.label_coords.output.labelled_coords,
            atlas_segmentation_in_native=bids(
                root=deriv_root,
                datatype="anat",
                suffix="dseg",
                desc="ATLAS",
                space="T1w",
                extension=".nii.gz",
                **inputs["post_ct"].wildcards,
            ),
            atlas_segmentation_in_mni=str(
                Path(workflow.basedir).parent.parent
                / "resources/atlases/tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii"
            ),
            atlas_labels=str(
                Path(workflow.basedir).parent.parent
                / "resources/atlases/tpl-MNI152NLin2009cSym_atlas-CerebA_dseg.tsv"
            ),
            native_dseg=get_smriprep_dseg,
            native_prob_seg_GM=get_smriprep_probseg("GM"),
            native_prob_seg_WM=get_smriprep_probseg("WM"),
            native_prob_seg_CSF=get_smriprep_probseg("CSF"),
        output:
            atlas_labelled_t1w_contactseg=bids(
                root=deriv_root,
                datatype="ieeg",
                suffix="coords",
                session="post",
                space="T1w",
                desc="anatomical_labels",
                extension=".fcsv",
                **inputs["post_ct"].wildcards,
            ),
            csv_file=bids(
                root=deriv_root,
                datatype="ieeg",
                suffix="coords",
                session="post",
                desc="anatomical_labels",
                extension=".csv",
                **inputs["post_ct"].wildcards,
            ),
        params:
            fuzzy_dist=2,
            GWmatter_labels=config["SMRIPREP_DIR"],
        script:
            "../scripts/lookup_atlas_labels.py"
