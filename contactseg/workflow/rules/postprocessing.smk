rule register_contacts:
    input:
        in_im=bids(
            root=config["output_dir"],
            suffix="dseg.nii.gz",
            desc="contacts_nnUNet",
            **inputs["post_ct"].wildcards,
        ),
        ref_im=bids(
            root=config["output_dir"],
            suffix="T1w",
            desc="n4",
            datatype="anat",
            session="pre",
            extension=".nii.gz",
            **inputs["pre_t1w"].wildcards,
        ),
        transform_matrix=get_reg_matrix(),
    output:
        out_im=bids(
            root=config["output_dir"],
            suffix="dseg.nii.gz",
            space="T1w",
            desc="contacts_nnUNet",
            **inputs["post_ct"].wildcards,
        ),
    script:
        "../scripts/apply_registration.py"


rule contacts_qc:
    input:
        ct_img=get_registered_ct_image(),
        t1w_img=bids(
            root=config["output_dir"],
            suffix="T1w",
            desc="n4",
            datatype="anat",
            session="pre",
            extension=".nii.gz",
            **inputs["pre_t1w"].wildcards,
        ),
        contact_fcsv_labelled=bids(
            root=config["output_dir"],
            datatype="slicer_fcsv",
            suffix="labelled_contactseg.fcsv",
            **inputs["post_ct"].wildcards,
        ),
        planned_fcsv=bids(
            root=config["bids_dir"],
            suffix="planned",
            extension=".fcsv",
            **inputs["post_ct"].wildcards,
        ),
    output:
        html=bids(
            root=config["output_dir"],
            datatype="qc",
            desc="contactseg",
            suffix="qc.html",
            **inputs["post_ct"].wildcards,
        ),
    script:
        "../scripts/contacts_qc.py"


if not config.get("SMRIPREP_DIR"):

    rule nonlinear_t1_to_mni:
        input:
            t1_in_mni_img=bids(
                root=config["output_dir"],
                datatype="atlas",
                space="MNI",
                suffix="T1w.nii.gz",
                **inputs["pre_t1w"].wildcards,
            ),
            templateflow_paths="resources/templateflow_template.txt",
        output:
            forward_warp=bids(
                root=config["output_dir"],
                desc="t1_to_mni",
                suffix="Warp.nii.gz",
                datatype="atlas",
                **inputs["pre_t1w"].wildcards,
            ),
            inverse_warp=bids(
                root=config["output_dir"],
                desc="mni_to_t1",
                suffix="InverseWarp.nii.gz",
                datatype="atlas",
                **inputs["pre_t1w"].wildcards,
            ),
            affine_syn=bids(
                root=config["output_dir"],
                desc="t1_to_mni",
                suffix="Affine.mat",
                datatype="atlas",
                **inputs["pre_t1w"].wildcards,
            ),
        script:
            "../scripts/nonlinear_registration.py"


def get_forward_transforms(wildcards):
    smriprep_dir = config.get("SMRIPREP_DIR") or config.get("SMRIPREP-DIR")
    session = getattr(wildcards, "session", "pre")
    if smriprep_dir:
        return [
            f"{smriprep_dir}/sub-{wildcards.subject}/ses-{session}/anat/sub-{wildcards.subject}_ses-{session}_from-MNI152NLin2009cSym_to-T1w_mode-image_xfm.h5"
        ]
    else:
        return [
            bids(
                root=config["output_dir"],
                desc="t1_to_mni",
                suffix="Warp.nii.gz",
                datatype="atlas",
                **inputs["pre_t1w"].wildcards,
            ),
            bids(
                root=config["output_dir"],
                desc="t1_to_mni",
                suffix="Affine.mat",
                datatype="atlas",
                **inputs["pre_t1w"].wildcards,
            ),
            bids(
                root=config["output_dir"],
                datatype="atlas",
                desc="from_T1w-to-MNI",
                suffix="slicer.mat",
                **inputs["pre_t1w"].wildcards,
            ),
        ]


def get_inverse_transforms(wildcards):
    smriprep_dir = config.get("SMRIPREP_DIR") or config.get("SMRIPREP-DIR")
    session = getattr(wildcards, "session", "pre")
    if smriprep_dir:
        return [
            f"{smriprep_dir}/sub-{wildcards.subject}/ses-{session}/anat/sub-{wildcards.subject}_ses-{session}_from-MNI152NLin2009cSym_to-T1w_mode-image_xfm.h5"
        ]
    else:
        return [
            bids(
                root=config["output_dir"],
                datatype="atlas",
                desc="from_T1w-to-MNI",
                suffix="slicer.mat",
                **inputs["pre_t1w"].wildcards,
            ),
            bids(
                root=config["output_dir"],
                desc="t1_to_mni",
                suffix="Affine.mat",
                datatype="atlas",
                **inputs["pre_t1w"].wildcards,
            ),
            bids(
                root=config["output_dir"],
                desc="mni_to_t1",
                suffix="InverseWarp.nii.gz",
                datatype="atlas",
                **inputs["pre_t1w"].wildcards,
            ),
        ]


rule apply_full_transformation:
    input:
        coords=bids(
            root=config["output_dir"],
            suffix="labelled_contactseg.fcsv",
            datatype="slicer_fcsv",
            **inputs["post_ct"].wildcards,
        ),
        transforms=get_forward_transforms,
    output:
        output_coords=bids(
            root=config["output_dir"],
            suffix="mni_transformed_contactseg.fcsv",
            datatype="atlas",
            **inputs["post_ct"].wildcards,
        ),
    script:
        "../scripts/apply_full_transformation.py"


rule transform_atlas_to_t1:
    input:
        templateflow_paths="resources/templateflow_template.txt",
        t1_image=bids(
            root=config["output_dir"],
            suffix="T1w",
            desc="n4",
            datatype="anat",
            session="pre",
            extension=".nii.gz",
            **inputs["pre_t1w"].wildcards,
        ),
        transforms=get_inverse_transforms,
    output:
        atlas_in_t1=bids(
            root=config["output_dir"],
            suffix="atlas_in_t1_space.nii.gz",
            datatype="atlas",
            **inputs["post_ct"].wildcards,
        ),
    script:
        "../scripts/transform_atlas_to_t1.py"
