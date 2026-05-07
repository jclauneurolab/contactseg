rule register_contacts:
    input:
        in_im=bids(
            root=deriv_root,
            suffix="dseg.nii.gz",
            desc="contacts_nnUNet",
            session="post",
            **inputs["post_ct"].wildcards,
        ),
        ref_im=rules.n4biascorr.output.corrected_t1w,
        transform_matrix=get_reg_matrix(),
    output:
        out_im=bids(
            root=deriv_root,
            suffix="dseg.nii.gz",
            space="T1w",
            desc="contacts_nnUNet",
            session="post",
            **inputs["post_ct"].wildcards,
        ),
    script:
        "../scripts/apply_registration.py"


rule contacts_qc:
    input:
        ct_img=get_registered_ct_image(),
        t1w_img=rules.n4biascorr.output.corrected_t1w,
        contact_fcsv_labelled=bids(
            root=deriv_root,
            datatype="ieeg",           
            suffix="coords",
            session="post",
            space="T1w",
            extension=".fcsv",
            desc="labeled",
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
            root=deriv_root,
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
                root=deriv_root,
                datatype="anat",
                space="TEMPLATE",
                suffix="T1w",
                extension=".nii.gz",
                **inputs["pre_t1w"].wildcards,
            ),
            mni_template=str(
                Path(workflow.basedir).parent.parent
                / "resources/atlases/tpl-MNI152NLin2009cSym_res-1_T1w.nii.gz"
            ),
        output:
            forward_warp=bids(
                root=deriv_root,
                datatype="anat",
                desc="T1w2MNI",
                suffix="xfm",
                extension=".nii.gz",
                **inputs["pre_t1w"].wildcards,
            ),
            inverse_warp=bids(
                root=deriv_root,
                datatype="anat",
                desc="MNI2T1w",
                suffix="xfm",
                extension=".nii.gz",
                **inputs["pre_t1w"].wildcards,
            ),
            affine_syn=bids(
                root=deriv_root,
                datatype="anat",
                desc="nonlinear_T1w2MNI",
                suffix="xfm",
                extension=".mat",
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
                root=deriv_root,
                datatype="anat",
                desc="T1w2MNI",
                suffix="xfm",
                extension=".nii.gz",
                **inputs["pre_t1w"].wildcards,
            ),
            bids(
                root=deriv_root,
                datatype="anat",
                desc="nonlinear_T1w2MNI",
                suffix="xfm",
                extension=".mat",
                **inputs["pre_t1w"].wildcards,
            ),
            bids(
                root=deriv_root,
                datatype="anat",
                desc="from_T1w-to-TEMPLATE",
                suffix="xfm",
                extension=".mat",
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
                root=deriv_root,
                datatype="anat",
                desc="from_T1w-to-TEMPLATE",
                suffix="xfm",
                extension=".mat",
                **inputs["pre_t1w"].wildcards,
            ),
            bids(
                root=deriv_root,
                datatype="anat",
                desc="nonlinear_T1w2MNI",
                suffix="xfm",
                extension=".mat",
                **inputs["pre_t1w"].wildcards,
            ),
            bids(
                root=deriv_root,
                datatype="anat",
                desc="MNI2T1w",
                suffix="xfm",
                extension=".nii.gz",
                **inputs["pre_t1w"].wildcards,
            ),
        ]


rule apply_full_transformation:
    input:
        coords=bids(
            root=deriv_root,
            datatype="ieeg",
            suffix="coords",
            desc="labeled",
            session="post",
            space="T1w",
            extension=".fcsv",
            **inputs["post_ct"].wildcards,
        ),
        transforms=get_forward_transforms,
    output:
        output_coords=bids(
            root=deriv_root,
            datatype="ieeg",
            suffix="coords",
            extension=".fcsv",
            session="post",
            space="TEMPLATE",
            desc="labeled",
            **inputs["post_ct"].wildcards,
        ),
    script:
        "../scripts/apply_full_transformation.py"


rule transform_atlas_to_t1:
    input:
        atlas_image=str(
            Path(workflow.basedir).parent.parent
            / "resources/atlases/tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii"
        ),
        t1_image=rules.n4biascorr.output.corrected_t1w,
        transforms=get_inverse_transforms,
    output:
        atlas_in_t1=bids(
            root=deriv_root,
            datatype="anat",
            suffix="dseg",
            extension=".nii.gz",
            desc="ATLAS",
            space="T1w",
            **inputs["post_ct"].wildcards,
        ),
    script:
        "../scripts/transform_atlas_to_t1.py"
