rule n4biascorr:
    input:
        t1w=bids(
            root=config["bids_dir"],
            suffix="T1w",
            session="pre",
            run="02",
            datatype="anat",
            extension=".nii.gz",
            **inputs["pre_t1w"].wildcards,
        ),
    output:
        corrected_t1w=bids(
            root=config["output_dir"],
            suffix="T1w",
            desc="n4",
            datatype="anat",
            session="pre",
            extension=".nii.gz",
            **inputs["pre_t1w"].wildcards,
        ),
    script:
        "../scripts/n4_bias_corr.py"


rule get_registration_matrix:
    input:
        post_ct=bids(
            root=config["bids_dir"],
            suffix="ct",
            datatype="ct",
            session="post",
            acq="Electrode",
            extension=".nii.gz",
            **inputs["post_ct"].wildcards,
        ),
        fixed_t1w=bids(
            root=config["output_dir"],
            suffix="T1w",
            desc="n4",
            datatype="anat",
            session="pre",
            extension=".nii.gz",
            **inputs["pre_t1w"].wildcards,
        ),
    output:
        xfm_slicer=bids(
            root=config["output_dir"],
            datatype="registration",
            desc="from_ct-to_T1w",
            suffix="slicer.mat",
            **inputs["post_ct"].wildcards,
        ),
        xfm_ras=bids(
            root=config["output_dir"],
            datatype="registration",
            desc="from_ct-to_T1w",
            suffix="xfm.txt",
            **inputs["post_ct"].wildcards,
        ),
        out_im=bids(
            root=config["output_dir"],
            datatype="anat",
            session="post",
            space="T1w",
            suffix="ct.nii.gz",
            **inputs["post_ct"].wildcards,
        ),
    script:
        "../scripts/registration.py"

if not config.get("SMRIPREP_DIR"):

    rule get_t1w_to_mni_matrix:
        input:
            t1_img=bids(
                root=config["output_dir"],
                suffix="T1w",
                desc="n4",
                datatype="anat",
                session="pre",
                extension=".nii.gz",
                **inputs["pre_t1w"].wildcards,
            ),
            mni_template=target_t1w
        output:
            xfm_slicer=bids(
                root=config["output_dir"],
                datatype="atlas",
                desc="from_T1w-to-MNI",
                suffix="slicer.mat",
                **inputs["pre_t1w"].wildcards,
            ),
            xfm_ras=bids(
                root=config["output_dir"],
                datatype="atlas",
                desc="from_T1w-to-MNI",
                suffix="xfm.txt",
                **inputs["pre_t1w"].wildcards,
            ),
            out_im=bids(
                root=config["output_dir"],
                datatype="atlas",
                space="MNI",
                suffix="T1w.nii.gz",
                **inputs["pre_t1w"].wildcards,
            ),
        script:
            "../scripts/affine.py"

if config["manual_reg_matrix"]:

    rule register_ct:
        input:
            in_im=bids(
                root=config["bids_dir"],
                suffix="ct",
                datatype="ct",
                session="post",
                acq="Electrode",
                extension=".nii.gz",
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
                datatype="anat",
                session="post",
                desc="user_registration",
                space="T1w",
                suffix="ct.nii.gz",
                **inputs["post_ct"].wildcards,
            ),
        script:
            "../scripts/apply_registration.py"

 