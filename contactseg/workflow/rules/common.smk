def get_reg_matrix():
    if not config["manual_reg_matrix"]:
        return bids(
            root=config["output_dir"],
            datatype="registration",
            desc="from_ct-to_T1w",
            suffix="xfm.txt",
            **inputs["post_ct"].wildcards,
        )
    else:
        return bids(
            root=config["bids_dir"],
            suffix="xfm.txt",
            **inputs["post_ct"].wildcards,
        )


def get_registered_ct_image():
    if not config["manual_reg_matrix"]:
        return (
            bids(
                root=config["output_dir"],
                datatype="anat",
                session="post",
                space="T1w",
                suffix="ct.nii.gz",
                **inputs["post_ct"].wildcards,
            ),
        )

    else:
        return bids(
            root=config["output_dir"],
            datatype="anat",
            session="post",
            space="T1w",
            desc="user_registration",
            suffix="ct.nii.gz",
            **inputs["post_ct"].wildcards,
        )


def get_final_output():
    final = []
    final.extend(
        inputs["post_ct"].expand(
            bids(
                root=config["output_dir"],
                suffix="dseg.nii.gz",
                desc="contacts_nnUNet",
                **inputs["post_ct"].wildcards,
            )
        )
    )
    if config["label"]:
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=config["output_dir"],
                    datatype="ieeg",
                    space="T1w",
                    suffix="electrodes",
                    session="post",
                    extension=".tsv",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
    if config["transform"]:
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=config["output_dir"],
                    datatype="slicer_fcsv",
                    suffix="transformed_contactseg.fcsv",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=config["output_dir"],
                    space="T1w",
                    desc="contacts_nnUNet",
                    suffix="dseg.nii.gz",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
    if config["contacts_qc"]:
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=config["output_dir"],
                    datatype="qc",
                    desc="contactseg",
                    suffix="qc.html",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
    if config["atlas_labels"]:
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=config["output_dir"],
                    suffix="labelled_contactseg_mni.fcsv",
                    datatype="slicer_fcsv",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
    if config["affine_test"]:
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=config["output_dir"],
                    datatype="slicer_fcsv",
                    suffix="mni_transformed_contactseg.fcsv",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=config["output_dir"],
                    datatype="anat",
                    suffix="atlas_in_t1_space.nii.gz",
                    session="post",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
    return final
