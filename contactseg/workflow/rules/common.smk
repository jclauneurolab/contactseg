from pathlib import Path

deriv_root = str(Path(config["output_dir"]) / "derivatives" / "contactseg")

def get_reg_matrix():
    if not config["manual_reg_matrix"]:
        return bids(
            root=deriv_root,
            datatype="anat",
            desc="from_ct-to_T1w",
            suffix="xfm",
            extension=".txt",
            **inputs["post_ct"].wildcards,
        )
    else:
        return bids(
            root=config["bids_dir"],
            suffix="xfm",
            extension=".txt",
            **inputs["post_ct"].wildcards,
        )


def get_registered_ct_image():
    if not config["manual_reg_matrix"]:
        return (
            bids(
                root=deriv_root,
                datatype="anat",
                session="post",
                space="T1w",
                suffix="ct.nii.gz",
                **inputs["post_ct"].wildcards,
            ),
        )

    else:
        return bids(
            root=deriv_root,
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
                root=deriv_root,
                suffix="dseg.nii.gz",
                session="post",
                desc="contacts_nnUNet",
                **inputs["post_ct"].wildcards,
            )
        )
    )
    if config["label"]:
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=deriv_root,
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
                    root=deriv_root,
                    datatype="ieeg",
                    space="T1w",
                    suffix="coords",
                    desc="transformed",
                    session="post",
                    extension=".fcsv",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=deriv_root,
                    space="T1w",
                    desc="contacts_nnUNet",
                    session="post",
                    suffix="dseg.nii.gz",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
    if config["contacts_qc"]:
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=deriv_root,
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
                    root=deriv_root,
                    datatype="ieeg",
                    space="T1w",
                    suffix="coords",
                    desc="anatomical_labels",
                    session="post",
                    extension=".fcsv",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=deriv_root,
                    datatype="ieeg",
                    suffix="coords",
                    extension=".fcsv",
                    space="TEMPLATE",
                    session="post",
                    desc="labeled",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=deriv_root,
                    datatype="anat",
                    desc="ATLAS",
                    space="T1w",
                    suffix="dseg",
                    extension=".nii.gz",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
    return final
