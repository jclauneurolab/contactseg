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
                    suffix="atlas_labelled_contactseg.fcsv",
                    datatype="atlas",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=config["output_dir"],
                    datatype="atlas",
                    suffix="mni_transformed_contactseg.fcsv",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=config["output_dir"],
                    datatype="atlas",
                    suffix="atlas_in_t1_space.nii.gz",
                    **inputs["post_ct"].wildcards,
                )
            )
        )
    return final


def get_template_atlas(wildcards=None):

    template_name = config.get("template_flow") or config.get("template-flow", False)
    atlas_name = config.get("template_atlas") or config.get("template-atlas", False)
    if template_name and atlas_name:
        return str(Path(config["output_dir"]) / f"templateflow/tpl-{template_name}_atlas-{atlas_name}_dseg.nii.gz")
    return str(Path(workflow.basedir).parent.parent / "resources/atlases/tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii")

def get_template_atlas_tsv(wildcards=None):

    template_name = config.get("template_flow") or config.get("template-flow", False)
    atlas_name = config.get("template_atlas") or config.get("template-atlas", False)
    if template_name and atlas_name:
        return str(Path(config["output_dir"]) / f"templateflow/tpl-{template_name}_atlas-{atlas_name}_dseg.tsv")
    return str(Path(workflow.basedir).parent.parent / "resources/atlases/tpl-MNI152NLin2009cSym_atlas-CerebA_dseg.tsv")

   
def get_template_t1w(wildcards=None):

    template_name = config.get("template_flow") or config.get("template-flow", False)
    atlas_name = config.get("template_atlas") or config.get("template-atlas", False)
    """
    Returns the path to the T1w template image. 
    Uses the downloaded TemplateFlow image if specified, otherwise uses a local fallback.
    """
    
    if template_name:
        return str(Path(config["output_dir"]) / f"templateflow/tpl-{template_name}_T1w.nii.gz")
    
    # Fallback to your local default template if template-flow is not used
    return str(Path(workflow.basedir).parent.parent / "resources/templates/tpl-MNI152NLin2009cSym_res-1_T1w.nii.gz")