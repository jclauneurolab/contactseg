import ants


def nonlinear_registration(
    t1_in_mni_img, templateflow_paths, forward_warp, inverse_warp, affine_syn
):
    """
    Function that performs nonlinear registration using the SyN transformation.

    Parameters
    ----------
    t1_in_mni_img : str
        Path to the T1-weighted image in MNI space (moving image).
    mni_template : str
        Path to the MNI template image (fixed image).
    forward_warp : str
        Path to save the forward warp field (deformation field from T1 to MNI space).
    inverse_warp : str
        Path to save the inverse warp field (deformation field from MNI to T1 space).
    affine_syn : str
        Path to save the affine transformation file generated during the SyN registration.

    Returns
    -------
    None
    """
    # Read the file paths from the template_txt
    with open(templateflow_paths, "r") as f:
        lines = f.readlines()

    # Extract path for template
    template_path = lines[0].strip().split(":")[1].strip()

    fixed = ants.image_read(template_path)
    moving = ants.image_read(t1_in_mni_img)

    reg = ants.registration(fixed=fixed, moving=moving, type_of_transform="SyN")

    fwd_warp = ants.image_read(reg["fwdtransforms"][0])
    ants.image_write(fwd_warp, forward_warp)

    inv_warp = ants.image_read(reg["invtransforms"][1])
    ants.image_write(inv_warp, inverse_warp)

    affine_mat = ants.read_transform(reg["fwdtransforms"][1])
    ants.write_transform(affine_mat, affine_syn)


if __name__ == "__main__":
    nonlinear_registration(
        t1_in_mni_img=snakemake.input.t1_in_mni_img,
        templateflow_paths=snakemake.input.templateflow_paths,
        forward_warp=snakemake.output.forward_warp,
        inverse_warp=snakemake.output.inverse_warp,
        affine_syn=snakemake.output.affine_syn,
    )
