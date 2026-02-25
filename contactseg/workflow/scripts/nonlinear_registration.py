import ants 

def nonlinear_registration(t1_image, mni_template, forward_warp, inverse_warp, affine):

    fixed = ants.image_read(t1_image)
    moving = ants.image_read(mni_template)

    reg = ants.registration(
        fixed=fixed, 
        moving=moving,
        type_of_transform="SyN" 
    )

    ants.write_transform(reg["fwdtransforms"][0], forward_warp)
    ants.write_transform(reg["invtransforms"][0], inverse_warp)
    ants.write_transform(reg["fwdtransforms"][1], affine)

if __name__ == "__main__":
    nonlinear_registration(
        t1_image=snakemake.input.t1_image,
        mni_template=snakemake.input.mni_template,
        forward_warp=snakemake.output.forward_warp,
        inverse_warp=snakemake.output.inverse_warp,
        affine=snakemake.output.affine,
    )
