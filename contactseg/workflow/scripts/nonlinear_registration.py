import ants
import shutil

def nonlinear_registration(t1_img, mni_template, forward_warp, inverse_warp, affine):
    fixed = ants.image_read(mni_template)
    moving = ants.image_read(t1_img)

    reg = ants.registration(
        fixed=fixed,
        moving=moving,
        type_of_transform="SyN"
    )

    # reg["fwdtransforms"] and reg["invtransforms"] are lists of file paths
    shutil.copy(reg["fwdtransforms"][0], forward_warp)
    shutil.copy(reg["invtransforms"][0], inverse_warp)
    shutil.copy(reg["fwdtransforms"][1], affine)

if __name__ == "__main__":
    nonlinear_registration(
        t1_img=snakemake.input.t1_img,
        mni_template=snakemake.input.mni_template,
        forward_warp=snakemake.output.forward_warp,
        inverse_warp=snakemake.output.inverse_warp,
        affine=snakemake.output.affine,
    )