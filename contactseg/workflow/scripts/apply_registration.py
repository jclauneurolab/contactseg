import ants
import numpy as np
import nibabel as nib


def apply_registration(in_img, ref_img, transform_matrix, out_img, non_interpolated=False):

    #Load transformation matrix
    mat = np.loadtxt(transform_matrix)

    if non_interpolated:
        float_img = nib.load(in_img)
        new_affine = np.dot(mat, float_img.affine)

        float_img_trans = nib.Nifti1Image(float_img.get_fdata(),
                                          affine=new_affine,
                                          header=float_img.header,)
        float_img_trans.set_qform(float_img_trans.affine, code=1)

        nib.save(float_img_trans, out_img)
        return

        
    img = ants.image_read(in_img)
    ref = ants.image_read(ref_img)

    lps2ras = np.diag([-1, -1, 1, 1])
    ras2lps = np.diag([-1, -1, 1, 1])

    transform_lps = np.dot(ras2lps, np.dot(np.linalg.inv(mat), lps2ras))

    rotation_scale = transform_lps[:3, :3].flatten()
    translation = transform_lps[:3, 3]

    # create the transform from the .txt
    transform = ants.create_ants_transform(
        transform_type="AffineTransform", dimension=3
    )

    # set rotation and scaling
    transform.set_parameters(np.concatenate([rotation_scale, translation]))

    # apply the transform using the same image as reference
    transformed_img = ants.apply_ants_transform_to_image(
        transform=transform, image=img, reference=ref, interpolation="bspline"
    )

    ants.image_write(transformed_img, out_img)


if __name__ == "__main__":
    apply_registration(
        in_img=snakemake.input.in_im,
        ref_img=snakemake.input.ref_im,
        transform_matrix=snakemake.input.transform_matrix,
        out_img=snakemake.output.out_im,
        non_interpolated=bool(snakemake.params.get("non_interpolated", False))
    )
