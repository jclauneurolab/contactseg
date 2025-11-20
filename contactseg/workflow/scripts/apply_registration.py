import ants
import numpy as np


def apply_registration(in_img, ref_img, transform_matrix, out_img):

    img = ants.image_read(in_img)
    ref = ants.image_read(ref_img)
    mat = np.loadtxt(transform_matrix)

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
    )
