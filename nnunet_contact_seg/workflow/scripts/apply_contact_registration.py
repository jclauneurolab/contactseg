import ants
import numpy as np


def apply_contact_registration(in_img, transform_matrix, out_img):

    img = ants.image_read(in_img)
    mat = np.loadtxt(transform_matrix)

    rotation_scale = mat[:3, :3].flatten()
    translation = mat[:3, 3]

    # create the transform from the .txt
    transform = ants.create_ants_transform(
        transform_type="AffineTransform", dimension=3
    )

    # set rotation and scaling
    transform.set_parameters(np.concatenate([rotation_scale, translation]))

    # apply the transform using the same image as reference
    transformed_img = ants.apply_ants_transform_to_image(
        transform=transform, image=img, reference=img
    )

    ants.image_write(transformed_img, out_img)


if __name__ == "__main__":
    apply_contact_registration(
        in_img=snakemake.input.in_im,
        transform_matrix=snakemake.input.transform_matrix,
        out_img=snakemake.output.out_im,
    )
