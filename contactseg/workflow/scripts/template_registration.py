"""Register the subject T1w to a volumetric template with ANTsPy.

A single SyN registration is run with the template as the fixed image. Its
pieces are written out separately so that the rest of the workflow can compose
them in either direction:

    T1w -> template (images) : [warp, affine]
    template -> T1w (images) : [affine (inverted), invwarp]

Point sets travel the opposite way to images under the ANTs convention, so
contacts are pushed into template space with the second chain.
"""

import csv

import ants
import numpy as np


def antsmat2mat(transform, m_center):
    """
    Function that creates a transformation matrix
    from ANTs .mat output.
    Note, transformation matrix is in LPS format.

    Parameters
    ----------

    transform : numpy array
        parameters portion of output transformation
    m_center : numpy array
        fixted parameters portion of output transformation

    Returns
    -------
    mat : nd.array
        4x4 transformation matrix
    """

    # Reshaping the first 9 elements of afftransform
    # into a 3x3 matrix and adding the translation vector
    mat = np.hstack(
        (
            np.reshape(transform[:9], (3, 3)),
            np.array(transform[9:12]).reshape(3, 1),
        )
    )

    # Adding the last row to the matrix
    mat = np.vstack((mat, [0, 0, 0, 1]))

    # Calculating the offset
    m_translation = mat[:3, 3]
    m_offset = np.zeros(3)

    for i in range(3):
        m_offset[i] = m_translation[i] + m_center[i]
        for j in range(3):
            m_offset[i] -= mat[i, j] * m_center[j]

    # Updating the translation part of the matrix with the calculated offset
    mat[:3, 3] = m_offset

    d = np.array([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

    lps_inmatrix = np.linalg.inv(mat)

    ras_inmatrix = d @ lps_inmatrix @ d

    return ras_inmatrix


def template_registration(
    t1w, template, xfm_affine, warp, invwarp, xfm_ras, warped_t1w
):
    """
    Function that registers a subject T1w image to a template.

    Parameters
    ----------
    t1w : str
        Path to the subject T1w image (moving image).
    template : str
        Path to the template image (fixed image).
    xfm_affine : str
        Path to save the affine part of the registration (ANTs .mat).
    warp : str
        Path to save the forward deformation field (T1w to template).
    invwarp : str
        Path to save the inverse deformation field (template to T1w).
    xfm_ras : str
        Path to save the affine as a 4x4 RAS matrix, readable by 3D Slicer.
    warped_t1w : str
        Path to save the subject T1w resampled into template space.

    Returns
    -------
    None
    """

    fixed_image = ants.image_read(template)
    moving_image = ants.image_read(t1w)

    registration_result = ants.registration(
        fixed=fixed_image,
        moving=moving_image,
        type_of_transform="SyN",
        aff_metric="mattes",
        syn_metric="mattes",
    )

    # fwdtransforms is [warp, affine], invtransforms is [affine, invwarp]
    ants.image_write(ants.image_read(registration_result["fwdtransforms"][0]), warp)
    ants.image_write(ants.image_read(registration_result["invtransforms"][1]), invwarp)
    ants.image_write(registration_result["warpedmovout"], warped_t1w)

    transform = ants.read_transform(registration_result["fwdtransforms"][1])
    ants.write_transform(transform, xfm_affine)

    full_matrix = antsmat2mat(transform.parameters, transform.fixed_parameters)
    with open(xfm_ras, "w", newline="") as file:
        writer = csv.writer(file, delimiter=" ")
        for row in full_matrix:
            writer.writerow(row)


if __name__ == "__main__":
    template_registration(
        t1w=snakemake.input.t1w,
        template=snakemake.input.template,
        xfm_affine=snakemake.output.xfm_affine,
        warp=snakemake.output.warp,
        invwarp=snakemake.output.invwarp,
        xfm_ras=snakemake.output.xfm_ras,
        warped_t1w=snakemake.output.warped_t1w,
    )
