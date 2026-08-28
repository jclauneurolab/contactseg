"""Rigid registration of a derivatives anatomical to the contactseg T1w.

sMRIPrep and FreeSurfer are often run on a different acquisition than the one
contactseg localises contacts against -- a non-contrast T1w for the same
patient, say. Their transforms, segmentations and surfaces then live in that
image's scanner RAS frame, which is offset from the pre_t1w frame by however
much the patient moved between the two scans. This rule measures that offset
once, and the rest of the workflow applies it as a plain 4x4 in RAS.

The matrix is written in the same convention as the CT-to-T1w matrix the
workflow already produces: a 4x4 RAS transform taking points in the moving
(derivatives) image to points in the fixed (pre_t1w) image.
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


def rigid_registration(moving_image, fixed_image, xfm_ras, xfm_slicer, out_im):
    """
    Function that rigidly registers a derivatives anatomical to the T1w.

    Parameters
    ----------
    moving_image : str
        Path to the anatomical the derivatives were computed on.
    fixed_image : str
        Path to the contactseg reference T1w.
    xfm_ras : str
        Path to save the 4x4 RAS matrix, moving to fixed.
    xfm_slicer : str
        Path to save the ANTs transform.
    out_im : str
        Path to save the moving image resampled onto the T1w, for QC.

    Returns
    -------
    None
    """

    fixed = ants.image_read(fixed_image)
    moving = ants.image_read(moving_image)

    registration_result = ants.registration(
        fixed=fixed,
        moving=moving,
        type_of_transform="Rigid",
        aff_metric="mattes",
    )

    transform = ants.read_transform(registration_result["fwdtransforms"][0])
    full_matrix = antsmat2mat(transform.parameters, transform.fixed_parameters)

    ants.image_write(registration_result["warpedmovout"], out_im)
    ants.write_transform(transform, xfm_slicer)

    with open(xfm_ras, "w", newline="") as file:
        writer = csv.writer(file, delimiter=" ")
        for row in full_matrix:
            writer.writerow(row)

    displacement = np.linalg.norm(full_matrix[:3, 3])
    print(
        f"[bridge] {moving_image} -> {fixed_image}: "
        f"translation {displacement:.2f} mm"
    )


if __name__ == "__main__":
    rigid_registration(
        moving_image=snakemake.input.derivatives_anat,
        fixed_image=snakemake.input.t1w,
        xfm_ras=snakemake.output.xfm_ras,
        xfm_slicer=snakemake.output.xfm_slicer,
        out_im=snakemake.output.out_im,
    )
