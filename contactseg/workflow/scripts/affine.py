import ants
import csv 
import numpy as np


def antsmat2mat(afftransform, m_center):
    """
    Function that creates a transformation matrix
    from ANTs .mat output.
    Note, transformation matrix is in LPS format.

    Parameters
    ----------
    afftransform : numpy array
        parameters portion of output transformation
    m_center : numpy array
        fixed parameters portion of output transformation

    Returns
    -------
    mat : nd.array
        4x4 transformation matrix
    """

    # Reshaping the first 9 elements of afftransform
    # into a 3x3 matrix and adding the translation vector
    mat = np.hstack(
        (
            np.reshape(afftransform[:9], (3, 3)),
            np.array(afftransform[9:12]).reshape(3, 1),
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


def affine_registration(t1_img, mni_template, out_im, xfm_ras, xfm_slicer):
    """
    Function that performs affine registration.

    Parameters
    ----------
    t1_img : str
        Path to the T1-weighted image.
    mni_template : str
        Path to the MNI template image.
    out_im : str
        Path to save the resampled (registered) image.
    xfm_ras : str
        Path to save the output affine transformation matrix (4x4).
    xfm_slicer : str
        Path to save the ANTs affine transform file.

    Returns
    -------
    None
    """

    #Load images
    t1_img = ants.image_read(t1_img)
    mni_template = ants.image_read(mni_template)

    #Perform affine registration
    registration_result = ants.registration(
        fixed=mni_template,
        moving=t1_img,
        type_of_transform="AffineFast",
        grad_step=0.25,
        aff_iterations=[1000, 500, 250, 100],
        aff_sampling=32,
        affine_random_sampling_rate=0.25,
        aff_shrink_factors=[8, 4, 2, 1],
        aff_smoothing_sigmas=[3, 2, 1, 0],
    )

    # Get the registered (warped) moving image
    registered_image = registration_result["warpedmovout"]

    # Get the forward transformation
    transformation_file_path = registration_result["fwdtransforms"][0]

    # Load the transformation matrix directly
    transform = ants.read_transform(transformation_file_path)
    full_matrix = antsmat2mat(transform.parameters, transform.fixed_parameters)

    # Save the registered image
    ants.image_write(registered_image, out_im)
    ants.write_transform(transform, xfm_slicer)

    # Save the 4x4 transformation matrix to a file
    with open(xfm_ras, "w", newline="") as file:
        writer = csv.writer(file, delimiter=" ")
        for row in full_matrix:
            writer.writerow(row)

if __name__ == "__main__":
    affine_registration(
        t1_img=snakemake.input.t1_img,
        mni_template=snakemake.input.mni_template,
        xfm_ras=snakemake.output.xfm_ras,
        xfm_slicer=snakemake.output.xfm_slicer,
        out_im=snakemake.output.out_im,
    )