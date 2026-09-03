"""Warp contact coordinates from subject T1w space into template space.

ANTs transforms map points in the opposite direction to images, so the chain
that resamples a template image into subject space is the one that carries
subject points into template space.

When the transforms were computed by something else on a different acquisition
-- sMRIPrep run on a non-contrast T1w, say -- the contacts are first carried
into that image's frame with a plain 4x4 in RAS, the same way the workflow
already moves coordinates through the CT-to-T1w matrix. Composing the bridge
into the ANTs chain instead would work in principle, but keeping it out means
the nonlinear part stays the two-transform chain that is covered by tests.
"""

import ants
import nibabel as nib
import numpy as np
import pandas as pd

FCSV_HEADER_ROWS = 3


def read_fcsv(input_fcsv):
    """Read the data rows of a Slicer fcsv file."""
    return pd.read_csv(input_fcsv, skiprows=FCSV_HEADER_ROWS, header=None)


def write_fcsv(input_fcsv, output_fcsv, coords_df):
    """Write ``coords_df`` back out under the header of ``input_fcsv``."""
    with open(input_fcsv, "r") as file:
        header = [line for line in file if line.startswith("#")]

    with open(output_fcsv, "w") as file:
        file.writelines(header)
        coords_df.to_csv(file, index=False, header=False, mode="a", float_format="%.3f")


def warp_coords(input_fcsv, output_fcsv, transforms, invert_flags, xfm_ras=None):
    """
    Function that warps fcsv coordinates through an ANTs transform chain.

    Parameters
    ----------
    input_fcsv : str
        Path to the fcsv file holding the coordinates to warp (RAS).
    output_fcsv : str
        Path to save the warped coordinates, keeping the original header.
    transforms : list of str
        ANTs transform chain, ordered as it would be to resample an image in
        the opposite direction.
    invert_flags : list of bool
        Which entries of ``transforms`` to invert.
    xfm_ras : str, optional
        Path to a 4x4 RAS matrix taking the derivatives anatomical to the
        contactseg T1w. Its inverse is applied to the contacts before the ANTs
        chain, so that they enter the chain in the frame it was computed in.

    Returns
    -------
    None
    """

    coords_df = read_fcsv(input_fcsv)
    points = coords_df[[1, 2, 3]].copy()
    points.columns = ["x", "y", "z"]

    if xfm_ras:
        points[["x", "y", "z"]] = nib.affines.apply_affine(
            np.linalg.inv(np.loadtxt(xfm_ras)), points.to_numpy(float)
        )

    # fcsv coordinates are RAS, ANTs works in LPS
    points["x"] *= -1
    points["y"] *= -1

    warped = ants.apply_transforms_to_points(
        dim=3,
        points=points,
        transformlist=list(transforms),
        whichtoinvert=list(invert_flags),
    )
    warped.columns = ["x", "y", "z"]

    warped["x"] *= -1
    warped["y"] *= -1

    coords_df[1] = warped["x"].values
    coords_df[2] = warped["y"].values
    coords_df[3] = warped["z"].values

    write_fcsv(input_fcsv, output_fcsv, coords_df)


def first_or_none(value):
    """Return a single path from an input that may be an empty list."""
    if isinstance(value, (list, tuple)):
        return value[0] if value else None

    return value or None


if __name__ == "__main__":
    warp_coords(
        input_fcsv=snakemake.input.coords,
        output_fcsv=snakemake.output.warped_coords,
        transforms=snakemake.input.transforms,
        invert_flags=snakemake.params.invert_flags,
        xfm_ras=first_or_none(snakemake.input.xfm_ras),
    )
