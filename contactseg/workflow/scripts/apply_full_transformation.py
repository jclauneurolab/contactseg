import pandas as pd
import ants


def warp_contacts_to_mni(input_fcsv, output_fcsv, transforms):
    """
    Function that applies a full transformation to warp contact coordinates to MNI space.

    Parameters
    ----------
    input_fcsv : str
        Path to the input .fcsv file containing native space contact coordinates.
    output_fcsv : str
        Path to save the transformed .fcsv file with updated template space contact coordinates.
    affine : str
        Path to the affine transformation matrix (4x4) from T1w to MNI space.
    affine_syn : str
        Path to the affine transformation file generated during the SyN registration.
    forward_warp : str
        Path to the forward warp field (deformation field from T1w to MNI space).

    Returns
    -------
    None
    """

    contacts_df = pd.read_csv(input_fcsv, skiprows=3, header=None)
    points = contacts_df[[1, 2, 3]].copy()
    points.columns = ["x", "y", "z"]

    # Convert RAS to LPS
    points["x"] *= -1
    points["y"] *= -1

    transforms = snakemake.input.transforms
    if isinstance(transforms, str):
        transforms = [transforms]

    if len(transforms) == 3:
        inv_flags = [False, True, True]
    else:
        inv_flags = [False]

    # Apply Transform- note order is affine->affine_syn->forward_warp
    transformed_points = ants.apply_transforms_to_points(
        dim=3,
        points=points,
        transformlist=transforms,
        whichtoinvert=inv_flags,
    )
    transformed_points.columns = ["x", "y", "z"]

    # Convert LPS back to RAS
    transformed_points["x"] *= -1
    transformed_points["y"] *= -1

    # Save
    contacts_df[1] = transformed_points["x"].values
    contacts_df[2] = transformed_points["y"].values
    contacts_df[3] = transformed_points["z"].values

    # Write output with original Slicer header
    with open(input_fcsv, "r") as f:
        header = [line for line in f if line.startswith("#")]

    with open(output_fcsv, "w") as f:
        f.writelines(header)
        # mode='a' appends the data after the header
        contacts_df.to_csv(f, index=False, header=False, mode="a", float_format="%.3f")


if __name__ == "__main__":
    warp_contacts_to_mni(
        input_fcsv=snakemake.input.coords,
        transforms=snakemake.input.transforms,
        output_fcsv=snakemake.output.output_coords,
    )
