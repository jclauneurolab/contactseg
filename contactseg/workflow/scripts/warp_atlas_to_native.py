"""Pull a template-space atlas segmentation into subject T1w space."""

import ants


def warp_atlas_to_native(atlas_dseg, t1w, transforms, invert_flags, output_dseg):
    """
    Function that resamples an atlas segmentation into subject space.

    Parameters
    ----------
    atlas_dseg : str
        Path to the atlas segmentation, in template space.
    t1w : str
        Path to the subject T1w image, used as the reference grid.
    transforms : list of str
        ANTs transform chain resampling template images into subject space.
    invert_flags : list of bool
        Which entries of ``transforms`` to invert.
    output_dseg : str
        Path to save the atlas segmentation in subject space.

    Returns
    -------
    None
    """

    atlas = ants.image_read(atlas_dseg)
    reference = ants.image_read(t1w)

    atlas_in_native = ants.apply_transforms(
        fixed=reference,
        moving=atlas,
        transformlist=list(transforms),
        whichtoinvert=list(invert_flags),
        interpolator="nearestNeighbor",
    )

    ants.image_write(atlas_in_native, output_dseg)


if __name__ == "__main__":
    warp_atlas_to_native(
        atlas_dseg=snakemake.input.atlas_dseg,
        t1w=snakemake.input.t1w,
        transforms=snakemake.input.transforms,
        invert_flags=snakemake.params.invert_flags,
        output_dseg=snakemake.output.atlas_dseg,
    )
