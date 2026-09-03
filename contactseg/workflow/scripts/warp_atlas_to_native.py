"""Pull a template-space atlas segmentation into subject T1w space.

When the transforms came from a derivatives dataset computed on a different
acquisition, the ANTs chain lands the atlas in that image's frame. The bridging
rigid is then folded into the output affine rather than applied as a second
resampling, so the labels are interpolated once and only once.
"""

import ants
import nibabel as nib
import numpy as np


def warp_atlas_to_native(
    atlas_dseg, t1w, transforms, invert_flags, output_dseg, xfm_ras=None
):
    """
    Function that resamples an atlas segmentation into subject space.

    Parameters
    ----------
    atlas_dseg : str
        Path to the atlas segmentation, in template space.
    t1w : str
        Path to the reference image defining the output grid: the contactseg
        T1w, or the derivatives anatomical when a bridge is in use.
    transforms : list of str
        ANTs transform chain resampling template images into subject space.
    invert_flags : list of bool
        Which entries of ``transforms`` to invert.
    output_dseg : str
        Path to save the atlas segmentation in subject space.
    xfm_ras : str, optional
        Path to a 4x4 RAS matrix taking the derivatives anatomical to the
        contactseg T1w, composed into the output affine when given.

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

    if not xfm_ras:
        ants.image_write(atlas_in_native, output_dseg)
        return

    # fold the bridge into the affine rather than resampling a second time
    ants.image_write(atlas_in_native, output_dseg)
    warped = nib.load(output_dseg)
    affine = np.loadtxt(xfm_ras) @ warped.affine
    bridged = nib.Nifti1Image(np.rint(warped.get_fdata()).astype(np.int32), affine)
    bridged.set_qform(affine, code=1)
    bridged.set_sform(affine, code=1)
    nib.save(bridged, output_dseg)


def first_or_none(value):
    """Return a single path from an input that may be an empty list."""
    if isinstance(value, (list, tuple)):
        return value[0] if value else None

    return value or None


if __name__ == "__main__":
    warp_atlas_to_native(
        atlas_dseg=snakemake.input.atlas_dseg,
        t1w=snakemake.input.t1w,
        transforms=snakemake.input.transforms,
        invert_flags=snakemake.params.invert_flags,
        output_dseg=snakemake.output.atlas_dseg,
        xfm_ras=first_or_none(snakemake.input.xfm_ras),
    )
