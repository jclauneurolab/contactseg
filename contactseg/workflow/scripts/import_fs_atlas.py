"""Import a freesurfer segmentation without touching its labels.

The conformed volume freesurfer writes shares its scanner RAS frame with the
image freesurfer was run on, so the segmentation only has to change container,
not grid. Resampling is deliberately avoided: nearest-neighbour interpolation
of a categorical volume moves boundaries around for no gain, since the lookup
happens in world coordinates anyway.

When freesurfer was run on a different acquisition than the contactseg T1w, the
bridging transform is folded into the image affine rather than applied to the
voxels. A rigid transform composes exactly with an affine, so the labels are
still never interpolated -- the volume simply sits on a rotated grid, which the
lookup handles because it works in world coordinates.
"""

import nibabel as nib
import numpy as np


def import_fs_atlas(atlas_dseg, output_dseg, xfm_ras=None):
    """
    Function that rewrites a freesurfer segmentation as a NIfTI file.

    Parameters
    ----------
    atlas_dseg : str
        Path to the freesurfer segmentation (.mgz or .nii.gz).
    output_dseg : str
        Path to save the segmentation as NIfTI, on its original grid.
    xfm_ras : str, optional
        Path to a 4x4 RAS matrix taking the freesurfer anatomical to the
        contactseg T1w. Composed into the affine when given.

    Returns
    -------
    None
    """

    seg = nib.load(atlas_dseg)
    data = np.rint(np.asarray(seg.dataobj)).astype(np.int32)

    affine = seg.affine
    if xfm_ras:
        affine = np.loadtxt(xfm_ras) @ affine

    out = nib.Nifti1Image(data, affine)
    out.set_qform(affine, code=1)
    out.set_sform(affine, code=1)
    nib.save(out, output_dseg)


def first_or_none(value):
    """Return a single path from an input that may be an empty list."""
    if isinstance(value, (list, tuple)):
        return value[0] if value else None

    return value or None


if __name__ == "__main__":
    import_fs_atlas(
        atlas_dseg=snakemake.input.atlas_dseg,
        output_dseg=snakemake.output.atlas_dseg,
        xfm_ras=first_or_none(snakemake.input.xfm_ras),
    )
