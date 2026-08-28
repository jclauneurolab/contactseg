"""Import a freesurfer segmentation without touching its labels.

The conformed volume freesurfer writes shares its scanner RAS frame with the
T1w the contacts were registered to, so the segmentation only has to change
container, not grid. Resampling is deliberately avoided: nearest-neighbour
interpolation of a categorical volume moves boundaries around for no gain,
since the lookup happens in world coordinates anyway.
"""

import nibabel as nib
import numpy as np


def import_fs_atlas(atlas_dseg, output_dseg):
    """
    Function that rewrites a freesurfer segmentation as a NIfTI file.

    Parameters
    ----------
    atlas_dseg : str
        Path to the freesurfer segmentation (.mgz or .nii.gz).
    output_dseg : str
        Path to save the segmentation as NIfTI, on its original grid.

    Returns
    -------
    None
    """

    seg = nib.load(atlas_dseg)
    data = np.rint(np.asarray(seg.dataobj)).astype(np.int32)

    out = nib.Nifti1Image(data, seg.affine)
    out.set_qform(seg.affine, code=1)
    out.set_sform(seg.affine, code=1)
    nib.save(out, output_dseg)


if __name__ == "__main__":
    import_fs_atlas(
        atlas_dseg=snakemake.input.atlas_dseg,
        output_dseg=snakemake.output.atlas_dseg,
    )
