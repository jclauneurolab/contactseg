"""Map a surface parcellation into the volume, one hemisphere at a time.

Surface atlases (the Yale Brain Atlas, freesurfer's aparc) only exist on the
cortical ribbon, so they are painted into the volume with a ribbon-constrained
mapping between the white and pial surfaces and then merged across
hemispheres. The result is an ordinary segmentation in subject T1w space, which
the contact lookup then reads like any other atlas.
"""

import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np


def label_to_volume(label_gii, midthickness, white, pial, ref_vol, output_nii):
    """Run the workbench ribbon-constrained label to volume mapping."""
    subprocess.run(
        [
            "wb_command",
            "-label-to-volume-mapping",
            label_gii,
            midthickness,
            ref_vol,
            output_nii,
            "-ribbon-constrained",
            white,
            pial,
        ],
        check=True,
    )


def merge_hemispheres(hemi_niis, key_offsets, output_nii):
    """Combine per-hemisphere label volumes into one segmentation.

    Surface parcellations usually number their parcels per hemisphere, so both
    hemispheres arrive with the same keys. ``key_offsets`` shifts each
    hemisphere into its own range before merging -- with 1000 and 2000 the
    result carries freesurfer's own ctx-lh/ctx-rh indices. Labels already
    unique across hemispheres take an offset of zero.

    Where the two hemispheres claim the same voxel the first one wins, so the
    (rare) overlap along the midline resolves deterministically rather than by
    file order.
    """
    merged = None
    affine = None
    for hemi_nii, offset in zip(hemi_niis, key_offsets):
        img = nib.load(hemi_nii)
        data = np.rint(np.asarray(img.dataobj)).astype(np.int32)
        data = np.where(data > 0, data + offset, 0)
        if merged is None:
            merged, affine = data, img.affine
        else:
            merged = np.where(merged > 0, merged, data)

    out = nib.Nifti1Image(merged, affine)
    out.set_qform(affine, code=1)
    out.set_sform(affine, code=1)
    nib.save(out, output_nii)


def surface_atlas_to_volume(
    label_gii, midthickness, white, pial, ref_vol, key_offsets, output_nii
):
    """
    Function that maps a surface atlas into subject volume space.

    Parameters
    ----------
    label_gii : list of str
        Per-hemisphere GIFTI label files.
    midthickness, white, pial : list of str
        Per-hemisphere surfaces, in the same order as ``label_gii``.
    ref_vol : str
        Reference volume defining the output grid.
    key_offsets : list of int
        Value added to the labels of each hemisphere before merging.
    output_nii : str
        Path to save the merged segmentation.

    Returns
    -------
    None
    """

    with tempfile.TemporaryDirectory() as tmpdir:
        hemi_niis = []
        for i, label in enumerate(label_gii):
            hemi_nii = str(Path(tmpdir) / f"hemi-{i}_dseg.nii.gz")
            label_to_volume(
                label,
                midthickness[i],
                white[i],
                pial[i],
                ref_vol,
                hemi_nii,
            )
            hemi_niis.append(hemi_nii)

        merge_hemispheres(hemi_niis, key_offsets, output_nii)


if __name__ == "__main__":
    surface_atlas_to_volume(
        label_gii=snakemake.input.label_gii,
        midthickness=snakemake.input.midthickness,
        white=snakemake.input.white,
        pial=snakemake.input.pial,
        ref_vol=snakemake.input.ref_vol,
        key_offsets=snakemake.params.key_offsets,
        output_nii=snakemake.output.atlas_dseg,
    )
