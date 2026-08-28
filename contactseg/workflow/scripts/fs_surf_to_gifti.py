"""Convert a freesurfer surface to GIFTI, in the scanner RAS of the T1w.

Freesurfer stores surface coordinates in tkrRAS, which is offset from the
scanner RAS the contacts and the T1w live in by the c_ras of the conformed
volume. Skipping that shift is the classic way to end up with surfaces that sit
a centimetre or so away from the volume they were derived from, so the offset
is read from the reference volume header and applied here.
"""

import nibabel as nib
import numpy as np


def get_tkr_to_scanner(ref_vol):
    """Return the tkrRAS to scanner RAS transform of ``ref_vol``."""
    header = nib.load(ref_vol).header

    return header.get_vox2ras() @ np.linalg.inv(header.get_vox2ras_tkr())


def fs_surf_to_gifti(surf, ref_vol, output_gii, structure, apply_cras):
    """
    Function that writes a freesurfer surface as a GIFTI surface file.

    Parameters
    ----------
    surf : str
        Path to the freesurfer surface (e.g. ``lh.white``).
    ref_vol : str
        Path to the conformed volume carrying the c_ras offset.
    output_gii : str
        Path to save the GIFTI surface.
    structure : str
        Workbench structure name, e.g. ``CORTEX_LEFT``.
    apply_cras : bool
        Whether to shift the coordinates into scanner RAS. Spheres are left in
        their own frame.

    Returns
    -------
    None
    """

    vertices, faces = nib.freesurfer.read_geometry(surf)
    vertices = vertices.astype(np.float32)

    if apply_cras:
        xfm = get_tkr_to_scanner(ref_vol)
        vertices = nib.affines.apply_affine(xfm, vertices).astype(np.float32)

    meta = nib.gifti.GiftiMetaData()
    meta["AnatomicalStructurePrimary"] = structure
    meta["GeometricType"] = "Anatomical" if apply_cras else "Spherical"

    gii = nib.gifti.GiftiImage(meta=meta)
    gii.add_gifti_data_array(
        nib.gifti.GiftiDataArray(
            data=vertices,
            intent="NIFTI_INTENT_POINTSET",
            datatype="NIFTI_TYPE_FLOAT32",
            coordsys=nib.gifti.GiftiCoordSystem(dataspace=1, xformspace=1),
        )
    )
    gii.add_gifti_data_array(
        nib.gifti.GiftiDataArray(
            data=faces.astype(np.int32),
            intent="NIFTI_INTENT_TRIANGLE",
            datatype="NIFTI_TYPE_INT32",
        )
    )
    nib.save(gii, output_gii)


if __name__ == "__main__":
    fs_surf_to_gifti(
        surf=snakemake.input.surf,
        ref_vol=snakemake.input.ref_vol,
        output_gii=snakemake.output.surf_gii,
        structure=snakemake.params.structure,
        apply_cras=snakemake.params.apply_cras,
    )
