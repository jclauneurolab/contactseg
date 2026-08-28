"""Sample tissue probability maps at every contact.

Used when the tissue classes come from a probabilistic segmentation (sMRIPrep,
or the ANTs Atropos run inside this workflow) rather than from freesurfer label
indices. Probabilities are read with trilinear interpolation, so a contact
sitting between voxels gets a blend rather than the value of whichever voxel it
happened to round into.
"""

import nibabel as nib
import numpy as np
import pandas as pd

FCSV_HEADER_ROWS = 3


def sample_volume(img, coords):
    """Trilinearly sample ``img`` at world coordinates ``coords``."""
    data = np.asarray(img.dataobj, dtype=float)
    ijk = nib.affines.apply_affine(np.linalg.inv(img.affine), coords)

    low = np.floor(ijk).astype(int)
    frac = ijk - low

    values = np.zeros(len(coords))
    for corner in np.ndindex(2, 2, 2):
        index = low + np.array(corner)
        weight = np.prod(np.where(np.array(corner) == 1, frac, 1.0 - frac), axis=1)
        inside = np.all((index >= 0) & (index < np.array(data.shape)), axis=1)
        if not inside.any():
            continue
        sampled = np.zeros(len(coords))
        valid = index[inside]
        sampled[inside] = data[valid[:, 0], valid[:, 1], valid[:, 2]]
        values += weight * sampled

    return values


def lookup_tissue_labels(coords_fcsv, probseg, tissue_labels, output_tsv):
    """
    Function that reports tissue probabilities at each contact.

    Parameters
    ----------
    coords_fcsv : str
        Path to the contact coordinates, in the space of the probability maps.
    probseg : list of str
        Probability maps, one per entry of ``tissue_labels`` and in the same
        order.
    tissue_labels : list of str
        Tissue class names, e.g. ``["GM", "WM", "CSF"]``.
    output_tsv : str
        Path to save one row per contact.

    Returns
    -------
    pandas.DataFrame
    """

    contacts = pd.read_csv(coords_fcsv, skiprows=FCSV_HEADER_ROWS, header=None)
    coords = contacts[[1, 2, 3]].to_numpy(float)

    tissue = pd.DataFrame({"name": contacts[11].astype(str)})
    for label, probseg_file in zip(tissue_labels, probseg):
        tissue[f"p_{label}"] = sample_volume(nib.load(probseg_file), coords)

    columns = [f"p_{label}" for label in tissue_labels]
    tissue["tissue"] = (
        tissue[columns]
        .idxmax(axis=1)
        .str.replace("p_", "", regex=False)
        .where(tissue[columns].max(axis=1) > 0, "Unknown")
    )

    tissue.to_csv(output_tsv, sep="\t", index=False, float_format="%.3f")

    return tissue


if __name__ == "__main__":
    lookup_tissue_labels(
        coords_fcsv=snakemake.input.coords,
        probseg=snakemake.input.probseg,
        tissue_labels=snakemake.params.tissue_labels,
        output_tsv=snakemake.output.tissue,
    )
