"""Convert a freesurfer annotation to a GIFTI label file.

Keeps the colour table that ships with the annotation, so the parcellation
looks the same in workbench as it does in freeview.
"""

import nibabel as nib
import numpy as np


def annot_to_label_gii(annot, output_gii, structure):
    """
    Function that writes a freesurfer annotation as a GIFTI label file.

    Parameters
    ----------
    annot : str
        Path to the freesurfer annotation (e.g. ``lh.aparc.DKTatlas.annot``).
    output_gii : str
        Path to save the GIFTI label file.
    structure : str
        Workbench structure name, e.g. ``CORTEX_LEFT``.

    Returns
    -------
    None
    """

    labels, ctab, names = nib.freesurfer.read_annot(annot)
    names = [name.decode() if isinstance(name, bytes) else name for name in names]

    # read_annot returns -1 for unlabelled vertices
    labels = np.where(labels < 0, 0, labels).astype(np.int32)

    label_table = nib.gifti.GiftiLabelTable()
    for index, name in enumerate(names):
        red, green, blue, _, _ = ctab[index][:5]
        label = nib.gifti.GiftiLabel(
            key=index,
            red=red / 255.0,
            green=green / 255.0,
            blue=blue / 255.0,
            alpha=1.0,
        )
        label.label = name
        label_table.labels.append(label)

    meta = nib.gifti.GiftiMetaData()
    meta["AnatomicalStructurePrimary"] = structure

    gii = nib.gifti.GiftiImage(meta=meta, labeltable=label_table)
    gii.add_gifti_data_array(
        nib.gifti.GiftiDataArray(
            data=labels,
            intent="NIFTI_INTENT_LABEL",
            datatype="NIFTI_TYPE_INT32",
        )
    )
    nib.save(gii, output_gii)


if __name__ == "__main__":
    annot_to_label_gii(
        annot=snakemake.input.annot,
        output_gii=snakemake.output.label_gii,
        structure=snakemake.params.structure,
    )
