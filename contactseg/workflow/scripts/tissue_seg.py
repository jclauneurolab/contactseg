"""Three-class tissue segmentation of the T1w with ANTs Atropos.

Only used with --atlas_source antspy, where no external derivatives dataset is
available to supply tissue probabilities. sMRIPrep and freesurfer runs bring
their own, and this rule is not part of those workflows.
"""

import ants


def tissue_seg(t1w, tissue_labels, output_probseg):
    """
    Function that segments a T1w image into tissue probability maps.

    Parameters
    ----------
    t1w : str
        Path to the (bias corrected) subject T1w image.
    tissue_labels : list of str
        Tissue class names, ordered from the darkest to the brightest class in
        a T1w image, i.e. ``["CSF", "GM", "WM"]`` after sorting by intensity.
    output_probseg : list of str
        Paths to save the probability maps, in the order of ``tissue_labels``.

    Returns
    -------
    None
    """

    image = ants.image_read(t1w)
    mask = ants.get_mask(image)

    segmentation = ants.atropos(
        a=image,
        x=mask,
        i="kmeans[3]",
        m="[0.2,1x1x1]",
        c="[5,0]",
    )

    # atropos orders its classes by increasing intensity: CSF, GM, WM
    by_intensity = {"CSF": 0, "GM": 1, "WM": 2}
    for label, probseg_file in zip(tissue_labels, output_probseg):
        ants.image_write(
            segmentation["probabilityimages"][by_intensity[label]],
            probseg_file,
        )


if __name__ == "__main__":
    tissue_seg(
        t1w=snakemake.input.t1w,
        tissue_labels=snakemake.params.tissue_labels,
        output_probseg=snakemake.output.probseg,
    )
