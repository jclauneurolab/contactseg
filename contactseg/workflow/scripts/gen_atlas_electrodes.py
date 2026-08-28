"""Collect the per-atlas lookups into one electrodes table.

Writes a BIDS-style electrodes.tsv carrying, for every contact, its native and
template coordinates, one label column per atlas, and the tissue class. The
same table is written as an fcsv so the labelled contacts can be dropped
straight into 3D Slicer.
"""

import pandas as pd

FCSV_HEADER_ROWS = 3
FCSV_LABEL_COLUMN = 11


def read_coords(coords_fcsv):
    """Read contact names and coordinates from a Slicer fcsv file."""
    contacts = pd.read_csv(coords_fcsv, skiprows=FCSV_HEADER_ROWS, header=None)

    return pd.DataFrame(
        {
            "name": contacts[FCSV_LABEL_COLUMN].astype(str),
            "x": contacts[1],
            "y": contacts[2],
            "z": contacts[3],
        }
    )


def write_fcsv(coords_fcsv, output_fcsv, labels):
    """Copy an fcsv, replacing the label column with ``labels``."""
    with open(coords_fcsv, "r") as file:
        lines = file.readlines()

    header = [line for line in lines if line.startswith("#")]
    rows = [line for line in lines if not line.startswith("#") and line.strip()]

    with open(output_fcsv, "w") as file:
        file.writelines(header)
        for i, line in enumerate(rows):
            fields = line.strip().split(",")
            while len(fields) <= FCSV_LABEL_COLUMN:
                fields.append("")
            fields[FCSV_LABEL_COLUMN] = str(labels[i])
            file.write(",".join(fields) + "\n")


def gen_atlas_electrodes(
    coords_fcsv,
    template_coords_fcsv,
    atlas_label_files,
    tissue_file,
    atlases,
    template,
    output_tsv,
    output_fcsv,
):
    """
    Function that merges the atlas lookups into one electrodes table.

    Parameters
    ----------
    coords_fcsv : str
        Path to the contact coordinates in subject T1w space.
    template_coords_fcsv : str or None
        Path to the same contacts in template space, if computed.
    atlas_label_files : list of str
        Per-atlas label tables, in the order of ``atlases``.
    tissue_file : str or None
        Tissue probability table, when tissue came from a probabilistic
        segmentation rather than from the atlas indices.
    atlases : list of str
        Atlas names, used to prefix the per-atlas columns.
    template : str
        Template name, used to name the template coordinate columns.
    output_tsv : str
        Path to save the merged table.
    output_fcsv : str
        Path to save the same table as a Slicer fcsv.

    Returns
    -------
    pandas.DataFrame
    """

    electrodes = read_coords(coords_fcsv)

    if template_coords_fcsv:
        template_coords = read_coords(template_coords_fcsv)
        for axis in ("x", "y", "z"):
            electrodes[f"{template}_{axis}"] = template_coords[axis]

    primary_label = None
    for atlas, label_file in zip(atlases, atlas_label_files):
        labels = pd.read_csv(label_file, sep="\t")
        electrodes = electrodes.merge(
            labels.add_prefix(f"{atlas}_").rename(columns={f"{atlas}_name": "name"}),
            on="name",
            how="left",
        )
        if primary_label is None:
            primary_label = electrodes[f"{atlas}_structure"]

    if tissue_file:
        tissue = pd.read_csv(tissue_file, sep="\t")
        electrodes = electrodes.merge(tissue, on="name", how="left")

    electrodes.to_csv(
        output_tsv, sep="\t", index=False, float_format="%.3f", na_rep="n/a"
    )
    write_fcsv(coords_fcsv, output_fcsv, primary_label.fillna("None").tolist())

    return electrodes


def first_or_none(value):
    """Return a single path from an input that may be an empty list."""
    if isinstance(value, (list, tuple)):
        return value[0] if value else None

    return value or None


if __name__ == "__main__":
    gen_atlas_electrodes(
        coords_fcsv=snakemake.input.coords,
        template_coords_fcsv=first_or_none(snakemake.input.template_coords),
        atlas_label_files=snakemake.input.atlas_labels,
        tissue_file=first_or_none(snakemake.input.tissue_labels),
        atlases=snakemake.params.atlases,
        template=snakemake.params.template,
        output_tsv=snakemake.output.electrodes_tsv,
        output_fcsv=snakemake.output.electrodes_fcsv,
    )
