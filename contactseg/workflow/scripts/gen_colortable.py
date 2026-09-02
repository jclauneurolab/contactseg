"""Write colour tables so a segmentation can be read in 3D Slicer or ITK-SNAP.

A labelmap loaded on its own shows numbers, not anatomy. Both viewers take a
sidecar naming and colouring each index, in formats that differ only in
punctuation, so both are written: Slicer's colour table and ITK-SNAP's label
description file.

Colours come from the atlas when it carries them -- a GIFTI label table or a
lookup table with r/g/b columns -- and are otherwise generated from the label
index, so the same atlas always comes out the same colour.
"""

import colorsys

import pandas as pd

GOLDEN_RATIO = 0.618033988749895

ITKSNAP_HEADER = """\
# ITK-SnAP Label Description File
# File format:
# IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL
################################################
    0     0    0    0     0    0   0    "Clear Label"
"""


def label_color(key, row=None):
    """Return an RGB tuple in 0-255 for ``key``."""
    if row is not None and {"r", "g", "b"} <= set(row.index):
        return int(row["r"]), int(row["g"]), int(row["b"])

    hue = (key * GOLDEN_RATIO) % 1.0
    red, green, blue = colorsys.hsv_to_rgb(hue, 0.6, 0.9)

    return int(red * 255), int(green * 255), int(blue * 255)


def has_hemisphere(name):
    """Whether a structure name already says which hemisphere it is in.

    Some atlases number their parcels per hemisphere and leave the side out of
    the name (freesurfer's DK set: "superiortemporal"); others carry it (the
    Yale atlas: "L_TP1_B"). Prefixing the second kind gives "L_L_TP1_B", so the
    side is added only where the name is silent about it.
    """
    lowered = name.lower()

    return (
        lowered.startswith(
            (
                "l_",
                "r_",
                "lh_",
                "rh_",
                "lh.",
                "rh.",
                "left-",
                "right-",
                "ctx-lh",
                "ctx-rh",
                "wm-lh",
                "wm-rh",
            )
        )
        or "(left)" in lowered
        or "(right)" in lowered
    )


def read_lut(lut_file):
    """Read a dseg.tsv into rows of ``(index, name, (r, g, b))``."""
    table = pd.read_csv(lut_file, sep="\t")
    table.columns = table.columns.str.lower()

    rows = []
    for _, row in table.iterrows():
        index = int(row["label"])
        name = str(row["name"])
        if (
            "hemi" in table.columns
            and str(row["hemi"]) in ("L", "R")
            and not has_hemisphere(name)
        ):
            name = f"{row['hemi']}_{name}"
        rows.append((index, name.replace(" ", "_"), label_color(index, row)))

    return sorted(rows)


def write_slicer_colortable(rows, output_ctbl, atlas):
    """Write a colour table for the Slicer Colors module."""
    with open(output_ctbl, "w") as file:
        file.write(f"# Color table file {atlas}\n")
        file.write(f"# {len(rows) + 1} values\n")
        file.write("0 background 0 0 0 0\n")
        for index, name, (red, green, blue) in rows:
            file.write(f"{index} {name} {red} {green} {blue} 255\n")


def write_itksnap_labels(rows, output_txt):
    """Write an ITK-SNAP label description file."""
    with open(output_txt, "w") as file:
        file.write(ITKSNAP_HEADER)
        for index, name, (red, green, blue) in rows:
            file.write(
                f"{index:>5} {red:>5} {green:>4} {blue:>4} {1.0:>5} "
                f'{1:>4} {1:>3}    "{name}"\n'
            )


def gen_colortable(lut_file, atlas, output_ctbl, output_txt):
    """
    Function that writes the Slicer and ITK-SNAP colour tables for an atlas.

    Parameters
    ----------
    lut_file : str
        Path to the atlas dseg.tsv.
    atlas : str
        Atlas name, written into the Slicer header.
    output_ctbl : str
        Path to save the Slicer colour table.
    output_txt : str
        Path to save the ITK-SNAP label description file.

    Returns
    -------
    None
    """

    rows = read_lut(lut_file)
    write_slicer_colortable(rows, output_ctbl, atlas)
    write_itksnap_labels(rows, output_txt)


def first_or_none(value):
    """Return a single path from an input that may be an empty list."""
    if isinstance(value, (list, tuple)):
        return value[0] if value else None

    return value or None


if __name__ == "__main__":
    gen_colortable(
        lut_file=first_or_none(snakemake.input.lut),
        atlas=snakemake.wildcards.atlas,
        output_ctbl=snakemake.output.colortable,
        output_txt=snakemake.output.itksnap,
    )
