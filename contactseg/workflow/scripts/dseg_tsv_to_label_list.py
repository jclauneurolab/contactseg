"""Write a workbench label list from a BIDS dseg.tsv.

``wb_command -metric-label-import`` expects pairs of lines: a name, then the
key with an RGBA colour. Colours are taken from the lookup table when it has
r/g/b columns and generated deterministically from the key otherwise, so the
same atlas always comes out the same colour.
"""

import colorsys

import pandas as pd

GOLDEN_RATIO = 0.618033988749895


def label_color(key, table_row=None):
    """Return an RGBA tuple in 0-255 for ``key``."""
    if table_row is not None and {"r", "g", "b"} <= set(table_row.index):
        return (
            int(table_row["r"]),
            int(table_row["g"]),
            int(table_row["b"]),
            255,
        )

    hue = (key * GOLDEN_RATIO) % 1.0
    red, green, blue = colorsys.hsv_to_rgb(hue, 0.6, 0.9)

    return int(red * 255), int(green * 255), int(blue * 255), 255


def dseg_tsv_to_label_list(lut_file, output_txt):
    """
    Function that converts a dseg.tsv into a workbench label list file.

    Parameters
    ----------
    lut_file : str
        Path to a BIDS dseg.tsv with at least ``label`` and ``name`` columns.
    output_txt : str
        Path to save the label list.

    Returns
    -------
    None
    """

    table = pd.read_csv(lut_file, sep="\t")
    table.columns = table.columns.str.lower()

    with open(output_txt, "w") as file:
        for _, row in table.iterrows():
            key = int(row["label"])
            name = str(row["name"])
            if "hemi" in table.columns and str(row["hemi"]) in ("L", "R"):
                name = f"{name}_{row['hemi']}"
            red, green, blue, alpha = label_color(key, row)
            file.write(f"{name}\n{key} {red} {green} {blue} {alpha}\n")


if __name__ == "__main__":
    dseg_tsv_to_label_list(
        lut_file=snakemake.input.lut,
        output_txt=snakemake.output.label_list,
    )
