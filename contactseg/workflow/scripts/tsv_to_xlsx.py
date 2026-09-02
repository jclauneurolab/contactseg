"""Write a tab-separated table as an Excel workbook.

The tsv stays the machine-readable output; this is the copy that opens in
Excel without an import dialog mangling contact names into dates. Columns are
widened to their contents and the header row is frozen, since these tables are
read by scrolling down a long list of contacts.

It lives in its own rule with its own environment rather than being written
alongside the tsv: openpyxl is needed for nothing else, and adding it to the
shared analysis environment would change that environment's hash and re-run
every rule that uses it.
"""

import pandas as pd

MAX_COLUMN_WIDTH = 60


def column_width(series, header):
    """Return a column width that fits the header and its values."""
    widest = max((len(str(value)) for value in series), default=0)

    return min(max(widest, len(str(header))) + 2, MAX_COLUMN_WIDTH)


def tsv_to_xlsx(input_tsv, output_xlsx, sheet_name="labels"):
    """
    Function that converts a tsv table to an Excel workbook.

    Parameters
    ----------
    input_tsv : str
        Path to the tab-separated table.
    output_xlsx : str
        Path to save the workbook.
    sheet_name : str
        Name of the sheet to write. Excel allows 31 characters, so longer
        names are truncated rather than failing the write.

    Returns
    -------
    pandas.DataFrame
    """

    table = pd.read_csv(input_tsv, sep="\t")
    sheet_name = str(sheet_name)[:31]

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        table.to_excel(writer, sheet_name=sheet_name, index=False)
        sheet = writer.sheets[sheet_name]
        sheet.freeze_panes = "A2"
        for i, name in enumerate(table.columns, start=1):
            letter = sheet.cell(row=1, column=i).column_letter
            sheet.column_dimensions[letter].width = column_width(table[name], name)

    return table


if __name__ == "__main__":
    tsv_to_xlsx(
        input_tsv=snakemake.input.tsv,
        output_xlsx=snakemake.output.xlsx,
        sheet_name=snakemake.wildcards.get("atlas", "labels"),
    )
