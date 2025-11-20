import json

import pandas as pd


def elec_manufacturer(elec_type):
    if elec_type.startswith("RD1"):
        return "AdTech"
    elif elec_type.startswith("D08"):
        return "Dixi"
    else:
        return "other"


def elec_size(elec_type):
    if elec_type.startswith("RD1"):
        return 6.76
    elif elec_type.startswith("D08"):
        return 6
    else:
        return None


def elec_material(elec_type):
    if elec_type.startswith("RD1"):
        return "platinum"
    elif elec_type.startswith("D08"):
        return "platinum-iridium"
    else:
        return "NaN"


def gen_ieeg_tsv(input_fcsv, output_tsv):
    data = pd.read_csv(input_fcsv, skiprows=3, header=None)
    data[12] = data[12].astype(str)
    data["manufacturer"] = data[12].apply(elec_manufacturer)
    data["size"] = data[12].apply(elec_size)
    data["material"] = data[12].apply(elec_material)

    tsv_data = {
        "name": data[11],
        "x": data[1],
        "y": data[2],
        "z": data[3],
        "size": data["size"],
        "material": data["material"],
        "manufacturer": data["manufacturer"],
        "group": ["depth"] * len(data),
        "hemisphere": data[11].str[0],
    }

    tsv_df = pd.DataFrame(tsv_data)
    tsv_df.to_csv(output_tsv, sep="\t", index=False)


def gen_coordsystem_json(output_json_path, ref_img_path):
    json_data = {
        "iEEGCoordinateSystem": "ScanRAS",
        "iEEGCoordinateUnits": "mm",
        "iEEGCoordinateProcessingDescription": "SEEG contacts localized with contactseg",
        "iEEGCoordinateProcessingReference": "https://github.com/jclauneurolab/contactseg",
        "IntendedFor": f"{ref_img_path}",
    }
    with open(output_json_path, "w") as f:
        json.dump(json_data, f)


if __name__ == "__main__":
    ref_img = snakemake.input["ref_ct"]
    fcsv = snakemake.input["fcsv"]

    tsv_file = snakemake.output["electrodes_tsv"]
    json_file = snakemake.output["coordsystem_json"]

    gen_ieeg_tsv(input_fcsv=fcsv, output_tsv=tsv_file)
    gen_coordsystem_json(json_file, ref_img)
