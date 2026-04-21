#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to look up atlas labels for a list of points and save them as a Slicer FCSV file.
"""

import pandas as pd
import numpy as np
import nibabel as nib


def lookup_atlas_label(df_template, coords_columns, dseg_nii, df_atlas, fuzzy_dist=2):
    """
    Look up atlas labels for a list of points.

    Parameters:
    - df_template: DataFrame containing the points in atlas space to label.
    - coords_columns: List of column names for the coordinates (e.g., ['x', 'y', 'z']).
    - dseg_nii: NIfTI image of the atlas segmentation.
    - df_atlas: DataFrame containing atlas labels and their corresponding voxel values.
    - fuzzy_dist: Distance for fuzzy matching when assigning atlas labels.

    Returns:
    - out_df: DataFrame with atlas labels assigned to the points.
    """
    dseg_vol = dseg_nii.get_fdata().astype("int")
    dseg_affine = dseg_nii.affine
    coords = df_template[coords_columns].to_numpy()
    labelnames = []

    for i in range(len(coords)):
        coords_vx = np.round(
            nib.affines.apply_affine(np.linalg.inv(dseg_affine), coords[i, :])
        ).astype(int)
        used_fuzzy = 0
        try:
            voxel_value = dseg_vol[coords_vx[0], coords_vx[1], coords_vx[2]]
        except IndexError:
            voxel_value = 0

        if (fuzzy_dist is not None) and (voxel_value == 0):
            fuzzy_diameter = fuzzy_dist * 2 + 1
            distances_mat = np.zeros((fuzzy_diameter, fuzzy_diameter, fuzzy_diameter))
            for x in range(fuzzy_diameter):
                for y in range(fuzzy_diameter):
                    for z in range(fuzzy_diameter):
                        distances_mat[x, y, z] = np.linalg.norm(
                            np.array([fuzzy_dist, fuzzy_dist, fuzzy_dist])
                            - np.array([x, y, z])
                        )
            # check if the distances box will exceed image boundaries (super edgy case)
            trim = np.zeros((3, 2), dtype=int)
            for j in range(3):
                if coords_vx[j] - fuzzy_dist < 0:
                    trim[j, 0] = fuzzy_dist - coords_vx[j]
                if coords_vx[j] + fuzzy_dist + 1 > dseg_vol.shape[j]:
                    trim[j, 1] = coords_vx[j] + fuzzy_dist + 1 - dseg_vol.shape[j]

            assert np.all(trim >= 0), "Trim (" + str(trim) + ") should be non-negative"
            # get nearest voxel that is not zero, but less than specified voxels away
            selected_atlasdata = dseg_vol[
                (coords_vx[0] - fuzzy_dist + trim[0, 0]) : (
                    coords_vx[0] + fuzzy_dist + 1 - trim[0, 1]
                ),
                (coords_vx[1] - fuzzy_dist + trim[1, 0]) : (
                    coords_vx[1] + fuzzy_dist + 1 - trim[1, 1]
                ),
                (coords_vx[2] - fuzzy_dist + trim[2, 0]) : (
                    coords_vx[2] + fuzzy_dist + 1 - trim[2, 1]
                ),
            ]
            trimmed_distances = distances_mat[
                trim[0, 0] : (-1 * trim[0, 1]) if trim[0, 1] != 0 else None,
                trim[1, 0] : (-1 * trim[1, 1]) if trim[1, 1] != 0 else None,
                trim[2, 0] : (-1 * trim[2, 1]) if trim[2, 1] != 0 else None,
            ]

            distances = np.ma.masked_where(
                (selected_atlasdata == 0) | (trimmed_distances > fuzzy_dist),
                trimmed_distances,
            )
            nearest_voxel = np.unravel_index(np.argmin(distances), distances.shape)
            voxel_value = selected_atlasdata[nearest_voxel]
            nearest_distance = distances[nearest_voxel]
            used_fuzzy = 1

        if voxel_value > 0:
             # Dynamically find the correct column name for the ID
            possible_id_cols = ["index", "label", "id", "value"]
            id_col = next((col for col in possible_id_cols if col in df_atlas.columns), None)
            
            if not id_col:
                raise ValueError(f"Could not find an ID column in atlas TSV. Found: {list(df_atlas.columns)}")
            
            # Retrieve the region name and hemisphere
            region_info = df_atlas.loc[df_atlas[id_col] == voxel_value].iloc[0]
            if "hemi" in region_info.index and pd.notna(region_info["hemi"]):
                hemisphere = "Left" if region_info["hemi"] == "L" else "Right"
                region_name = f"{region_info['name']} ({hemisphere})"
            else:
                region_name = str(region_info['name'])
        else:
            region_name = np.nan
            voxel_value = np.nan

        labelnames.append(region_name)

    return labelnames


def white_vs_grey_label(coords_list, dseg_path):

    # Load the native space segmentation
    img = nib.load(dseg_path)
    data = img.get_fdata()

    # Get the inverse affine to map physical (x,y,z) to voxel (i,j,k)
    inv_affine = np.linalg.inv(img.affine)

    tissue_labels = []

    for x, y, z in coords_list:

        coord_ras = [x, y, z]

        # Transform physical coordinate to voxel index
        voxel_idx = nib.affines.apply_affine(inv_affine, coord_ras)
        i, j, k = np.round(voxel_idx).astype(int)

        # Check if the voxel is inside the image bounds
        if (
            (0 <= i < data.shape[0])
            and (0 <= j < data.shape[1])
            and (0 <= k < data.shape[2])
        ):
            label_val = data[i, j, k]

            # Standard ANTs/sMRIPrep dseg values: 1=CSF, 2=Grey Matter, 3=White Matter
            if label_val == 1:
                tissue_labels.append("Grey Matter")
            elif label_val == 2:
                tissue_labels.append("White Matter")
            elif label_val == 3:
                tissue_labels.append("CSF")
            else:
                tissue_labels.append("Other")
        else:
            tissue_labels.append("Out of Bounds")

    return tissue_labels

def tissue_probability(coords_list, prob_seg_GM, prob_seg_WM, prob_seg_CSF):

    prob_seg_GM = nib.load(prob_seg_GM)
    inv_affine = np.linalg.inv(prob_seg_GM.affine)
    prob_seg_GM = prob_seg_GM.get_fdata()
    prob_seg_WM = nib.load(prob_seg_WM).get_fdata()
    prob_seg_CSF = nib.load(prob_seg_CSF).get_fdata()

    tissue_probs = []
    for x, y, z in coords_list:

        coord_ras = [x, y, z]

        # Transform physical coordinate to voxel index
        voxel_idx = nib.affines.apply_affine(inv_affine, coord_ras)
        i, j, k = np.round(voxel_idx).astype(int)

        if (
            (0 <= i < prob_seg_GM.shape[0])
            and (0 <= j < prob_seg_GM.shape[1])
            and (0 <= k < prob_seg_GM.shape[2])
        ):
            prob_GM = prob_seg_GM[i, j, k]
            prob_WM = prob_seg_WM[i, j, k]
            prob_CSF = prob_seg_CSF[i, j, k]
            tissue_probs.append((prob_GM, prob_WM, prob_CSF))
    return tissue_probs

def write_fcsv_with_labels(input_fcsv_path, output_fcsv_path, atlas_labels):
    """
    Copy an FCSV file, replacing the 'label' column (index 11) with atlas labels.

    Parameters:
    - input_fcsv_path: Path to the source FCSV file.
    - output_fcsv_path: Path to write the atlas labelled FCSV file.
    - atlas_labels: List of label strings aligned to the data rows.
    """
    with open(input_fcsv_path, "r") as f:
        lines = f.readlines()

    header = [l for l in lines if l.startswith("#")]
    data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

    with open(output_fcsv_path, "w") as f:
        f.writelines(header)
        for i, line in enumerate(data_lines):
            parts = line.strip().split(",")
            while len(parts) < 12:
                parts.append("")
            parts[11] = str(atlas_labels[i])
            f.write(",".join(parts) + "\n")

FCSV_COLUMNS = [
    "id",
    "x",
    "y",
    "z",
    "ow",
    "ox",
    "oy",
    "oz",
    "vis",
    "sel",
    "lock",
    "label",
    "desc",
    "associatedNodeID",
]

def load_fcsv(path):
    return pd.read_csv(path, sep=",", comment="#", header=None, names=FCSV_COLUMNS)

if __name__ == "__main__":

    native_space = snakemake.params.native_space
    GWmatter_labels = snakemake.params.GWmatter_labels

    native_df = load_fcsv(snakemake.input.native_coords)
    mni_df = load_fcsv(snakemake.input.mni_coords)
    df_atlas = pd.read_csv(snakemake.input.atlas_labels, sep="\t")

    tissue_type = "Unknown"
     # Initialize defaults in case GWmatter_labels is False
    num_coords = len(native_df)
    tissue_type = ["Unknown"] * num_coords
    tissue_probs = [(np.nan, np.nan, np.nan)] * num_coords


    if snakemake.params.GWmatter_labels:
        tissue_type = white_vs_grey_label(
            native_df[["x", "y", "z"]].values, snakemake.input.native_dseg
        )
        tissue_probs = tissue_probability(
            native_df[["x", "y", "z"]].values,
            snakemake.input.native_prob_seg_GM,
            snakemake.input.native_prob_seg_WM,
            snakemake.input.native_prob_seg_CSF
        )

    if native_space:
        print("Using native space for atlas labeling.")
        dseg_nii = nib.load(snakemake.input.atlas_segmentation_in_native)
        active_df = native_df
    else:
        print("Using MNI space for atlas labeling.")
        dseg_nii = nib.load(snakemake.input.atlas_segmentation_in_mni)
        active_df = mni_df

    atlas_labels = lookup_atlas_label(
        df_template=active_df,
        coords_columns=["x", "y", "z"],
        dseg_nii=dseg_nii,
        df_atlas=df_atlas,
        fuzzy_dist=snakemake.params.fuzzy_dist,
    )

    write_fcsv_with_labels(
        snakemake.input.native_coords,
        snakemake.output.atlas_labelled_t1w_contactseg,
        atlas_labels,
    )

    # Output to a CSV file containing both Native and MNI coordinates
    csv_df = pd.DataFrame(
        {
            "Original_Label": native_df["label"],
            "Atlas_Label": atlas_labels,
            "Tissue_Type": tissue_type,
            "GM_Probability": [p[0] for p in tissue_probs],
            "WM_Probability": [p[1] for p in tissue_probs],
            "CSF_Probability": [p[2] for p in tissue_probs],
            "Native_X": native_df["x"],
            "Native_Y": native_df["y"],
            "Native_Z": native_df["z"],
            "MNI_X": mni_df["x"],
            "MNI_Y": mni_df["y"],
            "MNI_Z": mni_df["z"],
        }
    )
    csv_df.to_csv(snakemake.output.csv_file, index=False)
