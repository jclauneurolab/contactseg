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
    - df_template: DataFrame containing the points to label.
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
    # fuzzy_list = []
    # voxel_list = []

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
            # Retrieve the region name and hemisphere
            region_info = df_atlas.loc[df_atlas["label"] == voxel_value].iloc[0]
            hemisphere = "Left" if region_info["hemi"] == "L" else "Right"
            region_name = f"{region_info['name']} ({hemisphere})"
        else:
            region_name = np.nan
            voxel_value = np.nan

        # If the region name is nan, look at surrounding labels
        if pd.isna(region_name) and (
            pd.isna(voxel_value)
        ):  # Check if region_name is nan
            print(
                f"Region name is nan for index {i+1}. Searching surrounding labels..."
            )
            non_nan_labels = []
            non_nan_distances = []

            # Iterate over the surrounding voxels
            for x in range(selected_atlasdata.shape[0]):
                for y in range(selected_atlasdata.shape[1]):
                    for z in range(selected_atlasdata.shape[2]):
                        label = selected_atlasdata[x, y, z]
                        if label > 0:  # Only consider non-zero labels
                            distance = distances_mat[x, y, z]
                            non_nan_labels.append(label)
                            non_nan_distances.append(distance)

            if non_nan_labels:
                # Find the closest label based on distance
                closest_index = np.argmin(non_nan_distances)
                voxel_value = non_nan_labels[closest_index]
                region_info = df_atlas.loc[df_atlas["label"] == voxel_value].iloc[0]
                hemisphere = "Left" if region_info["hemi"] == "L" else "Right"
                region_name = f"{region_info['name']} ({hemisphere})"
                print(
                    f"Closest label found for index {i+1}: {region_name} "
                    f"at distance {non_nan_distances[closest_index]:.3f}"
                )
            else:
                print(f"No non-nan labels found for index {i+1}. Label remains nan.")
                region_name = np.nan

        labelnames.append(region_name)
        # fuzzy_list.append(used_fuzzy)
        # voxel_list.append(voxel_value)

    return labelnames


def write_fcsv_with_labels(input_fcsv_path, output_fcsv_path, atlas_labels):
    """
    Copy an FCSV file, replacing the 'label' column (index 11) with atlas labels.

    Parameters:
    - input_fcsv_path: Path to the source FCSV file.
    - output_fcsv_path: Path to write the labelled FCSV file.
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


if __name__ == "__main__":
    # Read the FCSV file, skipping metadata lines starting with '#'
    points_df = pd.read_csv(
        snakemake.input.mni_coords,
        sep=",",  # FCSV files are comma-separated
        comment="#",  # Skip lines starting with '#'
        header=None,  # No header row in the data section
    )

    # Assign column names (adjust based on your FCSV format)
    points_df.columns = [
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

    # Load atlas segmentation and labels
    dseg_nii = nib.load(snakemake.input.atlas_segmentation)
    df_atlas = pd.read_csv(snakemake.input.atlas_labels, sep="\t")

    # Look up atlas labels
    atlas_labels = lookup_atlas_label(
        df_template=points_df,
        coords_columns=["x", "y", "z"],
        dseg_nii=dseg_nii,
        df_atlas=df_atlas,
        fuzzy_dist=snakemake.params.fuzzy_dist,
    )

    write_fcsv_with_labels(
        snakemake.input.t1w_coords,
        snakemake.output.atlas_labelled_t1w_contactseg,
        atlas_labels,
    )
