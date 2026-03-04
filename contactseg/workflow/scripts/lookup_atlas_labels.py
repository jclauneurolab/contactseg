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
    dseg_vol = dseg_nii.get_fdata().astype('int')
    dseg_affine = dseg_nii.affine
    coords = df_template[coords_columns].to_numpy()
    labelnames = []
    fuzzy_list = []
    voxel_list = []

    for i in range(len(coords)):
        print(f"Processing electrode {i}: {df_template.iloc[i]['label']} at coords {coords[i]}")
        coords_vx = np.round(nib.affines.apply_affine(np.linalg.inv(dseg_affine), coords[i, :])).astype(int)
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
                            np.array([fuzzy_dist, fuzzy_dist, fuzzy_dist]) - np.array([x, y, z])
                        )

            trim = np.zeros((3, 2), dtype=int)
            for j in range(3):
                if coords_vx[j] - fuzzy_dist < 0:
                    trim[j, 0] = fuzzy_dist - coords_vx[j]
                if coords_vx[j] + fuzzy_dist + 1 > dseg_vol.shape[j]:
                    trim[j, 1] = coords_vx[j] + fuzzy_dist + 1 - dseg_vol.shape[j]

            assert np.all(trim >= 0), 'Trim (' + str(trim) + ') should be non-negative'

            selected_atlasdata = dseg_vol[
                (coords_vx[0] - fuzzy_dist + trim[0, 0]):(coords_vx[0] + fuzzy_dist + 1 - trim[0, 1]),
                (coords_vx[1] - fuzzy_dist + trim[1, 0]):(coords_vx[1] + fuzzy_dist + 1 - trim[1, 1]),
                (coords_vx[2] - fuzzy_dist + trim[2, 0]):(coords_vx[2] + fuzzy_dist + 1 - trim[2, 1])
            ]
            trimmed_distances = distances_mat[
                trim[0, 0]:(-1 * trim[0, 1]) if trim[0, 1] != 0 else None,
                trim[1, 0]:(-1 * trim[1, 1]) if trim[1, 1] != 0 else None,
                trim[2, 0]:(-1 * trim[2, 1]) if trim[2, 1] != 0 else None
            ]

            distances = np.ma.masked_where((selected_atlasdata == 0) | (trimmed_distances > fuzzy_dist), trimmed_distances)
            nearest_voxel = np.unravel_index(np.argmin(distances), distances.shape)
            voxel_value = selected_atlasdata[nearest_voxel]
            used_fuzzy = 1

        if voxel_value > 0:
            region_name = df_atlas.loc[df_atlas['label'] == voxel_value, 'name'].to_list()[0]
        else:
            region_name = np.nan
            voxel_value = np.nan

        labelnames.append(region_name)
        fuzzy_list.append(used_fuzzy)
        voxel_list.append(voxel_value)

    out_df = df_template.copy()
    out_df['atlas_label'] = labelnames
    out_df['fuzzy'] = fuzzy_list
    out_df['vox_val'] = voxel_list
    return out_df


def save_as_fcsv(df, output_path):
    """
    Save the labeled points as a Slicer FCSV file.

    Parameters:
    - df: DataFrame containing the labeled points.
    - output_path: Path to save the FCSV file.
    """
    with open(output_path, 'w') as f:
        # Write FCSV metadata header
        f.write("# Markups fiducial file version = 4.11\n")
        f.write("# CoordinateSystem = 0\n")
        f.write("# Columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n")

        # Write the data rows
        for i, row in df.iterrows():
            f.write(f"F-{i+1},{row['x']},{row['y']},{row['z']},0,0,0,1,1,1,0,{row['atlas_label']},,\n")


if __name__ == "__main__":
    # Read the FCSV file, skipping metadata lines starting with '#'
    points_df = pd.read_csv(
        snakemake.input.points,
        sep=',',  # FCSV files are comma-separated
        comment='#',  # Skip lines starting with '#'
        header=None  # No header row in the data section
    )

    # Assign column names (adjust based on your FCSV format)
    points_df.columns = ['id', 'x', 'y', 'z', 'ow', 'ox', 'oy', 'oz', 'vis', 'sel', 'lock', 'label', 'desc', 'associatedNodeID']

    # Load atlas segmentation and labels
    dseg_nii = nib.load(snakemake.input.atlas_segmentation)
    df_atlas = pd.read_csv(snakemake.input.atlas_labels, sep='\t')

    # Look up atlas labels
    labeled_points = lookup_atlas_label(
        df_template=points_df,
        coords_columns=['x', 'y', 'z'],
        dseg_nii=dseg_nii,
        df_atlas=df_atlas,
        fuzzy_dist=snakemake.params.fuzzy_dist
    )

    # Save the labeled points as an FCSV file
    save_as_fcsv(labeled_points, snakemake.output.labeled_points)
    print(f"Labeled points saved to {snakemake.output.labeled_points}")