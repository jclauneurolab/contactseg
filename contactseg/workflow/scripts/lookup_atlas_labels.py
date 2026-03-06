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
			# check if the distances box will exceed image boundaries (super edgy case)
            trim = np.zeros((3, 2), dtype=int)
            for j in range(3):
                if coords_vx[j] - fuzzy_dist < 0:
                    trim[j, 0] = fuzzy_dist - coords_vx[j]
                if coords_vx[j] + fuzzy_dist + 1 > dseg_vol.shape[j]:
                    trim[j, 1] = coords_vx[j] + fuzzy_dist + 1 - dseg_vol.shape[j]

            assert np.all(trim >= 0), 'Trim (' + str(trim) + ') should be non-negative'
			# get nearest voxel that is not zero, but less than specified voxels away
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
            # Retrieve the region name and hemisphere
            region_info = df_atlas.loc[df_atlas['label'] == voxel_value].iloc[0]
            hemisphere = "Left" if region_info['hemi'] == 'L' else "Right"
            region_name = f"{region_info['name']} ({hemisphere})"
        else:
            region_name = np.nan
            voxel_value = np.nan

        labelnames.append(region_name)
        fuzzy_list.append(used_fuzzy)
        voxel_list.append(voxel_value)

    out_df = df_template.copy()

    # Preserve the original electrode label under a new column called 'electrode'
    out_df['electrode'] = df_template['label']

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
        f.write("# Columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,electrode,desc,associatedNodeID\n")

        # Write the data rows
        for i, row in df.iterrows():
            f.write(f"{i+1},{row['x']},{row['y']},{row['z']},0,0,0,1,1,1,0,{row['atlas_label']},{row['electrode']},,\n")

def add_labels_to_fcsv(input_fcsv_path, output_fcsv_path, atlas_labels):
    with open(input_fcsv_path, 'r') as f_in:
        lines = f_in.readlines()

    header = [l for l in lines if l.startswith('#')]
    data_lines = [l for l in lines if not l.startswith('#') and l.strip()]

    with open(output_fcsv_path, 'w') as f_out:
        # Write header exactly as it was
        for line in header:
            f_out.write(line)

        for i, line in enumerate(data_lines):
            parts = line.strip().split(',')
            
            # parts[12] is the 'desc' (Description) column (index starts at 0)
            # We replace it with our atlas label
            atlas_val = str(atlas_labels[i]).replace('nan', 'Unknown')
            
            if len(parts) >= 13:
                parts[12] = atlas_val
            else:
                # If for some reason the row is short, we extend it
                while len(parts) < 13:
                    parts.append("")
                parts[12] = atlas_val
            
            # Rejoin the line with commas
            f_out.write(",".join(parts) + "\n")

if __name__ == "__main__":
    # Read the FCSV file, skipping metadata lines starting with '#'
    points_df = pd.read_csv(
        snakemake.input.mni_coords,
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

    # Extract the produced labels dynamically
    produced_labels = labeled_points['atlas_label'].to_numpy()
    input_fcsv_path = snakemake.input.t1w_coords
    output_fcsv_path = snakemake.output.atlas_labelled_t1w_contactseg

    # Save the labeled points as an FCSV file in MNI space
    save_as_fcsv(labeled_points, snakemake.output.atlas_labelled_mni_contactseg)

    add_labels_to_fcsv(
        input_fcsv_path=input_fcsv_path,
        output_fcsv_path=output_fcsv_path,
        atlas_labels=produced_labels
    )