import pandas as pd
import ants

def warp_contacts_to_mni(input_fcsv, output_fcsv, affine, affine_syn, forward_warp):
    contacts_df = pd.read_csv(input_fcsv, skiprows=3, header=None)
    points = contacts_df[[1, 2, 3]].copy()
    points.columns = ['x', 'y', 'z']

    # Convert RAS to LPS
    points['x'] *= -1
    points['y'] *= -1

    # Apply Transform
    transformed_points = ants.apply_transforms_to_points(
        dim=3,
        points=points,
        transformlist=[forward_warp, affine_syn, affine],
        whichtoinvert=[False, True, True] 
    )

    # Convert LPS back to RAS
    transformed_points['x'] *= -1
    transformed_points['y'] *= -1

    # Save
    contacts_df[1] = transformed_points['x'].values
    contacts_df[2] = transformed_points['y'].values
    contacts_df[3] = transformed_points['z'].values
    
    # Write output with original Slicer header
    with open(input_fcsv, 'r') as f:
        header = [line for line in f if line.startswith("#")]
        
    with open(output_fcsv, 'w') as f:
        f.writelines(header)
        # mode='a' appends the data after the header
        contacts_df.to_csv(f, index=False, header=False, mode='a', float_format="%.3f")

if __name__ == "__main__":
    warp_contacts_to_mni(
        input_fcsv=snakemake.input.coords,
        output_fcsv=snakemake.output.output_coords,
        affine=snakemake.input.affine,
        forward_warp=snakemake.input.forward_warp,
        affine_syn=snakemake.input.affine_syn
    )