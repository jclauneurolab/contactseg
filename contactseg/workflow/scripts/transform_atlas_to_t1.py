import ants

def transform_atlas_to_t1(atlas_image, t1_image, inverse_warp, affine, output_atlas_t1):
    # Read the atlas image and T1 image
    atlas = ants.image_read(atlas_image)
    t1 = ants.image_read(t1_image)

    # Apply the inverse transformations (MNI → T1)
    atlas_in_t1 = ants.apply_transforms(
        fixed=t1,  
        moving=atlas, 
        transformlist=[inverse_warp, affine],  
        whichtoinvert=[False, True],  
        interpolator="linear"
    )

    # Save the transformed atlas image
    ants.image_write(atlas_in_t1, output_atlas_t1)


if __name__ == "__main__":
    transform_atlas_to_t1(
        atlas_image=snakemake.input.atlas_image,
        t1_image=snakemake.input.t1_image,
        inverse_warp=snakemake.input.inverse_warp,
        affine=snakemake.input.affine,
        output_atlas_t1=snakemake.output.atlas_in_t1
    )