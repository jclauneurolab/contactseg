import ants


def transform_atlas_to_t1(atlas_image, t1_image, transforms, output_atlas_t1):
    # Read the atlas image and T1 image
    atlas = ants.image_read(atlas_image)
    t1 = ants.image_read(t1_image)

    transforms = snakemake.input.transforms
    if isinstance(transforms, str):
        transforms = [transforms]

    if len(transforms) == 3:
        inv_flags = [True, True, False]
    else:
        inv_flags = [False]

    # Apply the inverse transformations (MNI → T1)
    atlas_in_t1 = ants.apply_transforms(
        fixed=t1,
        moving=atlas,
        transformlist=transforms,
        whichtoinvert=inv_flags,
        interpolator="nearestNeighbor",
    )

    # Save the transformed atlas image
    ants.image_write(atlas_in_t1, output_atlas_t1)


if __name__ == "__main__":
    transform_atlas_to_t1(
        atlas_image=snakemake.input.atlas_image,
        t1_image=snakemake.input.t1_image,
        transforms=snakemake.input.transforms,
        output_atlas_t1=snakemake.output.atlas_in_t1,
    )
