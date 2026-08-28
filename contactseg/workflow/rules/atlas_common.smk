from functools import partial

from contactseg.workflow.lib import atlas as atlas_lib

# ---- atlas selection -------------------------------------------------------


def get_atlases_by_kind(space, atlas_type):
    """Return the requested atlases matching ``space`` and ``atlas_type``.

    ``space`` is ``"template"`` (defined on the configured volumetric
    template), ``"fsaverage"`` or ``"native"``.
    """
    selected = []
    for atlas in config["atlas"]:
        entry = atlas_lib.get_atlas_entry(config, atlas)
        if space == "template":
            matches = entry["space"] == config["template"]
        else:
            matches = entry["space"] == space
        if matches and entry["type"] == atlas_type:
            selected.append(atlas)

    return selected


def get_surface_atlases():
    """Atlases that are defined on a surface, in any space."""
    return get_atlases_by_kind("fsaverage", "surface") + get_atlases_by_kind(
        "native", "surface"
    )


def needs_template_reg():
    """True if any output requires the subject to template transforms."""
    if get_atlases_by_kind("template", "volume"):
        return True

    # MNI coordinates are reported for every source that can provide them
    return config["atlas_source"] != "freesurfer"


def needs_surfaces():
    """True if surfaces have to be pulled in from freesurfer."""
    return bool(get_surface_atlases()) or config["map_atlas_surfaces"]


# ---- transforms ------------------------------------------------------------


def get_template_t1w():
    """Return the packaged reference image for the configured template."""
    template = config["template"]
    try:
        t1w = config["templates"][template]["T1w"]
    except KeyError as err:
        raise ValueError(
            f"template '{template}' is not defined in the templates config "
            "block"
        ) from err

    return str(atlas_lib.get_atlas_dir(workflow.basedir, config) / t1w)


def get_affine_xfm():
    """Affine part of the subject to template registration (ANTs, LPS)."""
    return bids(
        root=config["output_dir"],
        datatype="atlasreg",
        session="pre",
        desc="affine",
        suffix="xfm",
        extension=".mat",
        **{"from": "T1w", "to": config["template"]},
        **inputs["pre_t1w"].wildcards,
    )


def get_warp_xfm(inverse=False):
    """Deformation field of the subject to template registration."""
    entities = (
        {"from": config["template"], "to": "T1w"}
        if inverse
        else {"from": "T1w", "to": config["template"]}
    )

    return bids(
        root=config["output_dir"],
        datatype="atlasreg",
        session="pre",
        desc="warp",
        suffix="xfm",
        extension=".nii.gz",
        **entities,
        **inputs["pre_t1w"].wildcards,
    )


def get_smriprep_xfm(wildcards, inverse=False):
    """Composite transform written by sMRIPrep for this subject."""
    template = config["template"]
    from_space, to_space = (template, "T1w") if inverse else ("T1w", template)
    fname = (
        f"sub-{wildcards.subject}_from-{from_space}_to-{to_space}"
        "_mode-image_xfm.h5"
    )

    return atlas_lib.find_subject_file(
        atlas_lib.get_derivatives_dir(config, "smriprep_dir"),
        wildcards.subject,
        f"anat/{fname}",
        session=config["derivatives_session"],
    )


def get_transforms(direction):
    """Return the transform chain for ``direction``.

    ``direction`` follows the ANTs image convention: ``T1w_to_template``
    resamples a subject image into template space, ``template_to_T1w``
    resamples a template image (e.g. an atlas) into subject space. Point sets
    travel the opposite way, so contacts are pushed into template space with
    the ``template_to_T1w`` chain.
    """
    inverse = direction == "template_to_T1w"

    if config["atlas_source"] == "smriprep":
        return partial(get_smriprep_xfm, inverse=inverse)

    if inverse:
        return [get_affine_xfm(), get_warp_xfm(inverse=True)]

    return [get_warp_xfm(), get_affine_xfm()]


def get_invert_flags(direction):
    """Return the ``whichtoinvert`` flags matching ``get_transforms``."""
    if config["atlas_source"] == "smriprep":
        return [False]

    # the affine is stored in the subject to template direction, so it has to
    # be inverted when pulling the template into subject space
    return [True, False] if direction == "template_to_T1w" else [False, False]


# ---- coordinates -----------------------------------------------------------


def get_native_coords():
    """Contact coordinates in T1w space, labelled by electrode if available."""
    if config["label"]:
        return rules.label_coords.output.labelled_coords

    return rules.transform_coords.output.transformed_coords


def get_template_coords():
    """Contact coordinates warped into template space."""
    return bids(
        root=config["output_dir"],
        datatype="ieeg",
        session="post",
        space=config["template"],
        desc="contactseg",
        suffix="coords",
        extension=".fcsv",
        **inputs["post_ct"].wildcards,
    )


# ---- atlas files -----------------------------------------------------------


def get_atlas_dseg_in_native():
    """Atlas segmentation carried into subject T1w space by the workflow."""
    return bids(
        root=config["output_dir"],
        datatype="atlasreg",
        session="pre",
        space="T1w",
        atlas="{atlas}",
        suffix="dseg",
        extension=".nii.gz",
        **inputs["pre_t1w"].wildcards,
    )


def get_atlas_lut(wildcards):
    """Lookup table for ``wildcards.atlas``.

    Atlases that declare ``FreeSurferColorLUT`` get an empty input: the lookup
    script falls back to the aseg/aparc table it carries internally, so no
    freesurfer installation is needed to resolve label names.
    """
    entry = atlas_lib.get_atlas_entry(config, wildcards.atlas)
    if entry["lut"] == "FreeSurferColorLUT":
        return []

    return atlas_lib.get_atlas_file(
        workflow.basedir, config, wildcards.atlas, "lut"
    )


def get_key_offsets(atlas):
    """Per-hemisphere offsets applied before merging surface labels."""
    key_offset = atlas_lib.get_atlas_entry(config, atlas).get("key_offset", {})

    return [key_offset.get(hemi, 0) for hemi in config["hemi"]]


def get_fs_atlas_dseg(wildcards):
    """Segmentation shipped inside the freesurfer subject directory."""
    entry = atlas_lib.get_atlas_entry(config, wildcards.atlas)

    return atlas_lib.get_freesurfer_file(
        config,
        wildcards.subject,
        entry["dseg"],
        session=config["derivatives_session"],
    )


def get_template_atlas_dseg(wildcards):
    """Packaged atlas segmentation, in template space."""
    return atlas_lib.get_atlas_file(
        workflow.basedir, config, wildcards.atlas, "dseg"
    )


# ---- freesurfer surfaces ---------------------------------------------------


def get_fs_surface(key):
    """Return an input function for a freesurfer surface file."""

    def _get_surface(wildcards):
        return atlas_lib.get_freesurfer_file(
            config,
            wildcards.subject,
            config["surface_files"][key],
            session=config["derivatives_session"],
            hemi=wildcards.hemi,
        )

    return _get_surface


def get_fs_annot(wildcards):
    """Freesurfer annotation file backing a native surface atlas."""
    entry = atlas_lib.get_atlas_entry(config, wildcards.atlas)

    return atlas_lib.get_freesurfer_file(
        config,
        wildcards.subject,
        entry["label"],
        session=config["derivatives_session"],
        hemi=wildcards.hemi,
    )


def get_surf_gii(surfname, space="T1w"):
    """Subject surface, converted to gifti."""
    return bids(
        root=config["output_dir"],
        datatype="surf",
        session="pre",
        space=space,
        hemi="{hemi}",
        suffix=f"{surfname}.surf.gii",
        **inputs["pre_t1w"].wildcards,
    )


def get_label_gii(space="T1w"):
    """Atlas labels on the subject surface."""
    return bids(
        root=config["output_dir"],
        datatype="surf",
        session="pre",
        space=space,
        hemi="{hemi}",
        atlas="{atlas}",
        suffix="dseg",
        extension=".label.gii",
        **inputs["pre_t1w"].wildcards,
    )


# ---- tissue segmentation ---------------------------------------------------


def get_smriprep_probseg(wildcards):
    """Tissue probability maps written by sMRIPrep for this subject."""
    return [
        atlas_lib.find_subject_file(
            atlas_lib.get_derivatives_dir(config, "smriprep_dir"),
            wildcards.subject,
            f"anat/sub-{wildcards.subject}_label-{label}_probseg.nii.gz",
            session=config["derivatives_session"],
        )
        for label in config["tissue_labels"]
    ]


def get_tissue_probseg():
    """GM/WM/CSF probability maps, from sMRIPrep or computed with ANTs.

    With ``--atlas_source freesurfer`` no probability maps are used: tissue
    classes are derived from the aparc+aseg indices by the lookup script.
    """
    if config["atlas_source"] == "smriprep":
        return get_smriprep_probseg

    return rules.tissue_seg.output.probseg


# ---- targets ---------------------------------------------------------------


def get_atlas_labels_output():
    """Final targets added to rule all when --atlas_labels is set."""
    final = []
    for extension in (".tsv", ".fcsv"):
        final.extend(
            inputs["post_ct"].expand(
                bids(
                    root=config["output_dir"],
                    datatype="ieeg",
                    session="post",
                    space="T1w",
                    desc="atlas",
                    suffix="electrodes",
                    extension=extension,
                    **inputs["post_ct"].wildcards,
                )
            )
        )
    if config["map_atlas_surfaces"]:
        final.extend(
            inputs["pre_t1w"].expand(
                get_label_gii(),
                atlas=config["atlas"],
                hemi=config["hemi"],
            )
        )
    if config["export_nrrd"]:
        final.extend(
            inputs["pre_t1w"].expand(
                bids(
                    root=config["output_dir"],
                    datatype="atlasreg",
                    session="pre",
                    space="T1w",
                    atlas="{atlas}",
                    suffix="dseg.seg.nrrd",
                    **inputs["pre_t1w"].wildcards,
                ),
                atlas=config["atlas"],
            )
        )

    return final
