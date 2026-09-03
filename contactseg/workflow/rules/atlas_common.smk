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


# ---- derivatives bridge ----------------------------------------------------


def bridged_sources():
    """Derivatives datasets whose frame has to be aligned to the T1w.

    sMRIPrep and FreeSurfer are frequently run on a different acquisition than
    the one contactseg localises against -- a non-contrast T1w for the same
    patient, say. Everything they produce then sits in that image's scanner RAS
    frame, offset from the pre_t1w by however much the patient moved between
    the two scans. Each dataset in use gets its own bridge, because a run that
    takes transforms from sMRIPrep and surfaces from FreeSurfer is reading two
    frames, not one.
    """
    if config["derivatives_reg"] != "rigid":
        return []

    sources = []
    if config["atlas_source"] in ("smriprep", "freesurfer"):
        sources.append(config["atlas_source"])
    if needs_surfaces() and "freesurfer" not in sources:
        sources.append("freesurfer")

    return sources


def find_derivatives_anat(source, subject):
    """Path to the anatomical ``source`` was computed on, for ``subject``.

    The image can live inside the derivatives dataset (its own reference
    image, which is by definition the frame its transforms are expressed in)
    or in the raw BIDS dataset (the acquisition that was fed to it). Which one
    is used is set per source in the ``derivatives_anat`` config block.
    """
    session = config["derivatives_session"]
    root, relpath = atlas_lib.get_derivatives_anat_spec(config, source)
    relpath = relpath.format(subject=subject, session=session)

    if root == "bids":
        return atlas_lib.find_subject_file(
            config["bids_dir"], subject, relpath, session=session
        )

    if source == "freesurfer":
        return atlas_lib.get_freesurfer_file(
            config, subject, relpath, session=session
        )

    return atlas_lib.find_subject_file(
        atlas_lib.get_derivatives_dir(config, "smriprep_dir"),
        subject,
        relpath,
        session=session,
    )


def get_derivatives_anat(wildcards):
    """Input function keyed on the ``deriv`` wildcard."""
    return find_derivatives_anat(wildcards.deriv, wildcards.subject)


def derivatives_anat_for(source):
    """Return an input function for a fixed derivatives dataset."""

    def _get_anat(wildcards):
        return find_derivatives_anat(source, wildcards.subject)

    return _get_anat


def get_bridge_xfm(source=None, extension=".txt"):
    """4x4 RAS matrix taking a derivatives anatomical to the T1w.

    Returns an empty list when ``source`` needs no bridge, so the rules that
    consume it can declare it as an optional input. ``source`` defaults to the
    dataset the transforms come from; pass ``"freesurfer"`` for the rules that
    read surfaces or segmentations from the subject directory.
    """
    source = source or config["atlas_source"]
    if source not in bridged_sources():
        return []

    return bids(
        root=config["output_dir"],
        datatype="atlasreg",
        session="pre",
        mode="image",
        suffix="xfm",
        extension=extension,
        **{"from": source, "to": "T1w"},
        **inputs["pre_t1w"].wildcards,
    )


def get_bridge_xfm_pattern(extension=".txt"):
    """The same path with the derivatives dataset left as a wildcard."""
    return bids(
        root=config["output_dir"],
        datatype="atlasreg",
        session="pre",
        mode="image",
        suffix="xfm",
        extension=extension,
        **{"from": "{deriv}", "to": "T1w"},
        **inputs["pre_t1w"].wildcards,
    )


def get_atlas_reference():
    """Reference grid the template atlas is resampled onto.

    With a bridge in play the ANTs chain only reaches the derivatives
    anatomical, so that image is the reference and the bridge is folded into
    the output affine afterwards.
    """
    if config["atlas_source"] in bridged_sources():
        return derivatives_anat_for(config["atlas_source"])

    return rules.n4biascorr.output.corrected_t1w


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


def get_surface_atlas_lut():
    """Lookup table written alongside a surface atlas mapped into the volume."""
    return bids(
        root=config["output_dir"],
        datatype="atlasreg",
        session="pre",
        space="T1w",
        atlas="{atlas}",
        suffix="dseg",
        extension=".tsv",
        **inputs["pre_t1w"].wildcards,
    )


def get_atlas_labels(extension=".tsv"):
    """Per-atlas contact labels.

    These live in ``ses-pre/atlasreg`` beside the segmentation they were read
    out of, rather than with the contacts: the atlas, its lookup table and its
    colour tables are all there, and the labels are that segmentation sampled
    at the contact positions. The merged electrodes table stays in
    ``ses-post/ieeg``, where the contacts themselves are.
    """
    return bids(
        root=config["output_dir"],
        datatype="atlasreg",
        session="pre",
        space="T1w",
        atlas="{atlas}",
        desc="atlas",
        suffix="labels",
        extension=extension,
        **inputs["post_ct"].wildcards,
    )


def get_atlas_lut(wildcards):
    """Lookup table for ``wildcards.atlas``.

    A surface atlas names its own parcels in the GIFTI label table, so the
    table generated when it was mapped into the volume is used rather than the
    config entry: that keeps the names right for parcellations that are not
    freesurfer's Desikan-Killiany set.

    Atlases that declare ``FreeSurferColorLUT`` otherwise get an empty input,
    and the lookup script falls back to the aseg/aparc table it carries
    internally, so no freesurfer installation is needed to resolve names.
    """
    if wildcards.atlas in get_surface_atlases():
        return expand(
            get_surface_atlas_lut(), atlas=wildcards.atlas, allow_missing=True
        )

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


def get_fsaverage_label_gii():
    """Atlas labels on the fsaverage surface, as a gifti.

    Published surface atlases are usually distributed as freesurfer
    annotations, which are converted once and shared across subjects.
    """
    return bids(
        root=config["output_dir"],
        datatype="atlasreg",
        space="fsaverage",
        hemi="{hemi}",
        atlas="{atlas}",
        suffix="dseg",
        extension=".label.gii",
    )


def get_fsaverage_surf_gii(surfname="{surfname}"):
    """An fsaverage surface, as a gifti."""
    return bids(
        root=config["output_dir"],
        datatype="atlasreg",
        space="fsaverage",
        hemi="{hemi}",
        suffix=f"{surfname}.surf.gii",
    )


def get_colortable(extension=".ctbl"):
    """Viewer colour table written beside an atlas segmentation."""
    return bids(
        root=config["output_dir"],
        datatype="atlasreg",
        session="pre",
        space="T1w",
        atlas="{atlas}",
        desc="colors",
        suffix="dseg",
        extension=extension,
        **inputs["pre_t1w"].wildcards,
    )


def get_colortable_atlases():
    """Atlases whose labels come with a lookup table we can colour from.

    aparc+aseg is excluded: it uses freesurfer's own indices, and both Slicer
    and freesurfer already ship FreeSurferColorLUT.txt for it.
    """
    return [
        atlas
        for atlas in config["atlas"]
        if atlas in get_surface_atlases()
        or atlas_lib.get_atlas_entry(config, atlas)["lut"]
        != "FreeSurferColorLUT"
    ]


def atlas_ships_annot(atlas):
    """True when the packaged atlas labels are a freesurfer annotation."""
    entry = atlas_lib.get_atlas_entry(config, atlas)

    return str(entry.get("label", "")).endswith(".annot")


def get_atlas_label_source():
    """Input function for the labels ``resample_atlas_to_subject`` resamples.

    An atlas shipped as a gifti is read where it is; one shipped as an
    annotation goes through the conversion rule first.
    """

    def _get_labels(wildcards):
        if atlas_ships_annot(wildcards.atlas):
            return expand(
                get_fsaverage_label_gii(),
                atlas=wildcards.atlas,
                hemi=wildcards.hemi,
            )

        return atlas_lib.get_atlas_file(
            workflow.basedir,
            config,
            wildcards.atlas,
            "label",
            hemi=wildcards.hemi,
        )

    return _get_labels


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
    final = inputs["post_ct"].expand(
        bids(
            root=config["output_dir"],
            datatype="ieeg",
            session="post",
            space="T1w",
            desc="atlas",
            suffix="electrodes",
            extension=".tsv",
            **inputs["post_ct"].wildcards,
        )
    )
    final.extend(
        inputs["post_ct"].expand(
            get_atlas_labels(extension=".xlsx"),
            atlas=config["atlas"],
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
    if get_colortable_atlases():
        for extension in (".ctbl", ".txt"):
            final.extend(
                inputs["pre_t1w"].expand(
                    get_colortable(extension=extension),
                    atlas=get_colortable_atlases(),
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
