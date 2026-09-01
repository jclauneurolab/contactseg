def get_fs_surface_file(wildcards):
    """Freesurfer surface backing ``wildcards.surfname``."""
    return atlas_lib.get_freesurfer_file(
        config,
        wildcards.subject,
        config["surface_files"][wildcards.surfname],
        session=config["derivatives_session"],
        hemi=wildcards.hemi,
    )


def get_fs_reference_vol(wildcards):
    """Volume whose header carries the c_ras offset of the surfaces."""
    return atlas_lib.get_freesurfer_file(
        config,
        wildcards.subject,
        config["surface_reference"],
        session=config["derivatives_session"],
    )


def get_packaged_annot(wildcards):
    """Packaged atlas annotation, on the fsaverage surface."""
    return atlas_lib.get_atlas_file(
        workflow.basedir, config, wildcards.atlas, "label", hemi=wildcards.hemi
    )


def get_fsaverage_surf(wildcards):
    """An fsaverage surface, from the freesurfer installation."""
    return atlas_lib.find_fsaverage_file(
        config,
        config["fsaverage_files"][wildcards.surfname],
        hemi=wildcards.hemi,
    )


rule fs_surf_to_gifti:
    input:
        surf=get_fs_surface_file,
        ref_vol=get_fs_reference_vol,
        xfm_ras=get_bridge_xfm("freesurfer"),
    output:
        surf_gii=get_surf_gii("{surfname}"),
    wildcard_constraints:
        surfname="|".join(config["surface_files"].keys()),
    group:
        "subj"
    conda:
        "../envs/analysis.yaml"
    params:
        structure=lambda wildcards: config["structure_types"][wildcards.hemi],
        # spheres stay in their own frame, anatomical surfaces are shifted
        # from freesurfer tkrRAS into the scanner RAS of the T1w
        apply_cras=lambda wildcards: wildcards.surfname != "sphere",
    script:
        "../scripts/fs_surf_to_gifti.py"


rule create_midthickness:
    input:
        white=expand(
            get_surf_gii("{surfname}"), surfname="white", allow_missing=True
        ),
        pial=expand(
            get_surf_gii("{surfname}"), surfname="pial", allow_missing=True
        ),
    output:
        midthickness=expand(
            get_surf_gii("{surfname}"),
            surfname="midthickness",
            allow_missing=True,
        ),
    group:
        "subj"
    conda:
        "../envs/surface.yaml"
    shell:
        "wb_command -surface-average {output.midthickness}"
        " -surf {input.white} -surf {input.pial}"


if get_atlases_by_kind("native", "surface"):

    rule annot_to_label_gii:
        input:
            annot=get_fs_annot,
        output:
            label_gii=get_label_gii(),
        wildcard_constraints:
            atlas="|".join(get_atlases_by_kind("native", "surface")),
        group:
            "subj"
        conda:
            "../envs/analysis.yaml"
        params:
            structure=lambda wildcards: config["structure_types"][
                wildcards.hemi
            ],
        script:
            "../scripts/annot_to_label_gii.py"


if get_atlases_by_kind("fsaverage", "surface"):

    rule fsaverage_surf_to_gifti:
        input:
            surf=get_fsaverage_surf,
        output:
            surf_gii=get_fsaverage_surf_gii(),
        wildcard_constraints:
            surfname="|".join(config["fsaverage_files"].keys()),
        group:
            "subj"
        conda:
            "../envs/analysis.yaml"
        params:
            structure=lambda wildcards: config["structure_types"][
                wildcards.hemi
            ],
            # fsaverage surfaces are only ever used as a resampling reference,
            # never overlaid on a subject volume, so they stay in their own
            # frame
            apply_cras=False,
        script:
            "../scripts/fs_surf_to_gifti.py"

    rule atlas_annot_to_label_gii:
        input:
            annot=get_packaged_annot,
        output:
            label_gii=get_fsaverage_label_gii(),
        wildcard_constraints:
            atlas="|".join(
                atlas
                for atlas in get_atlases_by_kind("fsaverage", "surface")
                if atlas_ships_annot(atlas)
            ),
        group:
            "subj"
        conda:
            "../envs/analysis.yaml"
        params:
            structure=lambda wildcards: config["structure_types"][
                wildcards.hemi
            ],
        script:
            "../scripts/annot_to_label_gii.py"

    rule resample_atlas_to_subject:
        input:
            label_gii=get_atlas_label_source(),
            atlas_sphere=expand(
                get_fsaverage_surf_gii(),
                surfname="sphere",
                allow_missing=True,
            ),
            atlas_area=expand(
                get_fsaverage_surf_gii(),
                surfname="pial",
                allow_missing=True,
            ),
            subject_sphere=expand(
                get_surf_gii("{surfname}"),
                surfname="sphere",
                allow_missing=True,
            ),
            subject_area=expand(
                get_surf_gii("{surfname}"),
                surfname="pial",
                allow_missing=True,
            ),
        output:
            label_gii=get_label_gii(),
        wildcard_constraints:
            atlas="|".join(get_atlases_by_kind("fsaverage", "surface")),
        group:
            "subj"
        conda:
            "../envs/surface.yaml"
        params:
            structure=lambda wildcards: config["structure_types"][
                wildcards.hemi
            ],
        shell:
            # ADAP_BARY_AREA with -area-surfs weights each parcel by the
            # cortical area it actually covers, so small parcels are not lost
            # where fsaverage and the subject differ in surface area. Plain
            # BARYCENTRIC ignores that, which matters for a parcellation this
            # fine.
            "wb_command -label-resample {input.label_gii} {input.atlas_sphere}"
            " {input.subject_sphere} ADAP_BARY_AREA {output.label_gii}"
            " -area-surfs {input.atlas_area} {input.subject_area} && "
            "wb_command -set-structure {output.label_gii} {params.structure}"


if get_surface_atlases():

    rule map_surface_atlas_to_volume:
        input:
            label_gii=expand(
                get_label_gii(), hemi=config["hemi"], allow_missing=True
            ),
            midthickness=expand(
                get_surf_gii("{surfname}"),
                surfname="midthickness",
                hemi=config["hemi"],
                allow_missing=True,
            ),
            white=expand(
                get_surf_gii("{surfname}"),
                surfname="white",
                hemi=config["hemi"],
                allow_missing=True,
            ),
            pial=expand(
                get_surf_gii("{surfname}"),
                surfname="pial",
                hemi=config["hemi"],
                allow_missing=True,
            ),
            ref_vol=rules.n4biascorr.output.corrected_t1w,
        output:
            atlas_dseg=get_atlas_dseg_in_native(),
            lut=get_surface_atlas_lut(),
        wildcard_constraints:
            atlas="|".join(get_surface_atlases()),
        group:
            "subj"
        conda:
            "../envs/surface.yaml"
        params:
            key_offsets=lambda wildcards: get_key_offsets(wildcards.atlas),
            hemis=config["hemi"],
        script:
            "../scripts/surface_atlas_to_volume.py"


if config["map_atlas_surfaces"] and get_atlases_by_kind(
    "template", "volume"
) + get_atlases_by_kind("native", "volume"):

    rule dseg_tsv_to_label_list:
        input:
            lut=get_atlas_lut,
        output:
            label_list=bids(
                root=config["output_dir"],
                datatype="atlasreg",
                atlas="{atlas}",
                desc="wb",
                suffix="labellist",
                extension=".txt",
            ),
        group:
            "subj"
        conda:
            "../envs/analysis.yaml"
        script:
            "../scripts/dseg_tsv_to_label_list.py"

    rule map_volume_atlas_to_surface:
        input:
            atlas_dseg=get_atlas_dseg_in_native(),
            midthickness=expand(
                get_surf_gii("{surfname}"),
                surfname="midthickness",
                allow_missing=True,
            ),
            label_list=rules.dseg_tsv_to_label_list.output.label_list,
        output:
            label_gii=get_label_gii(),
        wildcard_constraints:
            atlas="|".join(
                get_atlases_by_kind("template", "volume")
                + get_atlases_by_kind("native", "volume")
            ),
        group:
            "subj"
        conda:
            "../envs/surface.yaml"
        params:
            structure=lambda wildcards: config["structure_types"][
                wildcards.hemi
            ],
            metric_gii=lambda wildcards, output: f"{output.label_gii}.shape.gii",
        shell:
            "wb_command -volume-to-surface-mapping {input.atlas_dseg}"
            " {input.midthickness} {params.metric_gii} -enclosing && "
            "wb_command -metric-label-import {params.metric_gii}"
            " {input.label_list} {output.label_gii} && "
            "wb_command -set-structure {output.label_gii} {params.structure} && "
            "rm -f {params.metric_gii}"
