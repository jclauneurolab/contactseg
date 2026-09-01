"""End-to-end tests for the atlas labelling scripts.

These run on synthetic volumes and surfaces, so they need no BIDS dataset, no
GPU and no atlas resources. Tests that need ANTsPy or connectome-workbench skip
themselves when those are not importable, which makes this file useful both as
a quick check after applying the branch and as CI.

    pytest tests/test_atlas_scripts.py -v
"""

import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).parent.parent / "contactseg" / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

FCSV_HEADER = (
    "# Markups fiducial file version = 4.11\n"
    "# CoordinateSystem = RAS\n"
    "# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,"
    "associatedNodeID\n"
)


def write_fcsv(path, contacts):
    """Write a minimal Slicer fcsv holding ``(name, x, y, z)`` tuples."""
    rows = [
        f"vtkMRMLMarkupsFiducialNode_{i},{x},{y},{z},0,0,0,1,1,1,0,{name},,"
        for i, (name, x, y, z) in enumerate(contacts)
    ]
    path.write_text(FCSV_HEADER + "\n".join(rows) + "\n")

    return path


def blob(center, shape=(48, 48, 48), radius=10):
    """A textured sphere, as something a registration can lock onto."""
    affine = np.eye(4)
    affine[:3, 3] = [-s / 2 for s in shape]
    grid = np.stack(np.meshgrid(*[np.arange(s) for s in shape], indexing="ij"), -1)
    world = nib.affines.apply_affine(affine, grid.reshape(-1, 3)).reshape(shape + (3,))
    dist = np.linalg.norm(world - np.array(center), axis=-1)
    data = np.where(dist < radius, 1.0, 0.0)
    data += 0.15 * np.where(dist < radius / 2, 1.0, 0.0)

    return nib.Nifti1Image(data.astype(np.float32), affine)


def centroid(path, threshold=0.5):
    """World-space centroid of the supra-threshold voxels of an image."""
    img = nib.load(str(path))
    indices = np.argwhere(np.asarray(img.dataobj) > threshold)

    return nib.affines.apply_affine(img.affine, indices.mean(0))


@pytest.fixture
def atlas_volume(tmp_path):
    """A two-block segmentation, left block 1 and right block 2."""
    vol = np.zeros((60, 60, 60), dtype=np.int32)
    vol[10:28, 20:40, 20:40] = 1
    vol[32:50, 20:40, 20:40] = 2
    affine = np.eye(4)
    affine[:3, 3] = [-30, -30, -30]

    path = tmp_path / "atlas_dseg.nii.gz"
    nib.save(nib.Nifti1Image(vol, affine), str(path))

    return path


@pytest.fixture
def lut(tmp_path):
    path = tmp_path / "lut.tsv"
    pd.DataFrame(
        {"label": [1, 2], "name": ["Block", "Block"], "hemi": ["L", "R"]}
    ).to_csv(path, sep="\t", index=False)

    return path


@pytest.fixture
def coords(tmp_path):
    return write_fcsv(
        tmp_path / "coords.fcsv",
        [
            ("LA1", -12.0, 0.0, 0.0),  # well inside the left block
            ("LA2", -3.0, 0.0, 0.0),  # on the medial boundary
            ("RB1", 12.0, 0.0, 0.0),  # well inside the right block
            ("RB2", 0.0, 45.0, 0.0),  # outside every label
        ],
    )


# ---- labelling -------------------------------------------------------------


def test_labels_follow_the_hemisphere_of_the_contact(
    tmp_path, atlas_volume, lut, coords
):
    from lookup_atlas_labels import lookup_atlas_labels

    labels = lookup_atlas_labels(
        coords, atlas_volume, lut, tmp_path / "labels.tsv", 1.0, 1.0, 3.0, 10.0
    ).set_index("name")

    assert labels.loc["LA1", "structure"] == "Block (Left)"
    assert labels.loc["RB1", "structure"] == "Block (Right)"


def test_a_contact_outside_every_label_is_unknown(tmp_path, atlas_volume, lut, coords):
    from lookup_atlas_labels import lookup_atlas_labels

    labels = lookup_atlas_labels(
        coords, atlas_volume, lut, tmp_path / "labels.tsv", 1.0, 1.0, 3.0, 10.0
    ).set_index("name")

    assert labels.loc["RB2", "structure"] == "Unknown"
    assert labels.loc["RB2", "nearest_structure"] == "None"


def test_a_contact_on_a_boundary_reports_lower_confidence(
    tmp_path, atlas_volume, lut, coords
):
    """The whole point of the Gaussian vote: an edge contact should say so."""
    from lookup_atlas_labels import lookup_atlas_labels

    labels = lookup_atlas_labels(
        coords, atlas_volume, lut, tmp_path / "labels.tsv", 1.0, 1.0, 3.0, 10.0
    ).set_index("name")

    assert labels.loc["LA2", "probability"] < labels.loc["LA1", "probability"]
    assert labels.loc["LA2", "entropy"] > labels.loc["LA1", "entropy"]
    assert (
        labels.loc["LA2", "dist_to_boundary_mm"]
        < labels.loc["LA1", "dist_to_boundary_mm"]
    )


def test_a_wider_sigma_widens_the_distribution(tmp_path, atlas_volume, lut, coords):
    """Distrusting the registration has to cost confidence, not change it."""
    from lookup_atlas_labels import lookup_atlas_labels

    tight = lookup_atlas_labels(
        coords, atlas_volume, lut, tmp_path / "a.tsv", 1.0, 0.5, 3.0, 10.0
    ).set_index("name")
    loose = lookup_atlas_labels(
        coords, atlas_volume, lut, tmp_path / "b.tsv", 1.0, 4.0, 3.0, 10.0
    ).set_index("name")

    assert loose.loc["LA2", "probability"] < tight.loc["LA2", "probability"]


def test_freesurfer_indices_resolve_without_a_lookup_table(tmp_path, coords):
    """No lut file means the built-in aseg/aparc table, and tissue classes."""
    from lookup_atlas_labels import lookup_atlas_labels

    vol = np.zeros((60, 60, 60), dtype=np.int32)
    vol[10:28, 20:40, 20:40] = 1001  # ctx-lh-bankssts
    vol[32:50, 20:40, 20:40] = 41  # Right-Cerebral-White-Matter
    affine = np.eye(4)
    affine[:3, 3] = [-30, -30, -30]
    dseg = tmp_path / "aparcaseg.nii.gz"
    nib.save(nib.Nifti1Image(vol, affine), str(dseg))

    labels = lookup_atlas_labels(
        coords, dseg, None, tmp_path / "labels.tsv", 1.0, 1.0, 3.0, 10.0
    ).set_index("name")

    assert labels.loc["LA1", "structure"] == "ctx-lh-bankssts"
    assert labels.loc["LA1", "tissue"] == "GM"
    assert labels.loc["RB1", "tissue"] == "WM"
    assert labels.loc["LA1", "p_GM"] == pytest.approx(1.0)


# ---- tissue and merging ----------------------------------------------------


def test_tissue_probabilities_are_sampled_at_the_contact(tmp_path, coords):
    from lookup_tissue_labels import lookup_tissue_labels

    affine = np.eye(4)
    affine[:3, 3] = [-30, -30, -30]
    probseg = []
    for i, name in enumerate(["GM", "WM", "CSF"]):
        data = np.full((60, 60, 60), 0.1 * (i + 1), dtype=np.float32)
        path = tmp_path / f"label-{name}_probseg.nii.gz"
        nib.save(nib.Nifti1Image(data, affine), str(path))
        probseg.append(path)

    tissue = lookup_tissue_labels(
        coords, probseg, ["GM", "WM", "CSF"], tmp_path / "tissue.tsv"
    )

    # the first three contacts sit inside the maps, the fourth is outside them
    assert list(tissue["tissue"]) == ["CSF", "CSF", "CSF", "Unknown"]
    assert tissue["p_GM"].iloc[0] == pytest.approx(0.1)
    assert tissue["p_GM"].iloc[3] == 0.0


def test_merged_table_prefixes_columns_per_atlas(tmp_path, atlas_volume, lut, coords):
    from gen_atlas_electrodes import gen_atlas_electrodes
    from lookup_atlas_labels import lookup_atlas_labels

    labels = tmp_path / "labels.tsv"
    lookup_atlas_labels(coords, atlas_volume, lut, labels, 1.0, 1.0, 3.0, 10.0)

    merged = gen_atlas_electrodes(
        coords_fcsv=coords,
        template_coords_fcsv=coords,
        atlas_label_files=[labels],
        tissue_file=None,
        atlases=["CerebrA"],
        template="MNI152NLin2009cSym",
        output_tsv=tmp_path / "electrodes.tsv",
        output_fcsv=tmp_path / "electrodes.fcsv",
    )

    assert "CerebrA_structure" in merged.columns
    assert "MNI152NLin2009cSym_x" in merged.columns
    assert len(merged) == 4

    # the fcsv carries the structure in the label column, for Slicer
    written = (tmp_path / "electrodes.fcsv").read_text().splitlines()
    assert written[3].split(",")[11] == "Block (Left)"


# ---- freesurfer surfaces ---------------------------------------------------


@pytest.fixture
def fs_surface(tmp_path):
    """A four-vertex surface plus a conformed volume with a c_ras offset."""
    affine = np.eye(4)
    affine[:3, 3] = [-16 + 5, -16 - 3, -16 + 2]
    orig = tmp_path / "orig.mgz"
    nib.save(nib.MGHImage(np.zeros((32, 32, 32), dtype=np.uint8), affine), orig)

    vertices = np.array(
        [[0, 0, 0], [10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=np.float32
    )
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], np.int32)
    surf = tmp_path / "lh.white"
    nib.freesurfer.write_geometry(str(surf), vertices, faces)

    return surf, orig, vertices


def test_surfaces_are_shifted_out_of_tkr_into_scanner_ras(tmp_path, fs_surface):
    """Forgetting c_ras puts the surfaces a centimetre off the volume."""
    from fs_surf_to_gifti import fs_surf_to_gifti

    surf, orig, vertices = fs_surface
    out = tmp_path / "lh.white.surf.gii"
    fs_surf_to_gifti(str(surf), str(orig), str(out), "CORTEX_LEFT", True)

    header = nib.load(str(orig)).header
    expected = nib.affines.apply_affine(
        header.get_vox2ras() @ np.linalg.inv(header.get_vox2ras_tkr()),
        vertices,
    )
    gii = nib.load(str(out))

    assert np.allclose(gii.agg_data("NIFTI_INTENT_POINTSET"), expected)
    assert gii.meta["AnatomicalStructurePrimary"] == "CORTEX_LEFT"


def test_spheres_are_left_in_their_own_frame(tmp_path, fs_surface):
    from fs_surf_to_gifti import fs_surf_to_gifti

    surf, orig, vertices = fs_surface
    out = tmp_path / "lh.sphere.surf.gii"
    fs_surf_to_gifti(str(surf), str(orig), str(out), "CORTEX_LEFT", False)

    assert np.allclose(nib.load(str(out)).agg_data("NIFTI_INTENT_POINTSET"), vertices)


def test_annotations_keep_their_names_and_colours(tmp_path):
    from annot_to_label_gii import annot_to_label_gii

    labels = np.array([0, 1, 1, 2], dtype=np.int32)
    ctab = np.array(
        [[25, 5, 25, 0, 0], [100, 20, 30, 0, 0], [10, 200, 40, 0, 0]],
        dtype=np.int32,
    )
    ctab[:, 4] = ctab[:, 0] + ctab[:, 1] * 2**8 + ctab[:, 2] * 2**16
    annot = tmp_path / "lh.aparc.annot"
    nib.freesurfer.write_annot(
        str(annot),
        labels,
        ctab,
        ["unknown", "parcelA", "parcelB"],
        fill_ctab=False,
    )

    out = tmp_path / "lh.aparc.label.gii"
    annot_to_label_gii(str(annot), str(out), "CORTEX_LEFT")
    gii = nib.load(str(out))

    assert list(gii.agg_data()) == [0, 1, 1, 2]
    assert [label.label for label in gii.labeltable.labels] == [
        "unknown",
        "parcelA",
        "parcelB",
    ]


def test_hemispheres_are_offset_apart_before_merging(tmp_path):
    """Two hemispheres numbered 1..n must not collide in the merged volume."""
    from surface_atlas_to_volume import merge_hemispheres

    affine = np.eye(4)
    hemis = []
    for i in range(2):
        data = np.zeros((10, 10, 10), dtype=np.int32)
        data[i * 5 : i * 5 + 4] = 3
        path = tmp_path / f"hemi{i}.nii.gz"
        nib.save(nib.Nifti1Image(data, affine), str(path))
        hemis.append(str(path))

    out = tmp_path / "merged.nii.gz"
    merge_hemispheres(hemis, [1000, 2000], str(out))
    merged = np.asarray(nib.load(str(out)).dataobj)

    assert set(np.unique(merged)) == {0, 1003, 2003}


# ---- exports ---------------------------------------------------------------


def test_segmentation_exports_named_segments_for_slicer(tmp_path, atlas_volume, lut):
    nrrd = pytest.importorskip("nrrd")
    from atlas_dseg_to_nrrd import atlas_dseg_to_nrrd

    out = tmp_path / "atlas.seg.nrrd"
    atlas_dseg_to_nrrd(str(atlas_volume), str(lut), str(out), "ras")
    _, header = nrrd.read(str(out))

    assert header["space"] == "right-anterior-superior"
    assert header["Segment0_Name"] == "L_Block"
    assert header["Segment1_Name"] == "R_Block"


def test_label_list_is_written_in_workbench_pairs(tmp_path, lut):
    from dseg_tsv_to_label_list import dseg_tsv_to_label_list

    out = tmp_path / "labellist.txt"
    dseg_tsv_to_label_list(str(lut), str(out))
    lines = out.read_text().splitlines()

    assert lines[0] == "Block_L"
    assert lines[1].split()[0] == "1"
    assert len(lines[1].split()) == 5  # key r g b a


def test_freesurfer_volumes_are_imported_without_resampling(tmp_path, atlas_volume):
    from import_fs_atlas import import_fs_atlas

    out = tmp_path / "imported.nii.gz"
    import_fs_atlas(str(atlas_volume), str(out))

    source, imported = nib.load(str(atlas_volume)), nib.load(str(out))
    assert np.allclose(source.affine, imported.affine)
    assert np.array_equal(np.asarray(source.dataobj), np.asarray(imported.dataobj))


# ---- registration ----------------------------------------------------------


@pytest.fixture
def registered_pair(tmp_path):
    """Register a blob at +10mm to the same blob at the origin."""
    pytest.importorskip("ants")
    from template_registration import template_registration

    shift = np.array([10.0, 0.0, 0.0])
    nib.save(blob(shift), str(tmp_path / "t1w.nii.gz"))
    nib.save(blob([0, 0, 0]), str(tmp_path / "template.nii.gz"))

    template_registration(
        str(tmp_path / "t1w.nii.gz"),
        str(tmp_path / "template.nii.gz"),
        str(tmp_path / "affine.mat"),
        str(tmp_path / "warp.nii.gz"),
        str(tmp_path / "invwarp.nii.gz"),
        str(tmp_path / "affine.txt"),
        str(tmp_path / "t1w_in_template.nii.gz"),
    )

    return tmp_path, shift


def test_atlas_is_pulled_onto_the_subject(registered_pair):
    """The image direction: template space in, subject space out."""
    from warp_atlas_to_native import warp_atlas_to_native

    tmp_path, shift = registered_pair
    template = nib.load(str(tmp_path / "template.nii.gz"))
    labels = np.where(np.asarray(template.dataobj) > 0.5, 5, 0).astype(np.int32)
    nib.save(
        nib.Nifti1Image(labels, template.affine),
        str(tmp_path / "template_atlas.nii.gz"),
    )

    out = tmp_path / "atlas_in_native.nii.gz"
    warp_atlas_to_native(
        str(tmp_path / "template_atlas.nii.gz"),
        str(tmp_path / "t1w.nii.gz"),
        [str(tmp_path / "affine.mat"), str(tmp_path / "invwarp.nii.gz")],
        [True, False],
        str(out),
    )

    assert np.allclose(centroid(out), shift, atol=1.0)


def test_contacts_are_pushed_into_template_space(registered_pair):
    """The point direction, which is the opposite of the image direction."""
    from warp_coords import warp_coords

    tmp_path, shift = registered_pair
    native = write_fcsv(
        tmp_path / "native.fcsv",
        [
            ("C1", *shift),
            ("C2", shift[0] + 5, 4.0, -3.0),
        ],
    )

    out = tmp_path / "template.fcsv"
    warp_coords(
        str(native),
        str(out),
        [str(tmp_path / "affine.mat"), str(tmp_path / "invwarp.nii.gz")],
        [True, False],
    )
    warped = pd.read_csv(out, skiprows=3, header=None)[[1, 2, 3]].to_numpy()

    assert np.allclose(warped[0], [0, 0, 0], atol=1.0)
    assert np.allclose(warped[1], [5, 4, -3], atol=1.0)


# ---- derivatives bridge ----------------------------------------------------


def rigid_ras(angle_deg=20.0, translation=(3.0, 6.0, -2.0)):
    """A 4x4 RAS rigid transform, standing in for the measured bridge."""
    theta = np.deg2rad(angle_deg)
    matrix = np.eye(4)
    matrix[:3, :3] = [
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1],
    ]
    matrix[:3, 3] = translation

    return matrix


@pytest.fixture
def bridge(tmp_path):
    """A bridge matrix on disk, in the 4x4 RAS convention used by main."""
    matrix = rigid_ras()
    path = tmp_path / "from-freesurfer_to-T1w_xfm.txt"
    np.savetxt(path, matrix)

    return path, matrix


def test_bridge_moves_the_segmentation_without_touching_its_labels(
    tmp_path, atlas_volume, bridge
):
    """A rigid composes exactly with the affine, so no voxel is resampled."""
    from import_fs_atlas import import_fs_atlas

    xfm_path, matrix = bridge
    out = tmp_path / "bridged.nii.gz"
    import_fs_atlas(str(atlas_volume), str(out), xfm_ras=str(xfm_path))

    source, bridged = nib.load(str(atlas_volume)), nib.load(str(out))
    assert np.array_equal(np.asarray(source.dataobj), np.asarray(bridged.dataobj))
    assert np.allclose(bridged.affine, matrix @ source.affine)


def test_labels_survive_a_derivatives_frame_offset(
    tmp_path, atlas_volume, lut, coords, bridge
):
    """The point of the bridge: same labels whether or not the frames differ.

    The segmentation is put into a deliberately offset frame, as it would be if
    freesurfer had run on a different acquisition, and the bridge is applied.
    Every contact has to come back with the label it had when the frames
    already agreed.
    """
    from import_fs_atlas import import_fs_atlas
    from lookup_atlas_labels import lookup_atlas_labels

    aligned = lookup_atlas_labels(
        coords, atlas_volume, lut, tmp_path / "aligned.tsv", 1.0, 1.0, 3.0, 10.0
    )

    # push the segmentation into the derivatives frame
    xfm_path, matrix = bridge
    source = nib.load(str(atlas_volume))
    offset = tmp_path / "offset_dseg.nii.gz"
    nib.save(
        nib.Nifti1Image(
            np.asarray(source.dataobj), np.linalg.inv(matrix) @ source.affine
        ),
        str(offset),
    )

    bridged = tmp_path / "bridged_dseg.nii.gz"
    import_fs_atlas(str(offset), str(bridged), xfm_ras=str(xfm_path))
    recovered = lookup_atlas_labels(
        coords, bridged, lut, tmp_path / "bridged.tsv", 1.0, 1.0, 3.0, 10.0
    )

    assert list(recovered["structure"]) == list(aligned["structure"])

    # and without the bridge the labels are wrong, which is the failure mode
    unbridged = lookup_atlas_labels(
        coords, offset, lut, tmp_path / "unbridged.tsv", 1.0, 1.0, 3.0, 10.0
    )
    assert list(unbridged["structure"]) != list(aligned["structure"])


def test_bridge_is_applied_to_surfaces_after_the_cras_shift(
    tmp_path, fs_surface, bridge
):
    from fs_surf_to_gifti import fs_surf_to_gifti

    surf, orig, vertices = fs_surface
    xfm_path, matrix = bridge
    out = tmp_path / "lh.white.surf.gii"
    fs_surf_to_gifti(
        str(surf),
        str(orig),
        str(out),
        "CORTEX_LEFT",
        True,
        xfm_ras=str(xfm_path),
    )

    header = nib.load(str(orig)).header
    cras = header.get_vox2ras() @ np.linalg.inv(header.get_vox2ras_tkr())
    expected = nib.affines.apply_affine(matrix @ cras, vertices)

    assert np.allclose(
        nib.load(str(out)).agg_data("NIFTI_INTENT_POINTSET"),
        expected,
        atol=1e-4,
    )


def test_tissue_maps_are_sampled_in_the_derivatives_frame(tmp_path, coords, bridge):
    """Contacts are pulled back through the bridge before sampling."""
    from lookup_tissue_labels import lookup_tissue_labels

    xfm_path, matrix = bridge
    affine = np.eye(4)
    affine[:3, 3] = [-30, -30, -30]

    # a ramp along x, so a mis-sampled contact reads a different value
    ramp = np.tile(np.linspace(0, 1, 60, dtype=np.float32)[:, None, None], (1, 60, 60))
    probseg = []
    for name in ["GM", "WM", "CSF"]:
        path = tmp_path / f"label-{name}_probseg.nii.gz"
        nib.save(nib.Nifti1Image(ramp, affine), str(path))
        probseg.append(path)

    bridged = lookup_tissue_labels(
        coords,
        probseg,
        ["GM", "WM", "CSF"],
        tmp_path / "bridged.tsv",
        xfm_ras=str(xfm_path),
    )
    plain = lookup_tissue_labels(
        coords, probseg, ["GM", "WM", "CSF"], tmp_path / "plain.tsv"
    )

    # the contact is pulled back by the inverse bridge before sampling
    contact = np.array([-12.0, 0.0, 0.0])
    expected_x = nib.affines.apply_affine(np.linalg.inv(matrix), contact)[0]
    assert bridged["p_GM"].iloc[0] == pytest.approx((expected_x + 30) / 59, abs=0.02)
    assert bridged["p_GM"].iloc[0] != pytest.approx(plain["p_GM"].iloc[0])


def test_surface_atlas_names_come_from_the_gifti_not_a_hard_coded_list(
    tmp_path,
):
    """A parcellation that is not Desikan-Killiany still names its regions.

    The volume labels are the GIFTI keys plus a per-hemisphere offset, so the
    lookup table has to be built from the same two things rather than assumed.
    """
    from annot_to_label_gii import annot_to_label_gii
    from surface_atlas_to_volume import write_lookup_table

    label_files = []
    for hemi, names in (
        ("L", ["unknown", "G_temp_sup-Lateral", "S_calcarine"]),
        ("R", ["unknown", "G_temp_sup-Lateral", "S_calcarine"]),
    ):
        labels = np.array([0, 1, 2, 1], dtype=np.int32)
        ctab = np.array(
            [[25, 5, 25, 0, 0], [100, 20, 30, 0, 0], [10, 200, 40, 0, 0]],
            dtype=np.int32,
        )
        ctab[:, 4] = ctab[:, 0] + ctab[:, 1] * 2**8 + ctab[:, 2] * 2**16
        annot = tmp_path / f"{hemi}.annot"
        nib.freesurfer.write_annot(str(annot), labels, ctab, names, fill_ctab=False)
        gii = tmp_path / f"{hemi}.label.gii"
        annot_to_label_gii(str(annot), str(gii), f"CORTEX_{hemi}")
        label_files.append(str(gii))

    out = tmp_path / "dseg.tsv"
    write_lookup_table(label_files, [1000, 2000], ["L", "R"], str(out))
    table = pd.read_csv(out, sep="\t")

    # keys are offset per hemisphere, names are the atlas's own
    assert set(table["label"]) == {1001, 1002, 2001, 2002}
    assert table.loc[table["label"] == 2002, "name"].iloc[0] == "S_calcarine"
    assert list(table.loc[table["label"] < 2000, "hemi"].unique()) == ["L"]
    # background is not a region
    assert 1000 not in set(table["label"])


def test_colour_tables_are_written_for_both_viewers(tmp_path):
    """Slicer and ITK-SNAP take the same information, punctuated differently."""
    from gen_colortable import gen_colortable

    lut = tmp_path / "lut.tsv"
    pd.DataFrame(
        {
            "label": [1, 1002],
            "name": ["G temp sup", "S calcarine"],
            "hemi": ["L", "R"],
            "r": [220, 10],
            "g": [20, 200],
            "b": [30, 40],
        }
    ).to_csv(lut, sep="\t", index=False)

    ctbl = tmp_path / "colors.ctbl"
    itksnap = tmp_path / "colors.txt"
    gen_colortable(str(lut), "Yale", str(ctbl), str(itksnap))

    slicer = ctbl.read_text().splitlines()
    assert slicer[0].startswith("# Color table file Yale")
    assert slicer[2] == "0 background 0 0 0 0"
    # index, name, then the atlas's own colour
    assert slicer[3] == "1 L_G_temp_sup 220 20 30 255"
    assert slicer[4].startswith("1002 R_S_calcarine 10 200 40")

    snap = itksnap.read_text()
    assert snap.startswith("# ITK-SnAP Label Description File")
    assert '"L_G_temp_sup"' in snap
    assert "1002" in snap


def test_colours_fall_back_to_a_stable_palette(tmp_path):
    """An atlas whose table carries no colours still gets distinct ones."""
    from gen_colortable import gen_colortable

    lut = tmp_path / "lut.tsv"
    pd.DataFrame({"label": [1, 2], "name": ["a", "b"]}).to_csv(
        lut, sep="\t", index=False
    )

    ctbl = tmp_path / "colors.ctbl"
    gen_colortable(str(lut), "CerebrA", str(ctbl), str(tmp_path / "c.txt"))
    rows = [line.split() for line in ctbl.read_text().splitlines()[3:] if line]

    assert rows[0][2:5] != rows[1][2:5]
    # and the same index always gets the same colour
    gen_colortable(str(lut), "CerebrA", str(ctbl), str(tmp_path / "c.txt"))
    assert [line.split() for line in ctbl.read_text().splitlines()[3:] if line] == rows


def test_a_sphere_ignores_the_bridge_even_when_given_one(tmp_path, fs_surface, bridge):
    """A sphere is a registration space, not an anatomical one.

    The rule hands every surface the bridge transform, so the guard that keeps
    it off the spheres is the only thing preventing a subject sphere from being
    rotated out of correspondence with fsaverage -- which would misresample
    every parcel without any obvious symptom.
    """
    from fs_surf_to_gifti import fs_surf_to_gifti

    surf, orig, vertices = fs_surface
    xfm_path, _ = bridge
    out = tmp_path / "lh.sphere.surf.gii"
    fs_surf_to_gifti(
        str(surf),
        str(orig),
        str(out),
        "CORTEX_LEFT",
        apply_cras=False,
        xfm_ras=str(xfm_path),
    )

    assert np.allclose(nib.load(str(out)).agg_data("NIFTI_INTENT_POINTSET"), vertices)


# ---- snakemake glue --------------------------------------------------------


class FakeIO(dict):
    """Stands in for snakemake's Namedlist: attribute and .get() access."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as err:
            raise AttributeError(name) from err


def run_script(name, **sections):
    """Execute a workflow script's __main__ block with a fake snakemake.

    A dry run never executes a script, so the block that unpacks the snakemake
    object is the one part of a rule that no DAG test can reach. Scripts shared
    by two rules with different inputs are where that bites: reading an input
    the other rule does not declare fails only at runtime.
    """
    import runpy
    from types import SimpleNamespace

    fake = SimpleNamespace(**{key: FakeIO(value) for key, value in sections.items()})
    runpy.run_path(
        str(SCRIPTS / f"{name}.py"),
        init_globals={"snakemake": fake},
        run_name="__main__",
    )


def test_subject_surface_glue_runs(tmp_path, fs_surface, bridge):
    surf, orig, _ = fs_surface
    xfm_path, _ = bridge
    out = tmp_path / "sub.surf.gii"

    run_script(
        "fs_surf_to_gifti",
        input={
            "surf": str(surf),
            "ref_vol": str(orig),
            "xfm_ras": str(xfm_path),
        },
        output={"surf_gii": str(out)},
        params={"structure": "CORTEX_LEFT", "apply_cras": True},
    )

    assert out.exists()


def test_fsaverage_surface_glue_runs_without_ref_vol_or_bridge(tmp_path, fs_surface):
    """The fsaverage rule declares neither, and must not fail on their absence."""
    surf, _, vertices = fs_surface
    out = tmp_path / "fsaverage.surf.gii"

    run_script(
        "fs_surf_to_gifti",
        input={"surf": str(surf)},
        output={"surf_gii": str(out)},
        params={"structure": "CORTEX_LEFT", "apply_cras": False},
    )

    assert np.allclose(nib.load(str(out)).agg_data("NIFTI_INTENT_POINTSET"), vertices)


def test_annot_glue_runs(tmp_path):
    from annot_to_label_gii import annot_to_label_gii  # noqa: F401

    labels = np.array([0, 1, 1, 2], dtype=np.int32)
    ctab = np.array(
        [[25, 5, 25, 0, 0], [100, 20, 30, 0, 0], [10, 200, 40, 0, 0]],
        dtype=np.int32,
    )
    ctab[:, 4] = ctab[:, 0] + ctab[:, 1] * 2**8 + ctab[:, 2] * 2**16
    annot = tmp_path / "lh.annot"
    nib.freesurfer.write_annot(
        str(annot), labels, ctab, ["unknown", "a", "b"], fill_ctab=False
    )

    out = tmp_path / "lh.label.gii"
    run_script(
        "annot_to_label_gii",
        input={"annot": str(annot)},
        output={"label_gii": str(out)},
        params={"structure": "CORTEX_LEFT"},
    )

    assert out.exists()


def test_summary_reports_what_it_can_when_a_column_is_absent(capsys):
    """The concordance line is a diagnostic; losing it must not lose labels."""
    from lookup_atlas_labels import report_concordance

    report_concordance(pd.DataFrame({"structure": ["a", "Unknown"]}), "x.nii")
    printed = capsys.readouterr().out

    assert "labelled 1/2" in printed
    # no name column, so the hemisphere half is simply not reported
    assert "hemisphere match" not in printed


def test_hemisphere_agreement_is_reported_without_a_column(capsys):
    """The frame check survives the column being dropped from the table."""
    from lookup_atlas_labels import report_concordance

    labels = pd.DataFrame(
        {
            "name": ["LA1", "LA2", "RB1", "RB2"],
            "structure": [
                "ctx-lh-insula",
                "ctx-rh-insula",  # left contact, right structure
                "ctx-rh-insula",
                "ctx-rh-insula",
            ],
        }
    )
    report_concordance(labels, "x.nii")
    printed = capsys.readouterr().out

    assert "hemisphere match 3/4 (75%)" in printed
    assert "low concordance" in printed
    assert "hemi_match" not in list(labels.columns)


def test_labels_are_written_before_the_summary(tmp_path, atlas_volume, lut, coords):
    from lookup_atlas_labels import lookup_atlas_labels

    out = tmp_path / "labels.tsv"
    lookup_atlas_labels(coords, atlas_volume, lut, out, 1.0, 1.0, 3.0, 10.0)

    assert out.exists()
    assert len(pd.read_csv(out, sep="\t")) == 4
