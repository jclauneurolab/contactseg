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
    assert labels.loc[["LA1", "RB1"], "hemi_match"].all()


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

    assert labels.loc["LA2", "p_top"] < labels.loc["LA1", "p_top"]
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

    assert loose.loc["LA2", "p_top"] < tight.loc["LA2", "p_top"]


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
