"""Tests for the atlas resource helpers."""

import pytest

from contactseg.workflow.lib import atlas


@pytest.fixture
def config():
    return {
        "atlas_dir": "resources/atlases",
        "derivatives_session": "pre",
        "atlases": {
            "CerebrA": {
                "space": "MNI152NLin2009cSym",
                "type": "volume",
                "dseg": "tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii.gz",
                "lut": "tpl-MNI152NLin2009cSym_atlas-CerebrA_dseg.tsv",
            },
            "DKTatlas": {
                "space": "native",
                "type": "surface",
                "label": "label/{hemi_fs}.aparc.DKTatlas.annot",
                "lut": "FreeSurferColorLUT",
            },
        },
    }


def test_get_atlas_file_resolves_relative_to_the_workflow(config):
    resolved = atlas.get_atlas_file(
        "/opt/contactseg/workflow", config, "CerebrA", "lut"
    )

    assert resolved == (
        "/opt/contactseg/resources/atlases/"
        "tpl-MNI152NLin2009cSym_atlas-CerebrA_dseg.tsv"
    )


def test_hemi_placeholders_are_filled_with_freesurfer_names(config, tmp_path):
    (tmp_path / "sub-P001" / "label").mkdir(parents=True)
    config["freesurfer_dir"] = str(tmp_path)

    resolved = atlas.get_freesurfer_file(
        config, "P001", config["atlases"]["DKTatlas"]["label"], hemi="R"
    )

    assert resolved.endswith("sub-P001/label/rh.aparc.DKTatlas.annot")


def test_unknown_atlas_is_reported_by_name(config):
    with pytest.raises(ValueError, match="Glasser"):
        atlas.get_atlas_entry(config, "Glasser")


def test_missing_derivatives_dir_names_the_flag(config):
    config["smriprep_dir"] = False

    with pytest.raises(ValueError, match="smriprep_dir"):
        atlas.get_derivatives_dir(config, "smriprep_dir")


def test_session_layout_is_preferred_when_present(tmp_path):
    anat = tmp_path / "sub-P001" / "ses-pre" / "anat"
    anat.mkdir(parents=True)
    (anat / "xfm.h5").touch()

    resolved = atlas.find_subject_file(tmp_path, "P001", "anat/xfm.h5", session="pre")

    assert resolved == str(anat / "xfm.h5")


def test_sessionless_path_is_returned_when_nothing_exists(tmp_path):
    resolved = atlas.find_subject_file(tmp_path, "P001", "anat/xfm.h5", session="pre")

    assert resolved == str(tmp_path / "sub-P001" / "anat" / "xfm.h5")


def test_derivatives_anat_defaults_to_the_derivatives_root(config):
    config["derivatives_anat"] = {"freesurfer": "mri/orig.mgz"}

    assert atlas.get_derivatives_anat_spec(config, "freesurfer") == (
        "derivatives",
        "mri/orig.mgz",
    )


def test_derivatives_anat_can_point_into_the_bids_dataset(config):
    config["derivatives_anat"] = {
        "freesurfer": {
            "root": "bids",
            "path": "anat/sub-{subject}_ses-{session}_run-01_T1w.nii.gz",
        }
    }

    root, path = atlas.get_derivatives_anat_spec(config, "freesurfer")

    assert root == "bids"
    assert path.format(subject="P167", session="pre") == (
        "anat/sub-P167_ses-pre_run-01_T1w.nii.gz"
    )


def test_an_unknown_derivatives_anat_root_is_rejected(config):
    config["derivatives_anat"] = {
        "freesurfer": {"root": "somewhere", "path": "mri/orig.mgz"}
    }

    with pytest.raises(ValueError, match="bids"):
        atlas.get_derivatives_anat_spec(config, "freesurfer")


def test_a_missing_derivatives_anat_entry_says_what_to_add(config):
    config["derivatives_anat"] = {}

    with pytest.raises(ValueError, match="smriprep"):
        atlas.get_derivatives_anat_spec(config, "smriprep")
