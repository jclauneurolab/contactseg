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
