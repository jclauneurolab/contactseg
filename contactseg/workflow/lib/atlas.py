"""Helpers for locating atlas resources and external derivatives.

These are kept out of the ``.smk`` files so they can be unit tested and reused
by the workflow scripts.
"""

from pathlib import Path

# freesurfer names hemispheres lh/rh, BIDS (and workbench) use L/R
HEMI_TO_FS = {"L": "lh", "R": "rh"}

# candidate sub-directories searched inside a derivatives dataset, in order
DERIV_SUBDIRS = ("", "ses-{session}", "ses-pre")


def get_atlas_dir(basedir, config):
    """Return the directory holding the packaged atlas resources.

    Parameters
    ----------
    basedir : str
        ``workflow.basedir`` of the running workflow.
    config : dict
        Workflow configuration.

    Returns
    -------
    pathlib.Path
    """

    return Path(basedir).parent / config["atlas_dir"]


def get_atlas_entry(config, atlas):
    """Return the configuration block describing ``atlas``."""

    try:
        return config["atlases"][atlas]
    except KeyError as err:
        raise ValueError(
            f"atlas '{atlas}' is not defined in the atlases config block"
        ) from err


def get_atlas_file(basedir, config, atlas, key, hemi=None):
    """Return the packaged file for ``key`` of ``atlas``.

    Patterns may contain ``{hemi}`` and ``{hemi_fs}`` placeholders, which are
    filled in from ``hemi``.
    """

    entry = get_atlas_entry(config, atlas)
    if key not in entry:
        raise ValueError(f"atlas '{atlas}' does not define a '{key}' file")

    return str(get_atlas_dir(basedir, config) / _fill_hemi(entry[key], hemi))


def get_derivatives_dir(config, key):
    """Return a derivatives directory from config, or raise if unset.

    Parameters
    ----------
    config : dict
        Workflow configuration.
    key : str
        Either ``smriprep_dir`` or ``freesurfer_dir``.
    """

    deriv_dir = config.get(key)
    if not deriv_dir:
        raise ValueError(f"--{key} must be set to use this part of the workflow")

    return Path(deriv_dir)


def get_derivatives_anat_spec(config, source):
    """Return ``(root, path)`` for the anatomical ``source`` was computed on.

    The config entry is either a plain path, resolved inside the derivatives
    dataset, or a mapping that also names the root it belongs to::

        derivatives_anat:
          smriprep: "anat/sub-{subject}_desc-preproc_T1w.nii.gz"
          freesurfer:
            root: bids
            path: "anat/sub-{subject}_ses-{session}_run-01_T1w.nii.gz"

    Either way the path is relative to the subject (or subject/session)
    directory of that root, so the two forms read the same way.

    Parameters
    ----------
    config : dict
        Workflow configuration.
    source : str
        Either ``smriprep`` or ``freesurfer``.

    Returns
    -------
    tuple of (str, str)
        The root name (``bids`` or ``derivatives``) and the path template.
    """

    try:
        entry = config["derivatives_anat"][source]
    except KeyError as err:
        raise ValueError(
            f"no derivatives_anat entry for '{source}'; add one naming the "
            "image that dataset was computed on"
        ) from err

    if isinstance(entry, str):
        return "derivatives", entry

    root = entry.get("root", "derivatives")
    if root not in ("bids", "derivatives"):
        raise ValueError(
            f"derivatives_anat[{source}].root must be 'bids' or "
            f"'derivatives', not '{root}'"
        )

    return root, entry["path"]


def find_subject_file(root, subject, relpath, session=None):
    """Locate ``relpath`` for ``subject`` inside a derivatives dataset.

    Both ``sub-<subject>/`` and ``sub-<subject>/ses-<session>/`` layouts are
    searched, so the same helper works for session-less sMRIPrep output and for
    the session-wise layout used by longitudinal datasets. When nothing is
    found the session-less path is returned, so that snakemake reports a
    missing input rather than an empty string.

    Parameters
    ----------
    root : str or pathlib.Path
        Root of the derivatives dataset.
    subject : str
        Subject label, without the ``sub-`` prefix.
    relpath : str
        Path relative to the subject (or subject/session) directory.
    session : str, optional
        Session label to try first.

    Returns
    -------
    str
    """

    root = Path(root)
    default = root / f"sub-{subject}" / relpath

    for subdir in DERIV_SUBDIRS:
        if "{session}" in subdir and session is None:
            continue
        candidate = root / f"sub-{subject}" / subdir.format(session=session) / relpath
        if candidate.exists():
            return str(candidate)

    return str(default)


def get_freesurfer_dir(config, subject, session=None):
    """Return the freesurfer subject directory for ``subject``."""

    root = get_derivatives_dir(config, "freesurfer_dir")
    for name in (
        f"sub-{subject}_ses-{session}" if session else None,
        f"sub-{subject}",
        subject,
    ):
        if name is None:
            continue
        if (root / name).exists():
            return root / name

    return root / f"sub-{subject}"


def get_freesurfer_file(config, subject, relpath, session=None, hemi=None):
    """Return a file inside the freesurfer subject directory."""

    subject_dir = get_freesurfer_dir(config, subject, session=session)

    return str(subject_dir / _fill_hemi(relpath, hemi))


def _fill_hemi(pattern, hemi):
    """Fill ``{hemi}``/``{hemi_fs}`` placeholders in a file pattern."""

    if hemi is None:
        return pattern

    return pattern.format(hemi=hemi, hemi_fs=HEMI_TO_FS[hemi])
