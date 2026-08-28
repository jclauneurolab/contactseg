"""Label contacts from an atlas segmentation, with an uncertainty estimate.

A hard nearest-voxel lookup throws away the fact that a contact sits in a
finite volume of tissue, and that the coordinates and the segmentation each
carry error. Instead of a single label, the neighbourhood of every contact is
sampled with a Gaussian weight and reported as a probability distribution over
structures. The hard label is still reported, so the output is a superset of
the classic lookup.

The Gaussian sigma combines two independent sources of error::

    sigma_total = sqrt(sigma_contact**2 + sigma_reg**2)

``sigma_contact`` is the physical contact extent plus the localisation error,
``sigma_reg`` the residual error of the registration that brought the atlas and
the contacts into a common space. Raising ``sigma_reg`` propagates a less
trusted registration into wider, less confident distributions rather than
silently reporting a crisp wrong answer.

Labels are never interpolated: the lookup happens in world coordinates, in the
segmentation's own grid, whatever that grid is.
"""

import nibabel as nib
import numpy as np
import pandas as pd

FCSV_HEADER_ROWS = 3

# structures freesurfer writes into aseg, used when no lookup table is given
ASEG = {
    0: "Unknown",
    2: "Left-Cerebral-White-Matter",
    3: "Left-Cerebral-Cortex",
    4: "Left-Lateral-Ventricle",
    5: "Left-Inf-Lat-Vent",
    7: "Left-Cerebellum-White-Matter",
    8: "Left-Cerebellum-Cortex",
    10: "Left-Thalamus",
    11: "Left-Caudate",
    12: "Left-Putamen",
    13: "Left-Pallidum",
    14: "3rd-Ventricle",
    15: "4th-Ventricle",
    16: "Brain-Stem",
    17: "Left-Hippocampus",
    18: "Left-Amygdala",
    24: "CSF",
    26: "Left-Accumbens-area",
    28: "Left-VentralDC",
    30: "Left-vessel",
    31: "Left-choroid-plexus",
    41: "Right-Cerebral-White-Matter",
    42: "Right-Cerebral-Cortex",
    43: "Right-Lateral-Ventricle",
    44: "Right-Inf-Lat-Vent",
    46: "Right-Cerebellum-White-Matter",
    47: "Right-Cerebellum-Cortex",
    49: "Right-Thalamus",
    50: "Right-Caudate",
    51: "Right-Putamen",
    52: "Right-Pallidum",
    53: "Right-Hippocampus",
    54: "Right-Amygdala",
    58: "Right-Accumbens-area",
    60: "Right-VentralDC",
    62: "Right-vessel",
    63: "Right-choroid-plexus",
    72: "5th-Ventricle",
    77: "WM-hypointensities",
    78: "Left-WM-hypointensities",
    79: "Right-WM-hypointensities",
    80: "non-WM-hypointensities",
    85: "Optic-Chiasm",
    251: "CC_Posterior",
    252: "CC_Mid_Posterior",
    253: "CC_Central",
    254: "CC_Mid_Anterior",
    255: "CC_Anterior",
}

DK = [
    "unknown",
    "bankssts",
    "caudalanteriorcingulate",
    "caudalmiddlefrontal",
    "corpuscallosum",
    "cuneus",
    "entorhinal",
    "fusiform",
    "inferiorparietal",
    "inferiortemporal",
    "isthmuscingulate",
    "lateraloccipital",
    "lateralorbitofrontal",
    "lingual",
    "medialorbitofrontal",
    "middletemporal",
    "parahippocampal",
    "paracentral",
    "parsopercularis",
    "parsorbitalis",
    "parstriangularis",
    "pericalcarine",
    "postcentral",
    "posteriorcingulate",
    "precentral",
    "precuneus",
    "rostralanteriorcingulate",
    "rostralmiddlefrontal",
    "superiorfrontal",
    "superiorparietal",
    "superiortemporal",
    "supramarginal",
    "frontalpole",
    "temporalpole",
    "transversetemporal",
    "insula",
]

# label indices grouped into coarse tissue classes, following the freesurfer
# convention. Only used for segmentations that carry freesurfer indices.
FS_SUBCORTICAL_GM = (
    10,
    11,
    12,
    13,
    17,
    18,
    26,
    28,
    49,
    50,
    51,
    52,
    53,
    54,
    58,
    60,
)
FS_WM = (2, 41, 7, 46, 77, 78, 79, 5001, 5002)
FS_CSF = (4, 5, 14, 15, 24, 43, 44, 72, 31, 63)


def build_fs_lut():
    """Return the aseg + aparc label table freesurfer uses by default."""
    lut = dict(ASEG)
    for i, name in enumerate(DK):
        lut[1000 + i] = f"ctx-lh-{name}"
        lut[2000 + i] = f"ctx-rh-{name}"
        lut[3000 + i] = f"wm-lh-{name}"
        lut[4000 + i] = f"wm-rh-{name}"
    lut[5001] = "Left-UnsegmentedWhiteMatter"
    lut[5002] = "Right-UnsegmentedWhiteMatter"

    return lut


def first_or_none(value):
    """Return a single path from an input that may be an empty list."""
    if isinstance(value, (list, tuple)):
        return value[0] if value else None

    return value or None


def read_lut(lut_file):
    """Read a BIDS dseg.tsv lookup table into ``{index: name}``.

    Names are suffixed with the hemisphere when the table carries one, so that
    a left and a right structure sharing a name stay distinguishable.
    """
    if not lut_file:
        return build_fs_lut(), True

    table = pd.read_csv(lut_file, sep="\t")
    table.columns = table.columns.str.lower()

    lut = {}
    for _, row in table.iterrows():
        name = str(row["name"])
        if "hemi" in table.columns and str(row["hemi"]) in ("L", "R"):
            side = "Left" if row["hemi"] == "L" else "Right"
            name = f"{name} ({side})"
        lut[int(row["label"])] = name

    # background is rarely listed in a dseg.tsv, but every unlabelled contact
    # lands on it
    lut.setdefault(0, "Unknown")

    return lut, False


def tissue_class(index, freesurfer_lut):
    """Map a label index to a coarse tissue class."""
    index = int(index)
    if index == 0:
        return "Unknown"
    if not freesurfer_lut:
        return "Unknown"
    if 1000 <= index < 1100 or 2000 <= index < 2100 or index in (3, 42):
        return "GM"
    if index in (8, 47) or index in FS_SUBCORTICAL_GM:
        return "GM"
    if index in FS_WM or 251 <= index <= 255 or 3000 <= index < 5000:
        return "WM"
    if index in FS_CSF:
        return "CSF"

    return "Other"


def hemi_of(name):
    """Hemisphere implied by a structure name, or '?' if it has none."""
    lowered = name.lower()
    if lowered.startswith(("left-", "ctx-lh", "wm-lh")) or "(left)" in lowered:
        return "L"
    if lowered.startswith(("right-", "ctx-rh", "wm-rh")) or "(right)" in lowered:
        return "R"

    return "?"


def neighbourhood(vol, affine, xyz, radius):
    """Return the voxel indices and squared distances within ``radius`` mm."""
    center = (np.linalg.inv(affine) @ np.r_[xyz, 1.0])[:3]
    voxel_size = np.linalg.norm(affine[:3, :3], axis=0)
    radius_vox = np.ceil(radius / voxel_size).astype(int)
    corner = np.round(center).astype(int)

    low = np.maximum(corner - radius_vox, 0)
    high = np.minimum(corner + radius_vox + 1, vol.shape)
    if np.any(low >= high):
        return None, None

    grid = np.meshgrid(
        np.arange(low[0], high[0]),
        np.arange(low[1], high[1]),
        np.arange(low[2], high[2]),
        indexing="ij",
    )
    ijk = np.stack([axis.ravel() for axis in grid], axis=1)
    world = (affine @ np.c_[ijk, np.ones(len(ijk))].T).T[:, :3]
    dist_sq = ((world - xyz) ** 2).sum(1)

    keep = dist_sq <= radius**2
    if not keep.any():
        return None, None

    return ijk[keep], dist_sq[keep]


def gaussian_vote(vol, affine, xyz, sigma, n_sigma):
    """Gaussian-weighted label histogram around the world point ``xyz``.

    Weights sum to one over every sampled voxel, background included, so a
    contact at the edge of the brain keeps probability mass on Unknown instead
    of being renormalised into false confidence.
    """
    ijk, dist_sq = neighbourhood(vol, affine, xyz, n_sigma * sigma)
    if ijk is None:
        return {0: 1.0}

    weights = np.exp(-dist_sq / (2.0 * sigma**2))
    weights /= weights.sum()
    labels = vol[ijk[:, 0], ijk[:, 1], ijk[:, 2]]

    probs = {}
    for label, weight in zip(labels, weights):
        probs[int(label)] = probs.get(int(label), 0.0) + float(weight)

    return probs


def norm_entropy(probs):
    """Shannon entropy of ``probs``, normalised to [0, 1]."""
    values = np.array([value for value in probs.values() if value > 0])
    if len(values) <= 1:
        return 0.0

    return float(-(values * np.log(values)).sum() / np.log(len(values)))


def dist_to_boundary(vol, affine, xyz, ref_label, max_dist):
    """Distance to the nearest voxel labelled differently to ``ref_label``."""
    ijk, dist_sq = neighbourhood(vol, affine, xyz, max_dist)
    if ijk is None:
        return np.nan

    labels = vol[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
    differing = dist_sq[labels != ref_label]
    if not len(differing):
        return np.nan

    return float(np.sqrt(differing.min()))


def nearest_labelled(vol, affine, xyz, max_dist):
    """Nearest non-zero label and its distance, for contacts in background."""
    ijk, dist_sq = neighbourhood(vol, affine, xyz, max_dist)
    if ijk is None:
        return 0, np.nan

    labels = vol[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
    labelled = labels > 0
    if not labelled.any():
        return 0, np.nan

    closest = np.argmin(np.where(labelled, dist_sq, np.inf))

    return int(labels[closest]), float(np.sqrt(dist_sq[closest]))


def confidence(p_top, margin):
    """Coarse confidence rating from the top probability and its margin."""
    if p_top >= 0.80 and margin >= 0.60:
        return "High"
    if p_top >= 0.50 and margin >= 0.20:
        return "Medium"

    return "Low"


def lookup_atlas_labels(
    coords_fcsv,
    atlas_dseg,
    lut_file,
    output_tsv,
    sigma_contact,
    sigma_reg,
    n_sigma,
    max_dist,
):
    """
    Function that assigns atlas labels to every contact of an fcsv file.

    Parameters
    ----------
    coords_fcsv : str
        Path to the contact coordinates, in the same space as the atlas.
    atlas_dseg : str
        Path to the atlas segmentation.
    lut_file : str or None
        Path to a BIDS dseg.tsv lookup table. When empty, the freesurfer
        aseg/aparc table carried by this script is used.
    output_tsv : str
        Path to save one row per contact.
    sigma_contact : float
        Contact extent plus localisation error, in mm.
    sigma_reg : float
        Residual registration error, in mm.
    n_sigma : float
        Neighbourhood radius, in units of the combined sigma.
    max_dist : float
        Cap (mm) on the nearest-structure and boundary-distance searches.

    Returns
    -------
    pandas.DataFrame
    """

    lut, freesurfer_lut = read_lut(lut_file)
    sigma = float(np.hypot(sigma_contact, sigma_reg))

    contacts = pd.read_csv(coords_fcsv, skiprows=FCSV_HEADER_ROWS, header=None)
    coords = contacts[[1, 2, 3]].to_numpy(float)
    names = contacts[11].astype(str)

    seg = nib.load(atlas_dseg)
    vol = np.rint(np.asarray(seg.dataobj)).astype(int)
    affine = seg.affine

    # printed up front so a long run is distinguishable from a stuck one
    print(
        f"[lookup] {atlas_dseg}: {'x'.join(str(n) for n in vol.shape)} grid, "
        f"{len(coords)} contacts, sigma {sigma:.2f} mm",
        flush=True,
    )

    rows = []
    for i, xyz in enumerate(coords):
        ijk = np.rint(np.linalg.inv(affine) @ np.r_[xyz, 1.0])[:3].astype(int)
        in_bounds = all(0 <= ijk[d] < vol.shape[d] for d in range(3))
        hard = int(vol[tuple(ijk)]) if in_bounds else 0
        hard_name = lut.get(hard, f"idx{hard}")

        probs = gaussian_vote(vol, affine, xyz, sigma, n_sigma)
        ranked = sorted(probs.items(), key=lambda item: -item[1])
        top_index, p_top = ranked[0]
        second_index, p_second = ranked[1] if len(ranked) > 1 else (0, 0.0)

        tissue = {}
        for index, prob in probs.items():
            group = tissue_class(index, freesurfer_lut)
            tissue[group] = tissue.get(group, 0.0) + prob

        if hard == 0:
            near_index, near_dist = nearest_labelled(vol, affine, xyz, max_dist)
        else:
            near_index, near_dist = hard, 0.0

        rows.append(
            {
                "name": names.iloc[i],
                "structure": hard_name,
                "tissue": tissue_class(hard, freesurfer_lut),
                "top_structure": lut.get(top_index, f"idx{top_index}"),
                "p_top": p_top,
                "second_structure": lut.get(second_index, f"idx{second_index}"),
                "p_second": p_second,
                "margin": p_top - p_second,
                "entropy": norm_entropy(probs),
                "confidence": confidence(p_top, p_top - p_second),
                "p_GM": tissue.get("GM", 0.0) if freesurfer_lut else np.nan,
                "p_WM": tissue.get("WM", 0.0) if freesurfer_lut else np.nan,
                "p_CSF": tissue.get("CSF", 0.0) if freesurfer_lut else np.nan,
                "n_structures": len(probs),
                "dist_to_boundary_mm": (
                    dist_to_boundary(vol, affine, xyz, hard, max_dist)
                    if in_bounds
                    else np.nan
                ),
                "nearest_structure": (
                    lut.get(near_index, f"idx{near_index}") if near_index else "None"
                ),
                "nearest_dist_mm": near_dist,
                "hemi_match": hemi_of(hard_name) == str(names.iloc[i])[:1].upper(),
            }
        )

    labels = pd.DataFrame(rows)
    report_concordance(labels, atlas_dseg)
    labels.to_csv(output_tsv, sep="\t", index=False, float_format="%.3f")

    return labels


def report_concordance(labels, atlas_dseg):
    """Print how many contacts got a label, and whether hemispheres agree.

    A registration that silently failed, or a segmentation in the wrong frame,
    shows up here as a collapse in either number.
    """
    total = len(labels)
    if not total:
        return

    labelled = int((labels["structure"] != "Unknown").sum())
    hemi_matched = int(labels["hemi_match"].sum())
    print(
        f"[concordance] {atlas_dseg}: labelled {labelled}/{total} "
        f"({labelled / total:.0%}), hemisphere match {hemi_matched}/{total} "
        f"({hemi_matched / total:.0%})"
    )
    if labelled / total < 0.80 or hemi_matched / total < 0.80:
        print(
            "  !! low concordance -- check that the segmentation and the "
            "contacts share a coordinate frame"
        )


if __name__ == "__main__":
    lookup_atlas_labels(
        coords_fcsv=snakemake.input.coords,
        atlas_dseg=snakemake.input.atlas_dseg,
        lut_file=first_or_none(snakemake.input.lut),
        output_tsv=snakemake.output.labels,
        sigma_contact=snakemake.params.sigma_contact,
        sigma_reg=snakemake.params.sigma_reg,
        n_sigma=snakemake.params.n_sigma,
        max_dist=snakemake.params.max_search,
    )
