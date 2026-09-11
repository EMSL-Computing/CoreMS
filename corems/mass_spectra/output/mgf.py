"""SIRIUS-compatible MGF export for LC-MS mass features."""

import warnings
from pathlib import Path

import numpy as np


def sirius_charge(polarity: str) -> str:
    """Map LCMS polarity to a SIRIUS CHARGE token.

    Parameters
    ----------
    polarity : str
        ``'positive'`` or ``'negative'``.

    Returns
    -------
    str
        ``'1+'`` or ``'1-'``.
    """
    value = (polarity or "").strip().lower()
    if value == "positive":
        return "1+"
    if value == "negative":
        return "1-"
    raise ValueError(
        "Polarity not set for dataset, must be a either 'positive' or 'negative'"
    )


def format_ion_block(headers: dict, mz, intensity) -> str:
    """Serialize one MGF spectrum as ``BEGIN IONS`` ... ``END IONS``.

    Parameters
    ----------
    headers : dict
        MGF ``NAME=VALUE`` pairs, written in insertion order.
    mz, intensity : array_like
        Peak lists.

    Returns
    -------
    str
        One ion block including a trailing newline.
    """
    lines = ["BEGIN IONS"]
    for key, value in headers.items():
        if value is None:
            continue
        lines.append(f"{key}={value}")
    mz = np.asarray(mz, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    for mass, abund in zip(mz, intensity):
        lines.append(f"{mass} {abund}")
    lines.append("END IONS")
    return "\n".join(lines) + "\n"


def ensure_mgf_path(out_file_path) -> Path:
    """Return a Path, appending ``.mgf`` when no suffix is given."""
    path = Path(out_file_path)
    if path.suffix == "":
        path = path.with_suffix(".mgf")
    return path


def peak_mz_abundance(spectrum):
    """Return (mz, abundance) from processed peaks, else raw arrays.

    Parameters
    ----------
    spectrum : MassSpectrum or None
        Spectrum to read.

    Returns
    -------
    tuple of numpy.ndarray
        ``(mz, abundance)``. Empty arrays if ``spectrum`` is None or has no peaks.
    """
    if spectrum is None:
        return np.array([]), np.array([])
    mspeaks = getattr(spectrum, "mspeaks", None)
    if mspeaks:
        return np.asarray(spectrum.mz_exp, dtype=float), np.asarray(
            spectrum.abundance, dtype=float
        )
    mz = getattr(spectrum, "_mz_exp", None)
    ab = getattr(spectrum, "_abundance", None)
    if mz is None or ab is None:
        return np.array([]), np.array([])
    return np.asarray(mz, dtype=float), np.asarray(ab, dtype=float)


def _ms1_precursor_peak(feature):
    """Return the feature precursor as a one-peak MS1 list.

    Raises
    ------
    ValueError
        If the feature has no precursor m/z.
    """
    mz = getattr(feature, "mz", None)
    if mz is None:
        raise ValueError(
            "Mass feature is missing precursor m/z; cannot export MGF"
        )
    intensity = getattr(feature, "intensity", None)
    if intensity is None:
        intensity = 0.0
    return np.array([float(mz)]), np.array([float(intensity)])


def _scan_of_spectrum(feature, spec):
    scan = getattr(spec, "scan_number", None)
    if scan is not None:
        return scan
    for key, value in getattr(feature, "ms2_mass_spectra", {}).items():
        if value is spec:
            return key
    return None


def _ms2_spectra(feature, ms2_mode: str):
    if ms2_mode not in ("best", "all"):
        raise ValueError(f"ms2_mode must be 'best' or 'all', got {ms2_mode!r}")
    if ms2_mode == "best":
        spec = feature.best_ms2
        if spec is None:
            return []
        return [(_scan_of_spectrum(feature, spec), spec)]
    spectra = []
    for scan, spec in feature.ms2_mass_spectra.items():
        mz, _ = peak_mz_abundance(spec)
        if mz.size:
            spectra.append((scan, spec))
    return spectra


def _ordered_headers(headers: dict) -> dict:
    order = [
        "FEATURE_ID",
        "PEPMASS",
        "CHARGE",
        "MSLEVEL",
        "TITLE",
        "RTINSECONDS",
        "SCANS",
    ]
    out = {}
    for key in order:
        if key in headers:
            out[key] = headers[key]
    for key, value in headers.items():
        if key not in out:
            out[key] = value
    return out


def _shared_headers(feature, feature_id, polarity, sample_name=None, title=None):
    if title is None:
        if sample_name:
            title = f"{sample_name} feature {feature_id}"
        else:
            title = f"feature {feature_id}"
    rt = getattr(feature, "retention_time", None)
    rtinseconds = None if rt is None else float(rt) * 60.0
    return {
        "FEATURE_ID": feature_id,
        "PEPMASS": feature.mz,
        "CHARGE": sirius_charge(polarity),
        "TITLE": title,
        "RTINSECONDS": rtinseconds,
    }


def feature_mgf_text(
    feature,
    feature_id,
    polarity,
    *,
    sample_name=None,
    ms2_mode="best",
    title=None,
):
    """Return MGF text for one feature, or None if MS2 is missing.

    Parameters
    ----------
    feature : LCMSMassFeature
        Mass feature with associated MS2 spectra. MS1 is the precursor
        ``mz`` / ``intensity`` on the feature, not a scanned spectrum.
    feature_id : hashable
        Value written as ``FEATURE_ID`` (mass-feature id or consensus cluster).
    polarity : str
        ``'positive'`` or ``'negative'``.
    sample_name : str, optional
        Used in the default TITLE.
    ms2_mode : {'best', 'all'}, optional
        Which associated MS2 spectra to write.
    title : str, optional
        Override TITLE.

    Returns
    -------
    str or None
        One or more ion blocks, or None if the feature is incomplete.
    """
    mz1, ab1 = _ms1_precursor_peak(feature)
    ms2_list = _ms2_spectra(feature, ms2_mode)
    if not ms2_list:
        return None

    shared = _shared_headers(
        feature, feature_id, polarity, sample_name=sample_name, title=title
    )
    blocks = []
    ms1_headers = dict(shared)
    ms1_headers["MSLEVEL"] = 1
    ms1_headers["SCANS"] = getattr(feature, "apex_scan", None)
    blocks.append(format_ion_block(_ordered_headers(ms1_headers), mz1, ab1))
    for scan, spec in ms2_list:
        mz2, ab2 = peak_mz_abundance(spec)
        if mz2.size == 0:
            continue
        ms2_headers = dict(shared)
        ms2_headers["MSLEVEL"] = 2
        ms2_headers["SCANS"] = (
            scan if scan is not None else getattr(spec, "scan_number", None)
        )
        blocks.append(format_ion_block(_ordered_headers(ms2_headers), mz2, ab2))
    if len(blocks) < 2:
        return None
    return "".join(blocks)


def write_feature_records_to_mgf(
    records,
    polarity,
    out_file_path,
    *,
    ms2_mode="best",
    overwrite=False,
) -> Path:
    """Write selected feature records to an MGF file.

    Parameters
    ----------
    records : iterable of tuple
        ``(feature_id, feature, sample_name)`` triples. ``sample_name`` may be None.
    polarity : str
        ``'positive'`` or ``'negative'``.
    out_file_path : str or Path
        Output path. ``.mgf`` is appended when no suffix is given.
    ms2_mode : {'best', 'all'}, optional
        Which associated MS2 spectra to write.
    overwrite : bool, optional
        Replace an existing file. Default is False.

    Returns
    -------
    pathlib.Path
        Path of the written MGF file.

    Raises
    ------
    FileExistsError
        If the path exists and ``overwrite`` is False.
    ValueError
        If no complete MS1/MS2 pairs remain after skipping.
    """
    path = ensure_mgf_path(out_file_path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"MGF file already exists: {path}")

    chunks = []
    n_in = 0
    n_skip = 0
    for feature_id, feature, sample_name in records:
        n_in += 1
        text = feature_mgf_text(
            feature,
            feature_id,
            polarity,
            sample_name=sample_name,
            ms2_mode=ms2_mode,
        )
        if text is None:
            n_skip += 1
            continue
        chunks.append(text)

    if not chunks:
        raise ValueError(
            "No complete MS1/MS2 feature pairs to export "
            f"(skipped {n_skip} of {n_in} selected features)"
        )
    if n_skip:
        warnings.warn(
            f"Skipped {n_skip} of {n_in} features missing a usable MS2 spectrum",
            UserWarning,
        )
    path.write_text("".join(chunks), encoding="utf-8")
    return path


def iter_collection_mgf_records(collection, cluster_ids=None):
    """Build ``(cluster_id, mass_feature, sample_name)`` from consensus representatives.

    Parameters
    ----------
    collection : LCMSCollection
        Collection with consensus clustering and loaded representatives.
    cluster_ids : iterable, optional
        Cluster ids to export. Default is all representative clusters.

    Returns
    -------
    list of tuple
        ``(cluster_id, LCMSMassFeature, sample_name)``.

    Raises
    ------
    ValueError
        If representatives are missing, requested clusters are unknown, or
        representative objects are not loaded on the samples.
    """
    if not hasattr(collection, "get_representative_mass_features_for_all_clusters"):
        raise ValueError("Collection cannot provide consensus representatives")
    reps = collection.get_representative_mass_features_for_all_clusters()
    if reps is None or len(reps) == 0:
        raise ValueError(
            "No consensus representatives found. Run process_consensus_features() "
            "with load_representatives=True, add_ms2=True first."
        )
    if cluster_ids is None:
        selected = reps
    else:
        requested = list(cluster_ids)
        selected = reps[reps["cluster"].isin(requested)]
        found = set(selected["cluster"].tolist())
        missing = [c for c in requested if c not in found]
        if missing:
            raise ValueError(f"Consensus cluster id(s) not found: {missing}")

    records = []
    lcms_map = getattr(collection, "_lcms", {})
    for _, row in selected.iterrows():
        cluster = row["cluster"]
        sample_name = None
        if "sample_name" in selected.columns:
            val = row["sample_name"]
            if not (isinstance(val, (float, np.floating)) and np.isnan(val)):
                sample_name = val
        if sample_name is None:
            if "sample_id" not in selected.columns:
                raise ValueError(
                    "Cannot resolve sample_name from consensus representatives"
                )
            sample_name = collection.samples[int(row["sample_id"])]
        mf_id = row["mf_id"]
        sample = lcms_map.get(sample_name)
        mass_features = getattr(sample, "mass_features", {}) if sample is not None else {}
        if mf_id not in mass_features:
            try:
                mf_id = int(mf_id)
            except (TypeError, ValueError):
                mf_id = row["mf_id"]
        if sample is None or mf_id not in mass_features:
            raise ValueError(
                "Representative mass features have not been loaded into samples. "
                "Call process_consensus_features() with load_representatives=True, "
                "add_ms2=True before exporting MGF."
            )
        records.append((cluster, mass_features[mf_id], sample_name))
    return records
