"""Parquet sidecar cache for parsed spectral library DataFrames."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

# Bump when standardized cache columns or precursor policy change.
COREMS_LIB_CACHE_V = 1
CACHE_SUFFIX = ".corems-lib.parquet"
_META_KEY = b"corems_lib_cache"


def default_cache_path(source_path) -> Path:
    """Return default sidecar path next to the source library file."""
    return Path(source_path).with_name(Path(source_path).name + CACHE_SUFFIX)


def _source_fingerprint(source_path) -> dict:
    path = Path(source_path)
    stat = path.stat()
    return {
        "schema_version": COREMS_LIB_CACHE_V,
        "source_name": path.name,
        "source_size": stat.st_size,
        "source_mtime_ns": stat.st_mtime_ns,
    }


def is_valid(cache_path, source_path, extra_meta: dict | None = None) -> bool:
    """Return True if cache exists and matches source fingerprint + schema."""
    cache_path = Path(cache_path)
    source_path = Path(source_path)
    if not cache_path.is_file() or not source_path.is_file():
        return False
    try:
        import pyarrow.parquet as pq

        meta = pq.read_schema(cache_path).metadata or {}
        raw = meta.get(_META_KEY)
        if raw is None:
            return False
        stored = json.loads(raw.decode("utf-8"))
    except Exception:
        return False

    expected = _source_fingerprint(source_path)
    if extra_meta:
        expected.update(extra_meta)
    for key, value in expected.items():
        if stored.get(key) != value:
            return False
    return True


def pack_peaks_for_parquet(df: pd.DataFrame) -> pd.DataFrame:
    """Convert a peaks column of (N, 2) arrays into mz/intensity list columns."""
    out = df.copy()
    if "peaks" not in out.columns:
        return out
    mz_col = []
    intensity_col = []
    for peaks in out["peaks"]:
        if peaks is None or (isinstance(peaks, float) and np.isnan(peaks)):
            mz_col.append([])
            intensity_col.append([])
            continue
        arr = np.asarray(peaks, dtype=np.float64)
        if arr.size == 0:
            mz_col.append([])
            intensity_col.append([])
        elif arr.ndim == 2 and arr.shape[1] == 2:
            mz_col.append(arr[:, 0].tolist())
            intensity_col.append(arr[:, 1].tolist())
        else:
            raise ValueError("peaks entries must be array-like of shape (N, 2)")
    out["mz"] = mz_col
    out["intensity"] = intensity_col
    out = out.drop(columns=["peaks"])
    return out


def unpack_peaks_from_parquet(df: pd.DataFrame) -> pd.DataFrame:
    """Rebuild peaks (N, 2) arrays from mz/intensity list columns."""
    out = df.copy()
    if "mz" not in out.columns or "intensity" not in out.columns:
        return out
    peaks = []
    for mz, intensity in zip(out["mz"], out["intensity"]):
        mz_arr = np.asarray(mz if mz is not None else [], dtype=np.float64)
        i_arr = np.asarray(intensity if intensity is not None else [], dtype=np.float64)
        if mz_arr.size == 0:
            peaks.append(np.zeros((0, 2), dtype=np.float64))
        else:
            peaks.append(np.column_stack([mz_arr, i_arr]))
    out["peaks"] = peaks
    out = out.drop(columns=["mz", "intensity"])
    return out


def read(cache_path) -> pd.DataFrame:
    """Load a cached library DataFrame and restore peaks arrays."""
    df = pd.read_parquet(cache_path)
    return unpack_peaks_from_parquet(df)


def write(cache_path, df: pd.DataFrame, source_path, extra_meta: dict | None = None) -> None:
    """Write standardized library DataFrame as a parquet sidecar."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    cache_path = Path(cache_path)
    packed = pack_peaks_for_parquet(df)
    table = pa.Table.from_pandas(packed, preserve_index=False)
    meta = _source_fingerprint(source_path)
    if extra_meta:
        meta.update(extra_meta)
    existing = table.schema.metadata or {}
    existing = dict(existing)
    existing[_META_KEY] = json.dumps(meta).encode("utf-8")
    table = table.replace_schema_metadata(existing)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, cache_path)
