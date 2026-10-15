"""Greedy mass-feature rollup keeps parents without an all-pairs matrix."""

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import KDTree

from corems.mass_spectra.calc.lc_calc import PHCalculations


def _roll(df, dims, tol, relative, sort_by="persistence", **kwargs):
    return PHCalculations().roll_up_dataframe(
        df=df,
        sort_by=sort_by,
        dims=dims,
        tol=tol,
        relative=relative,
        **kwargs,
    )


def _block_pair_matrix(monkeypatch):
    def _boom(*_args, **_kwargs):
        raise AssertionError("sparse_distance_matrix was called")

    monkeypatch.setattr(
        "corems.mass_spectra.calc.lc_calc.KDTree.sparse_distance_matrix",
        _boom,
    )


def _claims(parent_row, child_row, tol, relative):
    """Float32 predicate. parent_row and child_row are float32 vectors."""
    for d, rel in enumerate(relative):
        diff = np.float32(np.abs(np.float32(child_row[d]) - np.float32(parent_row[d])))
        if rel:
            limit = np.float32(np.float32(tol[d]) * np.float32(parent_row[d]))
        else:
            limit = np.float32(tol[d])
        if diff > limit:
            return False
    return True


def _greedy_positions(values, tol, relative):
    """Independent greedy scan. values are float32 columns, row 0 is highest sort."""
    n = len(values[0])
    kept = []
    for i in range(n):
        child = [values[d][i] for d in range(len(values))]
        if any(
            _claims([values[d][p] for d in range(len(values))], child, tol, relative)
            for p in kept
        ):
            continue
        kept.append(i)
    return kept


def _legacy_rollup(df, dims, tol, relative, sort_by="persistence"):
    """Today's KDTree, upper triangle, and while-loop. Small frames only."""
    og_columns = df.columns.copy()
    og_index = df.index.name if df.index.name is not None else "index"
    work = df.reset_index(drop=False)
    work = work.sort_values(by=sort_by, ascending=False).reset_index(drop=True)
    all_values = [work[dim].to_numpy() for dim in dims]
    distances = None
    for i, col in enumerate(all_values):
        values = np.asarray(col, dtype=np.float32)
        tree = KDTree(values.reshape(-1, 1))
        max_tol = tol[i] * float(values.max()) if relative[i] else tol[i]
        sdm = tree.sparse_distance_matrix(tree, max_tol, output_type="coo_matrix")
        sdm = sparse.triu(sdm, k=1)
        if relative[i]:
            valid = sdm.data <= tol[i] * values[sdm.row]
            sdm = sparse.coo_matrix(
                (
                    np.ones(int(valid.sum()), dtype=np.uint8),
                    (sdm.row[valid], sdm.col[valid]),
                ),
                shape=(len(values), len(values)),
            )
        else:
            sdm.data = np.ones(len(sdm.data), dtype=np.uint8)
        distances = sdm if distances is None else distances.multiply(sdm)
    distances = distances.tocoo()
    pairs = (
        np.stack((distances.row, distances.col), axis=1)
        if distances.nnz
        else np.empty((0, 2), dtype=int)
    )
    pairs_df = pd.DataFrame(pairs, columns=["parent", "child"]).set_index("parent")
    to_drop = []
    while not pairs_df.empty:
        root_parents = np.setdiff1d(
            np.unique(pairs_df.index.values), np.unique(pairs_df.child.values)
        )
        children = pairs_df.loc[root_parents, "child"].unique()
        to_drop.extend(children)
        pairs_df = pairs_df.drop(index=children, errors="ignore")
        pairs_df = pairs_df.reset_index().set_index("child")
        pairs_df = pairs_df.drop(index=children)
        pairs_df = pairs_df.reset_index().set_index("parent")
    out = work.drop(index=np.array(to_drop, dtype=int))
    return out.set_index(og_index).sort_index()[og_columns]


def _legacy_frame():
    rows = []
    persistence = 10_000
    for group in range(40):
        mz = 200.0 + group * 0.1
        rt = 1.0 + group * 0.01
        rows.append((persistence, mz, rt, f"p{group}"))
        persistence -= 1
        for child in range(4):
            rows.append(
                (persistence, mz * (1 + 2e-6), rt + 0.10, f"c{group}_{child}")
            )
            persistence -= 1
        rows.append((persistence, mz, rt + 0.50, f"f{group}"))
        persistence -= 1
    df = pd.DataFrame(rows, columns=["persistence", "mz", "scan_time", "id"])
    df.index = np.arange(100, 100 + len(df))
    return df


def test_chain_keeps_distant_grandchild():
    df = pd.DataFrame(
        {
            "persistence": [3.0, 2.0, 1.0],
            "mz": [200.0000, 200.0004, 200.0008],
            "scan_time": [1.00, 1.20, 1.45],
            "id": ["A", "B", "C"],
        },
        index=[10, 20, 30],
    )
    out = _roll(df, ["mz", "scan_time"], [5e-6, 0.3], [True, False])
    assert list(out["id"]) == ["A", "C"]
    assert list(out.index) == [10, 30]


def test_parent_referenced_mz_tolerance():
    tol = np.float32(1e-3)
    parent = np.float32(400.0)
    parent_limit = np.float32(tol * parent)
    child_at_limit = np.float32(parent - parent_limit)
    child_limit = np.float32(tol * child_at_limit)
    diff = np.float32((parent_limit + child_limit) / np.float32(2))
    near = np.float32(parent - diff)
    far_diff = np.float32(parent_limit + np.float32(0.05))
    far = np.float32(parent - far_diff)
    assert diff <= parent_limit
    assert diff > np.float32(tol * near)
    assert far_diff > parent_limit
    # The dropped middle is inside tolerance of the far row. Only the parent may drop it.
    assert np.float32(np.abs(far - near)) <= np.float32(tol * near)
    df = pd.DataFrame(
        {
            "persistence": [3.0, 2.0, 1.0],
            "mz": [float(parent), float(near), float(far)],
            "scan_time": [1.0, 1.0, 1.0],
            "id": ["P", "N", "F"],
        }
    )
    out = _roll(df, ["mz", "scan_time"], [float(tol), 0.3], [True, False])
    assert list(out["id"]) == ["P", "F"]


def test_legacy_oracle_matches_and_skips_pair_matrix(monkeypatch):
    df = _legacy_frame()
    expected = _legacy_rollup(df, ["mz", "scan_time"], [5e-6, 0.3], [True, False])
    assert set(expected["id"]) == {f"p{i}" for i in range(40)} | {f"f{i}" for i in range(40)}
    _block_pair_matrix(monkeypatch)
    out = _roll(df, ["mz", "scan_time"], [5e-6, 0.3], [True, False])
    pd.testing.assert_frame_equal(out, expected)


def test_dense_retention_time_matches_greedy_scan(monkeypatch):
    n = 3000
    n_mz = 30
    df = pd.DataFrame(
        {
            "persistence": np.arange(n, 0, -1, dtype=float),
            "mz": 200.0 + (np.arange(n) % n_mz) * 0.01,
            "scan_time": 1.0 + (np.arange(n) % 100) * 0.002,
        }
    )
    _block_pair_matrix(monkeypatch)
    out = _roll(df, ["mz", "scan_time"], [5e-6, 0.3], [True, False])
    ordered = df.sort_values("persistence", ascending=False)
    values = [ordered[col].to_numpy(dtype=np.float32) for col in ("mz", "scan_time")]
    kept = _greedy_positions(values, [5e-6, 0.3], [True, False])
    assert len(out) == len(kept) == n_mz
    assert set(out.index) == set(ordered.index[kept])


def test_repeated_mz_invariants_and_allocation(monkeypatch):
    import tracemalloc

    n = 100_000
    n_mz = 200
    rng = np.random.default_rng(0)
    mz_id = np.arange(n) % n_mz
    df = pd.DataFrame(
        {
            "persistence": np.arange(n, 0, -1, dtype=np.int64),
            "mz": 200.0 + mz_id * 0.01,
            "scan_time": rng.uniform(0.0, 10.0, size=n),
            "mz_id": mz_id,
        }
    )
    _block_pair_matrix(monkeypatch)
    tracemalloc.start()
    out = _roll(df, ["mz", "scan_time"], [5e-6, 0.3], [True, False])
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < 256 * 1024 * 1024

    # 0.01 Da is wider than 5 ppm at these masses, so a parent shares mz_id.
    max_mz = float(df["mz"].max())
    assert 0.01 > 5e-6 * max_mz
    kept_ids = set(out.index)
    for _, group in df.groupby("mz_id"):
        group = group.sort_values("persistence", ascending=False)
        values = [
            group[col].to_numpy(dtype=np.float32) for col in ("mz", "scan_time")
        ]
        kept_pos = [
            i for i, label in enumerate(group.index) if label in kept_ids
        ]
        for i in kept_pos:
            for j in kept_pos:
                if group.iloc[j]["persistence"] <= group.iloc[i]["persistence"]:
                    continue
                parent = [values[d][j] for d in range(2)]
                child = [values[d][i] for d in range(2)]
                assert not _claims(parent, child, [5e-6, 0.3], [True, False])
        for i, label in enumerate(group.index):
            if label in kept_ids:
                continue
            child = [values[d][i] for d in range(2)]
            assert any(
                _claims(
                    [values[d][j] for d in range(2)],
                    child,
                    [5e-6, 0.3],
                    [True, False],
                )
                for j in kept_pos
                if group.iloc[j]["persistence"] > group.iloc[i]["persistence"]
            )


def test_threshold_argument_does_not_change_rows(monkeypatch):
    df = _legacy_frame()
    _block_pair_matrix(monkeypatch)
    plain = _roll(df, ["mz", "scan_time"], [5e-6, 0.3], [True, False])
    forced = _roll(
        df,
        ["mz", "scan_time"],
        [5e-6, 0.3],
        [True, False],
        memory_opt_threshold=10,
    )
    pd.testing.assert_frame_equal(plain, forced)


def test_index_columns_empty_and_one_row(monkeypatch):
    _block_pair_matrix(monkeypatch)
    named = pd.DataFrame(
        {
            "persistence": [2.0, 1.0],
            "mz": [300.0, 300.0002],
            "scan_time": [2.0, 2.05],
            "area": [5, 9],
        },
        index=pd.Index([5, 8], name="feat_id"),
    )
    out = _roll(named, ["mz", "scan_time"], [5e-6, 0.3], [True, False])
    assert out.index.name == "feat_id"
    assert list(out.index) == [5]
    assert list(out.columns) == ["persistence", "mz", "scan_time", "area"]
    assert int(out.iloc[0]["area"]) == 5

    one = named.iloc[[0]].copy()
    one_out = _roll(one, ["mz", "scan_time"], [5e-6, 0.3], [True, False])
    assert list(one_out.index) == [5]
    assert float(one_out.iloc[0]["mz"]) == 300.0

    empty = named.iloc[0:0].copy()
    empty_out = _roll(empty, ["mz", "scan_time"], [5e-6, 0.3], [True, False])
    assert len(empty_out) == 0
    assert list(empty_out.columns) == list(named.columns)


def test_one_absolute_dimension(monkeypatch):
    _block_pair_matrix(monkeypatch)
    df = pd.DataFrame(
        {
            "persistence": [5.0, 4.0, 3.0, 2.0, 1.0],
            "x": [0.0, 0.2, 0.4, 2.0, 2.2],
        }
    )
    out = _roll(df, ["x"], [0.5], [False])
    ordered = df.sort_values("persistence", ascending=False)
    values = [ordered["x"].to_numpy(dtype=np.float32)]
    kept = _greedy_positions(values, [0.5], [False])
    assert list(out.index) == list(ordered.index[kept])


def test_negative_tolerance_keeps_every_row(monkeypatch):
    _block_pair_matrix(monkeypatch)
    df = pd.DataFrame({"persistence": [2.0, 1.0], "x": [0.0, 0.1]})
    out = _roll(df, ["x"], [-1.0], [False])
    assert len(out) == 2


def test_zero_tolerance_is_float32_equality(monkeypatch):
    _block_pair_matrix(monkeypatch)
    df = pd.DataFrame(
        {
            "persistence": [3.0, 2.0, 1.0],
            "x": [np.float32(1.25), np.float32(1.25), np.float32(1.5)],
        }
    )
    out = _roll(df, ["x"], [0.0], [False])
    assert list(out.index) == [0, 2]


def test_nonfinite_row_is_kept_and_claims_nobody(monkeypatch):
    _block_pair_matrix(monkeypatch)
    df = pd.DataFrame(
        {
            "persistence": [3.0, 2.0, 1.0],
            "x": [np.nan, 1.0, 1.0],
        }
    )
    out = _roll(df, ["x"], [0.5], [False])
    assert list(out.index) == [0, 1]
