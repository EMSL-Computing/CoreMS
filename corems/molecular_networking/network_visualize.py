"""Visualization mixin for molecular networking objects.

This module provides :class:`NetworkVisualizeMixin`, which adds interactive
HTML network plotting via ``networkx`` and ``ipysigma`` for classes exposing:

- ``similarity_matrices``: ``dict[str, SimilarityMatrix]``
- ``_all_query_ids``: ``list[str]``
- ``_threshold_for(metric) -> float``
"""

from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
from typing import Any, Sequence

import pandas as pd


class NetworkVisualizeMixin:
    """Mixin that adds HTML network plotting for similarity matrices."""

    _DEFAULT_COLOR_MAP = {
        "query": "#e74c3c",
        "library": "#3498db",
        "unknown": "#95a5a6",
    }

    _DEFAULT_LIBRARY_LABEL_FIELDS = ("compound_name", "name", "spectra_id", "id")

    def _library_record_for_node(self, node_id: str) -> dict[str, Any] | None:
        """Return FE library record dict for a library node ID, if available."""
        fe_lib = getattr(self, "fe_lib", None)
        if fe_lib is None:
            return None

        try:
            lib_idx = int(node_id)
        except (TypeError, ValueError):
            return None

        try:
            record = fe_lib[lib_idx]
        except Exception:
            return None

        return record if isinstance(record, dict) else None

    @staticmethod
    def _as_node_attr_value(value: Any) -> Any:
        """Convert values to JSON-friendly node attributes for Sigma hover panels."""
        if value is None or isinstance(value, (str, int, float, bool)):
            return value
        return str(value)

    def _resolve_library_label(
        self,
        node_id: str,
        record: dict[str, Any] | None,
        library_label_field: str | Sequence[str] | None,
    ) -> str:
        """Resolve display label for a library node from FE record fields."""
        fallback = node_id[:20] + "..." if len(node_id) > 20 else node_id
        if not record:
            return fallback

        if library_label_field is None:
            candidate_fields: Sequence[str] = self._DEFAULT_LIBRARY_LABEL_FIELDS
        elif isinstance(library_label_field, str):
            candidate_fields = (library_label_field,)
        else:
            candidate_fields = library_label_field

        for field in candidate_fields:
            value = record.get(field)
            if value is None:
                continue
            value_str = str(value).strip()
            if value_str:
                return value_str

        return fallback

    def _edges_to_dataframe(
        self,
        metric: str,
        score_threshold: float,
        include_queries_only: bool,
        exclude_self: bool,
    ) -> pd.DataFrame:
        """Return filtered edge list as a DataFrame with id1/id2/score."""
        mat = self.similarity_matrices[metric]
        all_edges = mat.get_pairs_above_threshold(0.0)

        if not all_edges:
            return pd.DataFrame(columns=["id1", "id2", "score"])

        edges = pd.DataFrame(all_edges, columns=["id1", "id2", "score"])
        edges = edges[edges["score"] >= score_threshold]

        if exclude_self:
            edges = edges[edges["id1"] != edges["id2"]]

        if include_queries_only and self._all_query_ids:
            query_id_set = set(self._all_query_ids)

            # Keep all query edges plus library-library edges where both
            # library nodes are connected to at least one query node.
            q1 = edges["id1"].isin(query_id_set)
            q2 = edges["id2"].isin(query_id_set)
            ql_or_qq = q1 | q2

            query_edges = edges[ql_or_qq]
            if query_edges.empty:
                return pd.DataFrame(columns=["id1", "id2", "score"])

            lib_from_id1 = query_edges.loc[~query_edges["id1"].isin(query_id_set), "id1"]
            lib_from_id2 = query_edges.loc[~query_edges["id2"].isin(query_id_set), "id2"]
            query_connected_libs = set(pd.concat([lib_from_id1, lib_from_id2]).tolist())

            ll_mask = (
                ~q1
                & ~q2
                & edges["id1"].isin(query_connected_libs)
                & edges["id2"].isin(query_connected_libs)
            )

            edges = pd.concat([query_edges, edges[ll_mask]], ignore_index=True).drop_duplicates()

        return edges

    @staticmethod
    def _sparsify_top_k_edges(G, weight_attr: str, top_k: int):
        """Return an undirected top-k edge sparsified copy of *G*."""
        if top_k <= 0 or G.number_of_nodes() == 0:
            return G.copy()

        H = G.to_undirected().copy()
        selected_edges: set[tuple[str, str]] = set()

        for node in H.nodes:
            neighbors = []
            for nbr in H.neighbors(node):
                edge_data = H.get_edge_data(node, nbr, default={})
                weight = float(edge_data.get(weight_attr, 1.0))
                neighbors.append((nbr, weight))

            neighbors.sort(key=lambda t: (-t[1], str(t[0])))
            for nbr, _ in neighbors[:top_k]:
                edge_key = tuple(sorted((str(node), str(nbr))))
                selected_edges.add(edge_key)

        S = H.__class__()
        S.add_nodes_from(H.nodes(data=True))
        for u, v, data in H.edges(data=True):
            edge_key = tuple(sorted((str(u), str(v))))
            if edge_key in selected_edges:
                S.add_edge(u, v, **data)

        return S

    def _detect_weighted_communities(
        self,
        G,
        *,
        weight_attr: str,
        super_threshold: int,
        sparsify_top_k: int,
        recursive_split: bool,
    ) -> list[set]:
        """Detect weighted communities with optional recursive splitting."""

        def partition_graph(graph, depth: int = 0) -> list[set]:
            if graph.number_of_nodes() <= 1:
                return [set(graph.nodes)]

            communities = list(
                __import__("networkx").community.greedy_modularity_communities(
                    graph,
                    weight=weight_attr,
                )
            )
            if len(communities) <= 1:
                return [set(graph.nodes)]

            sorted_communities = sorted(communities, key=lambda c: (-len(c), min(str(n) for n in c)))
            out: list[set] = []

            for community in sorted_communities:
                community_set = set(community)
                if not recursive_split or len(community_set) <= super_threshold:
                    out.append(community_set)
                    continue

                sub = graph.subgraph(community_set).copy()
                if sparsify_top_k > 0:
                    sub = self._sparsify_top_k_edges(sub, weight_attr=weight_attr, top_k=sparsify_top_k)

                if sub.number_of_edges() == 0:
                    out.append(community_set)
                    continue

                sub_parts = partition_graph(sub, depth + 1)
                if len(sub_parts) <= 1:
                    out.append(community_set)
                    continue

                for part in sub_parts:
                    if part:
                        out.append(set(part))

            return out

        raw = partition_graph(G.to_undirected())
        normalized = [set(c) for c in raw if c]
        if not normalized:
            return [set(G.nodes)]

        return sorted(normalized, key=lambda c: (-len(c), min(str(n) for n in c)))

    @staticmethod
    def _compute_hierarchical_layout(G, communities: list[set], *, weight_attr: str, seed: int):
        """Compute two-level (community + local) weighted spring layout."""
        nx = __import__("networkx")

        if G.number_of_nodes() == 0:
            return {}

        if not communities:
            pos = nx.spring_layout(
                G,
                seed=seed,
                weight=weight_attr,
                k=2.0 / max(G.number_of_nodes() ** 0.5, 1),
            )
            return {n: {"x": float(x), "y": float(y)} for n, (x, y) in pos.items()}

        node_to_comm = {}
        for comm_idx, nodes in enumerate(communities):
            for node in nodes:
                node_to_comm[node] = comm_idx

        meta = nx.Graph()
        for comm_idx, nodes in enumerate(communities):
            meta.add_node(comm_idx, size=len(nodes))

        inter_weights: dict[tuple[int, int], float] = defaultdict(float)
        for u, v, data in G.edges(data=True):
            cu = node_to_comm.get(u)
            cv = node_to_comm.get(v)
            if cu is None or cv is None or cu == cv:
                continue
            edge = (cu, cv) if cu < cv else (cv, cu)
            inter_weights[edge] += float(data.get(weight_attr, 1.0))

        for (cu, cv), w in inter_weights.items():
            meta.add_edge(cu, cv, weight=w)

        # Build local layouts first to estimate community footprint radius.
        local_layouts: dict[int, dict] = {}
        community_radius: dict[int, float] = {}
        for comm_idx, nodes in enumerate(communities):
            sub = G.subgraph(nodes)
            if sub.number_of_nodes() == 1:
                only = next(iter(sub.nodes))
                local_pos = {only: (0.0, 0.0)}
            else:
                local_pos = nx.spring_layout(
                    sub,
                    seed=seed,
                    weight=weight_attr,
                    k=1.5 / max(sub.number_of_nodes() ** 0.5, 1),
                )

            max_radius = 0.0
            for x, y in local_pos.values():
                max_radius = max(max_radius, float((x * x + y * y) ** 0.5))
            if max_radius <= 0:
                max_radius = 1.0

            radius = 0.6 + 0.18 * (len(nodes) ** 0.5)
            local_layouts[comm_idx] = local_pos
            community_radius[comm_idx] = radius

        if meta.number_of_nodes() == 1:
            meta_pos = {0: (0.0, 0.0)}
        else:
            # Larger k spreads community centers to reduce overlap.
            meta_pos = nx.spring_layout(
                meta,
                seed=seed,
                weight="weight",
                k=6.0 / max(meta.number_of_nodes() ** 0.5, 1),
            )

        # Apply pairwise separation pass so center distances exceed
        # footprint-based minimum spacing.
        meta_xy = {idx: [float(x), float(y)] for idx, (x, y) in meta_pos.items()}
        comm_ids = list(meta_xy.keys())
        for _ in range(8):
            moved = False
            for i in range(len(comm_ids)):
                ci = comm_ids[i]
                xi, yi = meta_xy[ci]
                for j in range(i + 1, len(comm_ids)):
                    cj = comm_ids[j]
                    xj, yj = meta_xy[cj]
                    dx = xj - xi
                    dy = yj - yi
                    dist = (dx * dx + dy * dy) ** 0.5
                    min_dist = community_radius[ci] + community_radius[cj] + 0.8
                    if dist < min_dist:
                        moved = True
                        if dist <= 1e-9:
                            ux, uy = 1.0, 0.0
                        else:
                            ux, uy = dx / dist, dy / dist
                        push = 0.5 * (min_dist - max(dist, 1e-9))
                        meta_xy[ci][0] -= ux * push
                        meta_xy[ci][1] -= uy * push
                        meta_xy[cj][0] += ux * push
                        meta_xy[cj][1] += uy * push
            if not moved:
                break

        final_pos: dict = {}
        for comm_idx, nodes in enumerate(communities):
            local_pos = local_layouts[comm_idx]
            max_radius = 0.0
            for x, y in local_pos.values():
                max_radius = max(max_radius, float((x * x + y * y) ** 0.5))
            if max_radius <= 0:
                max_radius = 1.0

            cx, cy = meta_xy.get(comm_idx, [0.0, 0.0])
            radius = community_radius[comm_idx]
            for node, (x, y) in local_pos.items():
                final_pos[node] = {
                    "x": float(cx + radius * (x / max_radius)),
                    "y": float(cy + radius * (y / max_radius)),
                }

        return final_pos

    def _ensure_cluster_cache(self) -> dict[str, dict]:
        """Return mutable in-memory cluster-artifact cache.

        The cache is keyed by metric name and stores artifacts generated by
        :meth:`compute_network_clusters` (tables, parameters, and layout).
        """
        if not hasattr(self, "_network_clusters"):
            self._network_clusters = {}
        return self._network_clusters

    def _export_cluster_node_id(self, node_id: str, query_id_set: set[str]) -> str:
        """Return export-safe node ID for cluster file output.

        Query IDs are preserved.  Library IDs are mapped through
        ``MolecularNetwork._export_node_id`` when present, so FlashEntropy
        internal indices can be exported as stable library identifiers (for
        example ``spectra_id``).
        """
        export_fn = getattr(self, "_export_node_id", None)
        if callable(export_fn):
            try:
                return str(export_fn(node_id, query_id_set))
            except Exception:
                pass
        return str(node_id)

    @staticmethod
    def _apply_query_connectivity_filters(
        G,
        *,
        nx,
        query_id_set: set,
        drop_components_without_queries: bool,
        drop_nodes_without_query_connection: bool,
    ) -> None:
        """Apply query-connectivity pruning in-place on graph ``G``.

        Parameters
        ----------
        G
            NetworkX graph to mutate.
        nx
            Imported ``networkx`` module, passed in to avoid repeated imports.
        query_id_set : set
            Set of query node IDs.
        drop_components_without_queries : bool
            Drop connected components that contain no query nodes.
        drop_nodes_without_query_connection : bool
            Keep only query nodes and one-hop neighbors of query nodes.
        """
        if drop_components_without_queries and G.number_of_nodes() > 0:
            if G.is_directed():
                components = list(nx.weakly_connected_components(G))
            else:
                components = list(nx.connected_components(G))

            nodes_to_remove = []
            for component in components:
                if component.isdisjoint(query_id_set):
                    nodes_to_remove.extend(component)

            if nodes_to_remove:
                G.remove_nodes_from(nodes_to_remove)

        if drop_nodes_without_query_connection and G.number_of_nodes() > 0 and query_id_set:
            query_nodes_in_graph = [n for n in query_id_set if G.has_node(n)]
            if query_nodes_in_graph:
                keep_nodes = set(query_nodes_in_graph)
                for qnode in query_nodes_in_graph:
                    keep_nodes.update(G.neighbors(qnode))

                if len(keep_nodes) < G.number_of_nodes():
                    nodes_to_remove = [n for n in G.nodes if n not in keep_nodes]
                    G.remove_nodes_from(nodes_to_remove)

    def compute_network_clusters(
        self,
        metric: str = "entropy_similarity",
        *,
        include_queries_only: bool = True,
        score_threshold: float | None = None,
        max_edges: int | None = 500,
        max_nodes: int | None = None,
        exclude_self: bool = True,
        directed: bool = False,
        drop_components_without_queries: bool = True,
        drop_nodes_without_query_connection: bool = True,
        cluster_method: str = "weighted_greedy_modularity",
        cluster_weight_attr: str = "score",
        cluster_super_threshold: int = 400,
        cluster_sparsify_top_k: int = 8,
        cluster_recursive_split: bool = True,
        cluster_min_size: int = 2,
        compute_layout: bool = True,
        layout_seed: int = 42,
    ) -> dict[str, Any]:
        """Compute and cache clustering artifacts for one similarity metric.

        This stage is separate from plotting.  It builds a graph from the
        selected similarity matrix, applies query-connectivity filters,
        computes communities, and stores node/community/edge/layout tables in
        an in-memory artifact cache for later plotting or export.

        Parameters
        ----------
        metric : str
            Similarity metric key from ``self.similarity_matrices``.
        include_queries_only : bool
            Use query-focused edge filtering during edge extraction.
        score_threshold : float or None
            Minimum similarity score to include; ``None`` uses
            ``_threshold_for(metric)``.
        max_edges : int or None
            Keep top-N edges by score after filtering.
        max_nodes : int or None
            Keep edges whose endpoints are within highest-degree N nodes.
        exclude_self : bool
            Drop self-loop edges.
        directed : bool
            Build directed graph when ``True``.
        drop_components_without_queries : bool
            Remove connected components with no query nodes.
        drop_nodes_without_query_connection : bool
            Keep only query nodes and one-hop query neighbors.
        cluster_method : str
            Community method; currently ``"weighted_greedy_modularity"``.
        cluster_weight_attr : str
            Edge attribute used as clustering weight.
        cluster_super_threshold : int
            Threshold above which recursive splitting can apply.
        cluster_sparsify_top_k : int
            Top-k edge sparsification for oversized community refinement.
        cluster_recursive_split : bool
            Enable recursive refinement of oversized communities.
        cluster_min_size : int
            Reserved for future use; kept for API stability.
        compute_layout : bool
            Compute and cache hierarchical node coordinates.
        layout_seed : int
            Random seed for deterministic layout generation.

        Returns
        -------
        dict
            Summary with ``metric``, ``n_nodes``, ``n_edges``, ``n_clusters``,
            and ``score_threshold``.

        Raises
        ------
        KeyError
            If *metric* is unknown.
        ValueError
            If *cluster_method* is unsupported.
        """
        if metric not in self.similarity_matrices:
            raise KeyError(
                f"Unknown metric '{metric}'. Available: {list(self.similarity_matrices.keys())}"
            )

        if cluster_method != "weighted_greedy_modularity":
            raise ValueError(
                "Unsupported cluster_method. "
                "Currently supported: ['weighted_greedy_modularity']"
            )

        nx = __import__("networkx")

        if score_threshold is None:
            score_threshold = self._threshold_for(metric)

        edges = self._edges_to_dataframe(
            metric=metric,
            score_threshold=score_threshold,
            include_queries_only=include_queries_only,
            exclude_self=exclude_self,
        )

        if max_nodes is not None and not edges.empty:
            degree = pd.concat([edges["id1"], edges["id2"]]).value_counts().head(max_nodes)
            top_nodes = set(degree.index)
            edges = edges[edges["id1"].isin(top_nodes) & edges["id2"].isin(top_nodes)]

        if max_edges is not None and not edges.empty:
            edges = edges.nlargest(max_edges, "score")

        G = nx.DiGraph() if directed else nx.Graph()
        query_id_set = set(self._all_query_ids)

        if not edges.empty:
            node_ids = pd.concat([edges["id1"], edges["id2"]]).drop_duplicates().tolist()
            for node_id in node_ids:
                G.add_node(node_id, node_type=("query" if node_id in query_id_set else "library"))
            for row in edges.itertuples(index=False):
                G.add_edge(row.id1, row.id2, score=float(row.score))

        self._apply_query_connectivity_filters(
            G,
            nx=nx,
            query_id_set=query_id_set,
            drop_components_without_queries=drop_components_without_queries,
            drop_nodes_without_query_connection=drop_nodes_without_query_connection,
        )

        communities: list[set] = []
        if G.number_of_nodes() > 0:
            communities = self._detect_weighted_communities(
                G.to_undirected(),
                weight_attr=cluster_weight_attr,
                super_threshold=cluster_super_threshold,
                sparsify_top_k=cluster_sparsify_top_k,
                recursive_split=cluster_recursive_split,
            )

        node_rows = []
        community_rows = []
        node_to_cluster: dict[str, str] = {}

        for comm_idx, nodes in enumerate(communities):
            cluster_id = f"C{comm_idx:03d}"
            query_count = sum(1 for n in nodes if n in query_id_set)
            community_rows.append(
                {
                    "cluster_id": cluster_id,
                    "community_level": 0,
                    "size": len(nodes),
                    "query_count_in_community": query_count,
                }
            )
            for node in sorted(nodes, key=str):
                node_to_cluster[str(node)] = cluster_id
                node_rows.append(
                    {
                        "node_id": str(node),
                        "cluster_id": cluster_id,
                        "community_level": 0,
                        "community_size": len(nodes),
                        "query_count_in_community": query_count,
                        "node_type": G.nodes[node].get("node_type", "unknown"),
                    }
                )

        # Ensure isolated nodes still receive a cluster assignment.
        if G.number_of_nodes() > 0:
            for node in G.nodes:
                node_s = str(node)
                if node_s in node_to_cluster:
                    continue
                cluster_id = f"C{len(community_rows):03d}"
                node_to_cluster[node_s] = cluster_id
                query_count = 1 if node in query_id_set else 0
                community_rows.append(
                    {
                        "cluster_id": cluster_id,
                        "community_level": 0,
                        "size": 1,
                        "query_count_in_community": query_count,
                    }
                )
                node_rows.append(
                    {
                        "node_id": node_s,
                        "cluster_id": cluster_id,
                        "community_level": 0,
                        "community_size": 1,
                        "query_count_in_community": query_count,
                        "node_type": G.nodes[node].get("node_type", "unknown"),
                    }
                )

        edge_rows = []
        for u, v, data in G.edges(data=True):
            u_s = str(u)
            v_s = str(v)
            edge_rows.append(
                {
                    "id1": u_s,
                    "id2": v_s,
                    "score": float(data.get("score", 0.0)),
                    "cluster_id_1": node_to_cluster.get(u_s),
                    "cluster_id_2": node_to_cluster.get(v_s),
                    "intra_cluster": node_to_cluster.get(u_s) == node_to_cluster.get(v_s),
                }
            )

        layout_df = pd.DataFrame(columns=["node_id", "x", "y"])
        if compute_layout and G.number_of_nodes() > 0:
            layout = self._compute_hierarchical_layout(
                G,
                communities,
                weight_attr=cluster_weight_attr,
                seed=layout_seed,
            )
            layout_df = pd.DataFrame(
                [
                    {"node_id": str(node_id), "x": coords["x"], "y": coords["y"]}
                    for node_id, coords in layout.items()
                ]
            )

        artifact = {
            "metric": metric,
            "params": {
                "include_queries_only": include_queries_only,
                "score_threshold": score_threshold,
                "max_edges": max_edges,
                "max_nodes": max_nodes,
                "exclude_self": exclude_self,
                "directed": directed,
                "drop_components_without_queries": drop_components_without_queries,
                "drop_nodes_without_query_connection": drop_nodes_without_query_connection,
                "cluster_method": cluster_method,
                "cluster_weight_attr": cluster_weight_attr,
                "cluster_super_threshold": cluster_super_threshold,
                "cluster_sparsify_top_k": cluster_sparsify_top_k,
                "cluster_recursive_split": cluster_recursive_split,
                "cluster_min_size": cluster_min_size,
                "compute_layout": compute_layout,
                "layout_seed": layout_seed,
            },
            "node_table": pd.DataFrame(node_rows),
            "community_table": pd.DataFrame(community_rows),
            "edge_table": pd.DataFrame(edge_rows),
            "layout_table": layout_df,
            "schema_version": 1,
        }

        self._ensure_cluster_cache()[metric] = artifact

        return {
            "metric": metric,
            "n_nodes": G.number_of_nodes(),
            "n_edges": G.number_of_edges(),
            "n_clusters": int(artifact["community_table"].shape[0]),
            "score_threshold": float(score_threshold),
        }

    def drop_network_clusters(self, metric: str | None = None) -> None:
        """Drop cached clustering artifacts while preserving similarity results.

        Parameters
        ----------
        metric : str or None
            Metric cache key to remove.  When ``None``, clear all cached
            clustering artifacts.
        """
        cache = self._ensure_cluster_cache()
        if metric is None:
            cache.clear()
        else:
            cache.pop(metric, None)

    def save_network_clusters(self, output_dir: str, metric: str, run_id: str | None = None) -> dict[str, str]:
        """Persist computed cluster artifacts to CSV files.

        Files written include node, community, edge, layout, and manifest CSVs.
        Library node IDs are exported using stable identifiers when available
        (for example ``spectra_id``) instead of internal index IDs.

        Parameters
        ----------
        output_dir : str
            Output directory for CSV artifacts.
        metric : str
            Metric cache key to persist.
        run_id : str or None
            Optional suffix to disambiguate multiple runs.

        Returns
        -------
        dict of str to str
            Paths for ``nodes``, ``communities``, ``edges``, ``layout``,
            and ``manifest`` files.

        Raises
        ------
        RuntimeError
            If no cached clusters exist for *metric*.
        """
        cache = self._ensure_cluster_cache()
        artifact = cache.get(metric)
        if artifact is None:
            raise RuntimeError(
                f"No clusters available for metric '{metric}'. Call compute_network_clusters(...) first."
            )

        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        suffix = f"_{run_id}" if run_id else ""
        prefix = f"{metric}_clusters{suffix}"

        node_path = out_dir / f"{prefix}_nodes.csv"
        community_path = out_dir / f"{prefix}_communities.csv"
        edge_path = out_dir / f"{prefix}_edges.csv"
        manifest_path = out_dir / f"{prefix}_manifest.csv"
        layout_path = out_dir / f"{prefix}_layout.csv"

        query_id_set = set(getattr(self, "_all_query_ids", []) or [])

        node_export_df = artifact["node_table"].copy()
        if "node_id" in node_export_df.columns:
            node_export_df["node_id"] = node_export_df["node_id"].astype(str).map(
                lambda node_id: self._export_cluster_node_id(node_id, query_id_set)
            )

        edge_export_df = artifact["edge_table"].copy()
        if "id1" in edge_export_df.columns:
            edge_export_df["id1"] = edge_export_df["id1"].astype(str).map(
                lambda node_id: self._export_cluster_node_id(node_id, query_id_set)
            )
        if "id2" in edge_export_df.columns:
            edge_export_df["id2"] = edge_export_df["id2"].astype(str).map(
                lambda node_id: self._export_cluster_node_id(node_id, query_id_set)
            )

        layout_export_df = artifact["layout_table"].copy()
        if "node_id" in layout_export_df.columns:
            layout_export_df["node_id"] = layout_export_df["node_id"].astype(str).map(
                lambda node_id: self._export_cluster_node_id(node_id, query_id_set)
            )

        node_export_df.to_csv(node_path, index=False)
        artifact["community_table"].to_csv(community_path, index=False)
        edge_export_df.to_csv(edge_path, index=False)
        layout_export_df.to_csv(layout_path, index=False)

        pd.DataFrame(
            [
                {
                    "schema_version": artifact.get("schema_version", 1),
                    "metric": metric,
                    "run_id": run_id or "",
                    "params_json": json.dumps(artifact.get("params", {}), sort_keys=True),
                }
            ]
        ).to_csv(manifest_path, index=False)

        return {
            "nodes": str(node_path),
            "communities": str(community_path),
            "edges": str(edge_path),
            "layout": str(layout_path),
            "manifest": str(manifest_path),
        }

    def load_network_clusters(self, output_dir: str, metric: str, run_id: str | None = None) -> dict[str, Any]:
        """Load cluster artifacts from CSV files into in-memory cache.

        Parameters
        ----------
        output_dir : str
            Directory containing saved cluster CSV artifacts.
        metric : str
            Metric key to load.
        run_id : str or None
            Optional run suffix used when saving.

        Returns
        -------
        dict
            Summary with ``metric``, ``n_nodes``, ``n_edges``, and
            ``n_clusters``.

        Raises
        ------
        FileNotFoundError
            If required artifact files are missing.
        """
        out_dir = Path(output_dir)
        suffix = f"_{run_id}" if run_id else ""
        prefix = f"{metric}_clusters{suffix}"

        node_path = out_dir / f"{prefix}_nodes.csv"
        community_path = out_dir / f"{prefix}_communities.csv"
        edge_path = out_dir / f"{prefix}_edges.csv"
        manifest_path = out_dir / f"{prefix}_manifest.csv"
        layout_path = out_dir / f"{prefix}_layout.csv"

        if not node_path.exists() or not community_path.exists() or not edge_path.exists():
            raise FileNotFoundError(
                f"Missing cluster artifact files for metric '{metric}' in '{output_dir}'."
            )

        params = {}
        schema_version = 1
        if manifest_path.exists():
            mdf = pd.read_csv(manifest_path)
            if not mdf.empty:
                schema_version = int(mdf.iloc[0].get("schema_version", 1))
                params_json = mdf.iloc[0].get("params_json", "{}")
                try:
                    params = json.loads(params_json) if isinstance(params_json, str) else {}
                except Exception:
                    params = {}

        artifact = {
            "metric": metric,
            "params": params,
            "schema_version": schema_version,
            "node_table": pd.read_csv(node_path),
            "community_table": pd.read_csv(community_path),
            "edge_table": pd.read_csv(edge_path),
            "layout_table": pd.read_csv(layout_path) if layout_path.exists() else pd.DataFrame(columns=["node_id", "x", "y"]),
        }

        self._ensure_cluster_cache()[metric] = artifact
        return {
            "metric": metric,
            "n_nodes": int(artifact["node_table"].shape[0]),
            "n_edges": int(artifact["edge_table"].shape[0]),
            "n_clusters": int(artifact["community_table"].shape[0]),
        }

    def plot_network(
        self,
        metric: str = "entropy_similarity",
        out_path: str = "network.html",
        *,
        include_queries_only: bool = True,
        score_threshold: float | None = None,
        max_edges: int | None = 500,
        max_nodes: int | None = None,
        exclude_self: bool = True,
        directed: bool = False,
        color_map: dict[str, str] | None = None,
        sigma_options: dict | None = None,
        metadata: dict[str, dict] | None = None,
        library_label_field: str | Sequence[str] | None = None,
        library_node_attrs: Sequence[str] | None = None,
        drop_components_without_queries: bool = True,
        drop_nodes_without_query_connection: bool = True,
        layout_seed: int = 42,
        bypass_clustering: bool = False,
        separate_communities: bool = True,
    ) -> str:
        """Render network HTML from precomputed clusters, or raw graph when bypassing.

        Parameters
        ----------
        separate_communities : bool
            When plotting from precomputed clusters (``bypass_clustering=False``),
            remove inter-community edges so each community renders as a separate
            subnetwork.  Default ``True``.
        """
        if metric not in self.similarity_matrices:
            raise KeyError(
                f"Unknown metric '{metric}'. Available: {list(self.similarity_matrices.keys())}"
            )

        try:
            import networkx as nx
            from ipysigma import Sigma
        except ImportError as exc:
            raise ImportError(
                "networkx and ipysigma are required for plot_network(). "
                "Install with: pip install networkx ipysigma"
            ) from exc

        cache = self._ensure_cluster_cache()
        artifact = cache.get(metric)
        if not bypass_clustering and artifact is None:
            raise RuntimeError(
                f"No clusters available for metric '{metric}'. "
                "Call compute_network_clusters(...) first, or set bypass_clustering=True."
            )

        if score_threshold is None:
            score_threshold = self._threshold_for(metric)

        if bypass_clustering:
            edges = self._edges_to_dataframe(
                metric=metric,
                score_threshold=score_threshold,
                include_queries_only=include_queries_only,
                exclude_self=exclude_self,
            )
        else:
            edges = artifact["edge_table"][["id1", "id2", "score"]].copy()
            edges = edges[edges["score"] >= score_threshold]

        if max_nodes is not None and not edges.empty:
            degree = pd.concat([edges["id1"], edges["id2"]]).value_counts().head(max_nodes)
            top_nodes = set(degree.index)
            edges = edges[edges["id1"].isin(top_nodes) & edges["id2"].isin(top_nodes)]

        if max_edges is not None and not edges.empty:
            edges = edges.nlargest(max_edges, "score")

        G: nx.Graph = nx.DiGraph() if directed else nx.Graph()

        palette = dict(self._DEFAULT_COLOR_MAP)
        if color_map:
            palette.update(color_map)

        query_id_set = set(self._all_query_ids)

        if not edges.empty:
            node_ids = pd.concat([edges["id1"], edges["id2"]]).drop_duplicates().tolist()

            for node_id in node_ids:
                node_kind = "query" if node_id in query_id_set else "library"
                lib_record = self._library_record_for_node(node_id) if node_kind == "library" else None

                label = (
                    self._resolve_library_label(
                        node_id=node_id,
                        record=lib_record,
                        library_label_field=library_label_field,
                    )
                    if node_kind == "library"
                    else (node_id[:20] + "..." if len(node_id) > 20 else node_id)
                )

                attrs: dict = {
                    "label": label,
                    "node_type": node_kind,
                }

                if node_kind == "library" and lib_record and library_node_attrs:
                    for field in library_node_attrs:
                        if field in lib_record:
                            attrs[field] = self._as_node_attr_value(lib_record.get(field))

                if metadata:
                    attrs.update(metadata.get(node_id, {}))

                G.add_node(node_id, **attrs)

            for row in edges.itertuples(index=False):
                G.add_edge(row.id1, row.id2, score=float(row.score))

        if bypass_clustering:
            self._apply_query_connectivity_filters(
                G,
                nx=nx,
                query_id_set=query_id_set,
                drop_components_without_queries=drop_components_without_queries,
                drop_nodes_without_query_connection=drop_nodes_without_query_connection,
            )

        # Default node styling: query nodes are visually emphasised.
        # Keep this in graph attributes so callers can still override with
        # sigma_options if they prefer another style.
        for node in G.nodes:
            node_type = G.nodes[node].get("node_type", "unknown")
            if node_type == "query":
                G.nodes[node]["viz_size"] = 2
                G.nodes[node]["viz_border_color"] = "#000000"
                G.nodes[node]["viz_border_size"] = 3
            else:
                G.nodes[node]["viz_size"] = 1
                G.nodes[node]["viz_border_color"] = "#00000000"
                G.nodes[node]["viz_border_size"] = 0

        community_groups: list[set] = []
        if not bypass_clustering and artifact is not None and not artifact["node_table"].empty:
            node_table = artifact["node_table"].set_index("node_id")
            for node in list(G.nodes):
                node_s = str(node)
                if node_s not in node_table.index:
                    continue
                row = node_table.loc[node_s]
                G.nodes[node]["community_id"] = row.get("cluster_id")
                G.nodes[node]["community_level"] = int(row.get("community_level", 0))
                G.nodes[node]["community_size"] = int(row.get("community_size", 1))
                G.nodes[node]["query_count_in_community"] = int(row.get("query_count_in_community", 0))

            by_cluster = node_table.reset_index().groupby("cluster_id")["node_id"]
            community_groups = [set(group.tolist()) for _, group in by_cluster]

            if separate_communities:
                inter_edges = []
                for u, v in G.edges():
                    cu = G.nodes[u].get("community_id")
                    cv = G.nodes[v].get("community_id")
                    if cu is not None and cv is not None and cu != cv:
                        inter_edges.append((u, v))
                if inter_edges:
                    G.remove_edges_from(inter_edges)

        layout: dict | None = None
        if G.number_of_nodes() > 0:
            if (
                not bypass_clustering
                and artifact is not None
                and not separate_communities
                and not artifact["layout_table"].empty
            ):
                layout = {
                    str(row.node_id): {"x": float(row.x), "y": float(row.y)}
                    for row in artifact["layout_table"].itertuples(index=False)
                    if G.has_node(str(row.node_id))
                }
                if len(layout) != G.number_of_nodes():
                    pos = nx.spring_layout(
                        G,
                        seed=layout_seed,
                        weight="score",
                        k=2.0 / max(G.number_of_nodes() ** 0.5, 1),
                    )
                    for node, (x, y) in pos.items():
                        layout.setdefault(str(node), {"x": float(x), "y": float(y)})
            elif not bypass_clustering and artifact is not None and community_groups:
                present_groups = []
                for group in community_groups:
                    nodes_present = {n for n in group if G.has_node(n)}
                    if nodes_present:
                        present_groups.append(nodes_present)
                layout = self._compute_hierarchical_layout(
                    G,
                    present_groups,
                    weight_attr="score",
                    seed=layout_seed,
                )
            else:
                pos = nx.spring_layout(
                    G,
                    seed=layout_seed,
                    weight="score",
                    k=2.0 / max(G.number_of_nodes() ** 0.5, 1),
                )
                layout = {str(n): {"x": float(x), "y": float(y)} for n, (x, y) in pos.items()}

        clustered_mode = not bypass_clustering and artifact is not None and G.number_of_nodes() > 0
        write_kwargs: dict = {
            "layout": layout,
            "node_color": "community_id" if clustered_mode else "node_type",
            "node_size": "viz_size",
            "node_size_range": (8, 18),
            "node_border_color": "viz_border_color",
            "node_border_size": "viz_border_size",
            "node_border_size_range": (0, 4),
            "default_edge_type": "line",
            "default_edge_color": "#aaa",
            "fullscreen": True,
        }
        if not clustered_mode:
            write_kwargs["node_color_palette"] = palette
        if sigma_options:
            write_kwargs.update(sigma_options)

        out = Path(out_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        Sigma.write_html(G, str(out), **write_kwargs)

        return str(out)
