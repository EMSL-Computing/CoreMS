"""Visualization mixin for molecular networking objects.

This module provides :class:`NetworkVisualizeMixin`, which adds interactive
HTML network plotting via ``networkx`` and ``ipysigma`` for classes exposing:

- ``similarity_matrices``: ``dict[str, SimilarityMatrix]``
- ``_all_query_ids``: ``list[str]``
- ``_threshold_for(metric) -> float``
"""

from __future__ import annotations

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
    ) -> str:
        """Render similarity network to an interactive HTML file.

        Uses ``networkx`` as the graph container and ``ipysigma`` (backed by
        sigma.js / WebGL) for rendering.  The output is a self-contained HTML
        file suitable for direct browser use or embedding in a web application.

        Parameters
        ----------
        metric : str
            Similarity metric to visualize.  Must be a key in
            ``self.similarity_matrices``.
        out_path : str
            Output HTML file path.  Parent directories are created if needed.
        include_queries_only : bool
            When ``True`` (default), keep only edges where at least one
            endpoint is a query ID, plus library–library edges between nodes
            that are both connected to at least one query node.
        score_threshold : float or None
            Minimum edge score to include.  If ``None``, the per-metric
            threshold from ``_threshold_for(metric)`` is used.
        max_edges : int or None
            If set, retain only the top-N edges ranked by score after all
            other filters have been applied.  Default 500.
        max_nodes : int or None
            If set, restrict the graph to the N highest-degree nodes and keep
            only edges where both endpoints are in that set (coherent
            subgraph).  Applied before ``max_edges``.  Default ``None``
            (no node cap).
        exclude_self : bool
            Drop self-loop edges (id1 == id2).  Default ``True``.
        directed : bool
            When ``True``, build a ``networkx.DiGraph`` instead of an
            undirected ``Graph``.  Default ``False``.
        color_map : dict or None
            Override node fill colors.  Keys are ``"query"``, ``"library"``,
            and/or ``"unknown"``; values are CSS color strings.  Unspecified
            keys fall back to ``_DEFAULT_COLOR_MAP``.
        sigma_options : dict or None
            Extra keyword arguments forwarded verbatim to
            ``Sigma.write_html()``.  Use this to customise rendering, e.g.::

                sigma_options={
                    "node_size": "degree",          # size nodes by degree
                    "node_metrics": ["louvain"],     # compute + color by community
                    "node_color": "louvain",
                    "default_edge_type": "curve",
                    "node_size_range": (3, 20),
                }

            Any key that overlaps with a positional argument of this method
            (e.g. ``"node_color_palette"``) will be forwarded as-is, allowing
            full control over ipysigma's API.  See the ipysigma documentation
            for the complete list of accepted parameters.  Example overrides::

                sigma_options={
                    # colour nodes by a metadata field instead of query/library
                    "node_color": "my_field",
                    # switch back to curved edges
                    "default_edge_type": "curve",
                    # run ForceAtlas2 in-browser for N seconds instead of
                    # using the pre-computed spring layout
                    "layout": None,
                    "start_layout": 5,
                }

        metadata : dict or None
            Optional mapping of node ID → dict of extra attributes.  Each
            key–value pair in the inner dict is added as a node attribute on
            the NetworkX graph and will appear as a labeled row in ipysigma's
            hover panel.  Example::

                metadata = {
                    "spec_001": {"name": "Compound A", "mz": 312.1},
                }

        library_label_field : str or sequence of str, or None
            Field name(s) to use for library-node labels from the FlashEntropy
            record dict (``self.fe_lib[int(node_id)]``).  If a sequence is
            provided, fields are tried in order until a non-empty value is
            found.  If ``None`` (default), tries
            ``("compound_name", "name", "spectra_id", "id")`` before
            falling back to the library node ID.
        library_node_attrs : sequence of str or None
            Optional FlashEntropy record keys to copy onto each library node
            as node attributes (shown in ipysigma hover panels).  Example::

                library_node_attrs=("compound_name", "spectra_id", "precursor_mz")

            Any values that are not JSON-scalar types are stringified.
        drop_components_without_queries : bool
            When ``True`` (default), remove any connected component that does
            not contain at least one query spectrum.  This suppresses isolated
            library-only islands (for example, small doublet/triplet groups)
            from the rendered network.
        drop_nodes_without_query_connection : bool
            When ``True`` (default), keep only query nodes plus nodes that are
            directly adjacent (one hop) to at least one query node.  This
            prevents chained propagation through library-only paths.

        Returns
        -------
        str
            The absolute path to the saved HTML file.

        Raises
        ------
        KeyError
            If *metric* is not present in ``self.similarity_matrices``.
        ImportError
            If ``networkx`` or ``ipysigma`` are not installed.
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

        # ── Build NetworkX graph ───────────────────────────────────────────
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

                if node_kind == "library":
                    label = self._resolve_library_label(
                        node_id=node_id,
                        record=lib_record,
                        library_label_field=library_label_field,
                    )
                else:
                    label = node_id[:20] + "..." if len(node_id) > 20 else node_id

                # Start with display / type attributes
                attrs: dict = {
                    "label": label,
                    "node_type": node_kind,
                }

                if node_kind == "library" and lib_record and library_node_attrs:
                    for field in library_node_attrs:
                        if field in lib_record:
                            attrs[field] = self._as_node_attr_value(lib_record.get(field))

                # Spread any caller-supplied metadata as individual node
                # attributes so ipysigma surfaces them in its hover panel.
                if metadata:
                    attrs.update(metadata.get(node_id, {}))

                G.add_node(node_id, **attrs)

            # Edges are added unweighted — topology only.  The score is
            # retained as a data attribute for hover display but is not used
            # for visual weight or layout, so highly-connected nodes cluster
            # together purely by connection count.
            for row in edges.itertuples(index=False):
                G.add_edge(row.id1, row.id2, score=float(row.score))

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

        # ── Pre-compute layout ─────────────────────────────────────────────
        # Fruchterman-Reingold (spring) treats all edges as equal-weight
        # springs, which pulls highly-connected nodes toward their neighbours.
        # Baking positions into the file means the layout is deterministic and
        # renders instantly — no animation needed on load.
        #
        # k controls the natural spring length (higher → nodes spread more).
        # For sparse similarity graphs a larger k prevents overcrowding.
        layout: dict | None = None
        if G.number_of_nodes() > 0:
            pos = nx.spring_layout(G, seed=42, k=2.0 / max(G.number_of_nodes() ** 0.5, 1))
            layout = {n: {"x": float(x), "y": float(y)} for n, (x, y) in pos.items()}

        # ── Render to HTML ─────────────────────────────────────────────────
        # Defaults prioritise topology-based layout and clean visuals:
        #   - layout=<pre-computed>  deterministic positions from spring_layout
        #   - default_edge_type="line"  straight edges; less visual noise than
        #     curves when nodes are well-separated by the layout
        #   - node_size="degree"  larger nodes have more connections, making
        #     hubs immediately obvious
        #   - edge_size / edge_size_range omitted  all edges same thickness
        #     so high-degree hubs aren't obscured by thick-edge clutter
        #
        # Any of these can be overridden via sigma_options.
        write_kwargs: dict = {
            "layout": layout,
            "node_color": "node_type",
            "node_color_palette": palette,
            "node_size": "degree",
            "node_size_range": (4, 24),
            "default_edge_type": "line",
            "default_edge_color": "#aaa",
            "fullscreen": True,
        }
        if sigma_options:
            write_kwargs.update(sigma_options)

        out = Path(out_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        Sigma.write_html(G, str(out), **write_kwargs)

        return str(out)
