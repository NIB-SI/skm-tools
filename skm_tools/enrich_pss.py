"""
PSS network enrichment with CKN.

PART ONE — build & persist
    Enriches the curated PSS network with nodes and edges from CKN that
    are not already present in PSS, restricted to CKN content that
    touches at least one PSS node. Built once (e.g. whenever PSS or CKN
    is updated) and saved to disk.

PART TWO — enrichment analysis
    Two exports characterising what the enrichment introduced:
    - Export A: PSS pairs newly joined by a direct CKN edge
    - Export B: PSS pairs newly joined through exactly one CKN intermediate
    Plus `extract_pair_subnetwork` to visualise any one pair in context.

PART THREE — node neighbourhood inspection
    Pick any node (random or chosen) and inspect how much CKN enriched
    its immediate neighbourhood.

Main entry points
------------------
build_enriched_pss(pss, ckn)
    Build the enriched graph. Every node/edge carries an `origin`
    attribute: "PSS" (from the original PSS network) or "CKN" (added
    from CKN during enrichment). PSS edges are never replaced by CKN
    edges between the same two nodes.

save_enriched(enriched, out_dir) / load_enriched(in_dir)
    Persist / reload the enriched graph (GraphML + TSV node/edge tables
    + a diff table listing only what CKN contributed).

get_enrichment_summary(pss, enriched)
    Quick before/after node and edge counts.

visualize_enriched_network(enriched, network_name=...)
    Push the enriched graph (or a subgraph of it) to Cytoscape, coloured
    by `origin` (PSS = blue, CKN = orange). Used to sanity-check PART ONE.

find_new_direct_edges(enriched, pss, max_rank=None)
    Export A.

find_new_intermediate_pairs(enriched, pss, max_rank=None)
    Export B.

extract_pair_subnetwork(enriched, pss, a, b, intermediate=None, max_neighbours=None)
    Small subnetwork around an enriched pair, for Cytoscape inspection.

pick_random_node(enriched, node_type=None, seed=None)
    Pick a node to inspect.

extract_node_neighbourhood(enriched, node, radius=1, max_neighbours=None)
    Ego-graph around one node.

get_node_neighbourhood_enrichment(pss, enriched, node, radius=1, max_neighbours=None)
    How much CKN enriched this node's neighbourhood.
"""

from pathlib import Path
import networkx as nx
import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
#  Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _pss_has_edge(pss, u, v):
    """True iff PSS has an edge between u and v in either direction."""
    return pss.has_edge(u, v) or pss.has_edge(v, u)


def _clean_value(v):
    """Coerce an attribute value to a GraphML-friendly scalar."""
    if v is None:
        return ""
    if isinstance(v, (list, tuple, set)):
        return ",".join(str(x) for x in v)
    if isinstance(v, (str, int, float, bool)):
        return v
    return str(v)


def _clean_attrs(d):
    """Clean a dict of attributes for GraphML export."""
    return {k: _clean_value(v) for k, v in d.items() if v is not None}


def _enriched_has_ckn_edge(enriched, u, v):
    """True iff enriched already has an origin='CKN' edge between u and v."""
    for a, b in [(u, v), (v, u)]:
        if enriched.has_edge(a, b):
            for ed in enriched[a][b].values():
                if ed.get("origin") == "CKN":
                    return True
    return False


def _undirected_neighbours(G, n):
    """Successors ∪ predecessors for directed graphs; neighbours otherwise."""
    if G.is_directed():
        return set(G.successors(n)) | set(G.predecessors(n))
    return set(G.neighbors(n))


def _ckn_edge_ranks(enriched, u, v):
    """Return list of ranks for all origin='CKN' edges between u and v."""
    ranks = []
    for a, b in [(u, v), (v, u)]:
        if enriched.has_edge(a, b):
            for ed in enriched[a][b].values():
                if ed.get("origin") == "CKN":
                    r = ed.get("rank")
                    if r not in (None, ""):
                        try:
                            ranks.append(int(r))
                        except (ValueError, TypeError):
                            pass
    return ranks


# ─────────────────────────────────────────────────────────────────────────────
#  PART ONE — Build
# ─────────────────────────────────────────────────────────────────────────────

def build_enriched_pss(pss, ckn, verbose=True):
    """
    Build an enriched PSS network by adding CKN nodes and edges that touch
    at least one PSS node.

    Scope rule:
      A CKN edge (u, v) is added iff
      - u or v is a PSS node (both or one), OR
      - u and v are both CKN nodes that were *already* pulled in by the
        rule above (so local CKN connectivity around PSS is preserved).

    Every node and edge carries an `origin` attribute:
      "PSS"  — came from the original PSS network
      "CKN"  — added from CKN during enrichment

    A PSS edge is NEVER overwritten by a CKN edge. If CKN has an edge
    between two PSS nodes that already have a PSS edge, it is not added.

    Parameters
    ----------
    pss : nx.DiGraph | nx.MultiDiGraph
        Original PSS network.
    ckn : nx.DiGraph | nx.MultiDiGraph
        Full CKN network (edges should carry a `rank` attribute).
    verbose : bool
        Print progress.

    Returns
    -------
    nx.MultiDiGraph
        The enriched graph. Every node has `origin`, every edge has
        `origin` and (for CKN edges) `rank`.
    """
    enriched = nx.MultiDiGraph()

    # Step 1 — copy PSS nodes/edges
    if verbose:
        print("[1/3] Copying PSS nodes and edges …")
    for n, d in pss.nodes(data=True):
        enriched.add_node(n, origin="PSS", **_clean_attrs(d))
    for u, v, d in pss.edges(data=True):
        enriched.add_edge(u, v, origin="PSS", **_clean_attrs(d))

    n_pss_nodes = enriched.number_of_nodes()
    n_pss_edges = enriched.number_of_edges()
    if verbose:
        print(f"      PSS: {n_pss_nodes:,} nodes, {n_pss_edges:,} edges")

    # Step 2 — add CKN edges touching PSS (and the CKN nodes they bring in)
    if verbose:
        print("[2/3] Adding CKN edges that touch PSS …")
    pss_nodes = set(pss.nodes())
    ckn_added_nodes = set()
    n_ckn_edges_added = 0

    for u, v, d in ckn.edges(data=True):
        u_in_pss = u in pss_nodes
        v_in_pss = v in pss_nodes

        if not (u_in_pss or v_in_pss):
            continue

        # PSS keeps priority: skip if PSS already has an edge here
        if u_in_pss and v_in_pss and _pss_has_edge(pss, u, v):
            continue

        if not u_in_pss and u not in enriched:
            enriched.add_node(u, origin="CKN", **_clean_attrs(ckn.nodes[u]))
            ckn_added_nodes.add(u)
        if not v_in_pss and v not in enriched:
            enriched.add_node(v, origin="CKN", **_clean_attrs(ckn.nodes[v]))
            ckn_added_nodes.add(v)

        enriched.add_edge(u, v, origin="CKN", **_clean_attrs(d))
        n_ckn_edges_added += 1

    if verbose:
        print(f"      + {len(ckn_added_nodes):,} CKN nodes, "
              f"{n_ckn_edges_added:,} CKN edges (touching PSS)")

    # Step 3 — add CKN edges between two CKN nodes that both got added above
    if verbose:
        print("[3/3] Adding CKN edges between pulled-in CKN nodes …")
    n_ckn_internal = 0
    for u, v, d in ckn.edges(data=True):
        if u in ckn_added_nodes and v in ckn_added_nodes:
            enriched.add_edge(u, v, origin="CKN", **_clean_attrs(d))
            n_ckn_internal += 1

    if verbose:
        print(f"      + {n_ckn_internal:,} CKN-CKN edges")
        print("-" * 58)
        print(f"  Enriched: {enriched.number_of_nodes():,} nodes "
              f"({n_pss_nodes:,} PSS + {len(ckn_added_nodes):,} CKN)")
        print(f"  Enriched: {enriched.number_of_edges():,} edges "
              f"({n_pss_edges:,} PSS + "
              f"{n_ckn_edges_added + n_ckn_internal:,} CKN)")
        print("-" * 58)

    return enriched


# ─────────────────────────────────────────────────────────────────────────────
#  Save / load
# ─────────────────────────────────────────────────────────────────────────────

def save_enriched(enriched, out_dir):
    """
    Save enriched graph to disk.

    Files written in `out_dir`:
      enriched.graphml      — full graph
      enriched_nodes.tsv    — node id + origin + attributes
      enriched_edges.tsv    — source, target, origin, rank, …
      enrichment_diff.tsv   — only CKN-added nodes and edges
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    G = nx.MultiDiGraph()
    for n, d in enriched.nodes(data=True):
        G.add_node(n, **_clean_attrs(d))
    for u, v, d in enriched.edges(data=True):
        G.add_edge(u, v, **_clean_attrs(d))
    nx.write_graphml(G, out / "enriched.graphml")

    node_rows = [{"node": n, **d} for n, d in enriched.nodes(data=True)]
    pd.DataFrame(node_rows).to_csv(out / "enriched_nodes.tsv", sep="\t", index=False)

    edge_rows = [{"source": u, "target": v, **d}
                 for u, v, d in enriched.edges(data=True)]
    pd.DataFrame(edge_rows).to_csv(out / "enriched_edges.tsv", sep="\t", index=False)

    diff_nodes = [r for r in node_rows if r.get("origin") == "CKN"]
    diff_edges = [r for r in edge_rows if r.get("origin") == "CKN"]
    diff_rows = (
        [{"kind": "node", "source": r["node"], "target": "",
          **{k: v for k, v in r.items() if k != "node"}}
         for r in diff_nodes]
        + [{"kind": "edge", **r} for r in diff_edges]
    )
    pd.DataFrame(diff_rows).to_csv(out / "enrichment_diff.tsv", sep="\t", index=False)

    print(f"OK: Saved enriched network to {out}/")
    print(f"  enriched.graphml      ({enriched.number_of_nodes():,} nodes, "
          f"{enriched.number_of_edges():,} edges)")
    print(f"  enriched_nodes.tsv    ({len(node_rows):,} rows)")
    print(f"  enriched_edges.tsv    ({len(edge_rows):,} rows)")
    print(f"  enrichment_diff.tsv   ({len(diff_rows):,} rows: "
          f"{len(diff_nodes):,} CKN nodes + {len(diff_edges):,} CKN edges)")


def load_enriched(in_dir):
    """Load the enriched graph from the GraphML saved by `save_enriched`."""
    path = Path(in_dir) / "enriched.graphml"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found — run build_enriched_pss first.")
    G = nx.read_graphml(path)
    if not isinstance(G, nx.MultiDiGraph):
        M = nx.MultiDiGraph()
        M.add_nodes_from(G.nodes(data=True))
        M.add_edges_from(G.edges(data=True))
        G = M
    for u, v, k, d in G.edges(keys=True, data=True):
        if "rank" in d:
            try:
                d["rank"] = int(d["rank"])
            except (ValueError, TypeError):
                pass
    print(f"OK: Loaded enriched: {G.number_of_nodes():,} nodes, "
          f"{G.number_of_edges():,} edges")
    return G


# ─────────────────────────────────────────────────────────────────────────────
#  Quick stats (for testing PART ONE)
# ─────────────────────────────────────────────────────────────────────────────

def get_enrichment_summary(pss, enriched):
    """
    Before/after summary comparing the original PSS network to the
    enriched network.

    Returns
    -------
    dict with keys:
      pss_nodes, pss_edges, enriched_nodes, enriched_edges,
      ckn_nodes_added, ckn_edges_added
    """
    ckn_nodes_added = sum(1 for _, d in enriched.nodes(data=True)
                          if d.get("origin") == "CKN")
    ckn_edges_added = sum(1 for _, _, d in enriched.edges(data=True)
                          if d.get("origin") == "CKN")
    return {
        "pss_nodes": pss.number_of_nodes(),
        "pss_edges": pss.number_of_edges(),
        "enriched_nodes": enriched.number_of_nodes(),
        "enriched_edges": enriched.number_of_edges(),
        "ckn_nodes_added": ckn_nodes_added,
        "ckn_edges_added": ckn_edges_added,
    }


def print_enrichment_summary(pss, enriched):
    """Pretty-print the output of `get_enrichment_summary`."""
    s = get_enrichment_summary(pss, enriched)
    W = 24
    print("=" * (W + 26))
    print(f"{'':>{W}}  {'PSS (before)':>12}  {'Enriched (after)':>17}")
    print("-" * (W + 26))
    print(f"{'Nodes':>{W}}  {s['pss_nodes']:>12,}  {s['enriched_nodes']:>17,}")
    print(f"{'  of which CKN-added':>{W}}  {'—':>12}  {s['ckn_nodes_added']:>17,}")
    print(f"{'Edges':>{W}}  {s['pss_edges']:>12,}  {s['enriched_edges']:>17,}")
    print(f"{'  of which CKN-added':>{W}}  {'—':>12}  {s['ckn_edges_added']:>17,}")
    print("=" * (W + 26))
    return s


# ─────────────────────────────────────────────────────────────────────────────
#  Cytoscape visualisation (for testing PART ONE)
# ─────────────────────────────────────────────────────────────────────────────

def visualize_enriched_network(graph, network_name="enriched_pss", base_style="ckn"):
    """
    Send a graph to Cytoscape, coloured by `origin` (PSS = blue, CKN = orange).

    Use this on the full enriched network (small PSS only) or on an
    ego-graph / sample subgraph for large networks — Cytoscape struggles
    with tens of thousands of nodes.

    Requires Cytoscape running and py4cytoscape installed.

    Returns
    -------
    int
        Cytoscape SUID of the created network.
    """
    import py4cytoscape as p4c
    from .cytoscape_utils import apply_builtin_style

    G = nx.MultiDiGraph()
    for n, d in graph.nodes(data=True):
        G.add_node(n, **_clean_attrs(d))
    for u, v, k, d in graph.edges(keys=True, data=True):
        G.add_edge(u, v, key=k, **_clean_attrs(d))

    print(f"Sending to Cytoscape: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges …")

    suid = p4c.create_network_from_networkx(G, title=network_name)
    apply_builtin_style(suid, base_style)

    current_style = p4c.get_current_style(network=suid)
    origin_style = f"{current_style}-origin"
    if origin_style in p4c.get_visual_style_names():
        p4c.delete_visual_style(origin_style)
    p4c.copy_visual_style(current_style, origin_style)

    p4c.set_edge_color_mapping(
        "origin",
        table_column_values=["PSS", "CKN"],
        colors=["#2171B5", "#E6550D"],
        mapping_type="d",
        style_name=origin_style,
    )
    p4c.set_node_border_color_mapping(
        "origin",
        table_column_values=["PSS", "CKN"],
        colors=["#2171B5", "#E6550D"],
        mapping_type="d",
        style_name=origin_style,
    )
    p4c.set_node_border_width_mapping(
        "origin",
        table_column_values=["PSS", "CKN"],
        widths=[5, 2],
        mapping_type="d",
        style_name=origin_style,
    )

    p4c.set_visual_style(origin_style, network=suid)
    p4c.layout_network("force-directed", network=suid)
    p4c.fit_content(network=suid)

    print(f"OK: Network '{network_name}' created in Cytoscape (SUID: {suid})")
    return suid


# ─────────────────────────────────────────────────────────────────────────────
#  PART TWO — Enrichment analysis
# ─────────────────────────────────────────────────────────────────────────────

def find_new_direct_edges(enriched, pss, max_rank=None):
    """
    Export A — PSS pairs now directly connected by a CKN edge that was not
    in the original PSS.

    Parameters
    ----------
    enriched : nx.MultiDiGraph
    pss : nx.DiGraph | nx.MultiDiGraph
    max_rank : int or None
        If set, keep only pairs where the best CKN rank is <= max_rank.

    Returns
    -------
    pandas.DataFrame
        Columns: u, v, rank, direction
    """
    pss_nodes = set(pss.nodes())
    seen = set()
    rows = []

    for u, v, d in enriched.edges(data=True):
        if d.get("origin") != "CKN":
            continue
        if u not in pss_nodes or v not in pss_nodes:
            continue
        if _pss_has_edge(pss, u, v):
            continue

        pair = tuple(sorted([u, v]))
        if pair in seen:
            continue
        seen.add(pair)

        ranks = _ckn_edge_ranks(enriched, u, v)
        best = min(ranks) if ranks else None
        if max_rank is not None and (best is None or best > max_rank):
            continue

        has_uv = (enriched.has_edge(u, v)
                  and any(ed.get("origin") == "CKN" for ed in enriched[u][v].values()))
        has_vu = (enriched.has_edge(v, u)
                  and any(ed.get("origin") == "CKN" for ed in enriched[v][u].values()))
        direction = "both" if (has_uv and has_vu) else ("u->v" if has_uv else "v->u")

        rows.append({"u": u, "v": v, "rank": best, "direction": direction})

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["rank", "u", "v"]).reset_index(drop=True)
    return df


def find_new_intermediate_pairs(enriched, pss, max_rank=None):
    """
    Export B — PSS pairs connected through exactly ONE CKN intermediate,
    where no PSS edge existed between them before.

    Parameters
    ----------
    enriched : nx.MultiDiGraph
    pss : nx.DiGraph | nx.MultiDiGraph
    max_rank : int or None
        If set, keep only paths where max(rank_ux, rank_xv) <= max_rank.

    Returns
    -------
    pandas.DataFrame
        Columns: u, v, intermediate, rank_ux, rank_xv, path_rank
        One row per (unordered_pair, intermediate) triple.
    """
    from itertools import combinations

    pss_nodes = set(pss.nodes())
    ckn_nodes_in_enriched = {n for n, d in enriched.nodes(data=True)
                             if d.get("origin") == "CKN"}

    ckn_to_pss_neighbours = {x: set() for x in ckn_nodes_in_enriched}
    for x in ckn_nodes_in_enriched:
        for nb in _undirected_neighbours(enriched, x):
            if nb in pss_nodes:
                ckn_to_pss_neighbours[x].add(nb)

    seen = set()
    rows = []

    for x, partners in ckn_to_pss_neighbours.items():
        if len(partners) < 2:
            continue
        for a, b in combinations(sorted(partners), 2):
            if _pss_has_edge(pss, a, b):
                continue
            key = (a, b, x)
            if key in seen:
                continue
            seen.add(key)

            ranks_ax = _ckn_edge_ranks(enriched, a, x)
            ranks_xb = _ckn_edge_ranks(enriched, x, b)
            if not ranks_ax or not ranks_xb:
                continue
            rank_ax = min(ranks_ax)
            rank_xb = min(ranks_xb)
            path_rank = max(rank_ax, rank_xb)
            if max_rank is not None and path_rank > max_rank:
                continue

            rows.append({
                "u": a, "v": b, "intermediate": x,
                "rank_ux": rank_ax, "rank_xv": rank_xb, "path_rank": path_rank,
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["path_rank", "u", "v", "intermediate"]).reset_index(drop=True)
    return df


def extract_pair_subnetwork(enriched, pss, a, b, intermediate=None, max_neighbours=None):
    """
    Extract a small subnetwork around an enriched pair (a, b) for
    Cytoscape inspection.

    Node set:
      - a, b              -> role = "focus"
      - intermediate      -> role = "intermediate"   (if provided)
      - 1-hop PSS neighbours of a and b (from the ORIGINAL pss network)
                          -> role = "context"
        capped at `max_neighbours` per focus node if specified.

    Edge set:
      All edges from `enriched` between nodes in the chosen set. New CKN
      edges among the focus/intermediate nodes are tagged
      edge_status="new_from_ckn"; everything else is "pss_original" or
      "ckn_context".

    Returns
    -------
    nx.MultiDiGraph
    """
    def _pss_neighbours_capped(node):
        if node not in pss:
            return set()
        nbrs = _undirected_neighbours(pss, node)
        nbrs.discard(a)
        nbrs.discard(b)
        if intermediate is not None:
            nbrs.discard(intermediate)
        if max_neighbours is not None and len(nbrs) > max_neighbours:
            nbrs = set(sorted(nbrs)[:max_neighbours])
        return nbrs

    context_a = _pss_neighbours_capped(a)
    context_b = _pss_neighbours_capped(b)

    node_set = {a, b} | context_a | context_b
    if intermediate is not None:
        node_set.add(intermediate)
    node_set &= set(enriched.nodes())

    sub = enriched.subgraph(node_set).copy()

    for n in sub.nodes():
        if n in (a, b):
            sub.nodes[n]["role"] = "focus"
        elif intermediate is not None and n == intermediate:
            sub.nodes[n]["role"] = "intermediate"
        else:
            sub.nodes[n]["role"] = "context"

    focal = {a, b}
    if intermediate is not None:
        focal.add(intermediate)
    for u, v, k, d in sub.edges(keys=True, data=True):
        origin = d.get("origin", "PSS")
        if origin == "CKN" and u in focal and v in focal:
            sub.edges[u, v, k]["edge_status"] = "new_from_ckn"
        else:
            sub.edges[u, v, k]["edge_status"] = "pss_original" if origin == "PSS" else "ckn_context"

    return sub


def visualize_pair_subnetwork(graph, network_name="enriched_pair", base_style="ckn"):
    """
    Send a pair-subnetwork (from `extract_pair_subnetwork`) to Cytoscape:
    PSS/CKN colouring on nodes and edges, focus/intermediate/context node
    border widths, and new_from_ckn edges drawn extra thick.

    Requires Cytoscape running and py4cytoscape installed.
    """
    import py4cytoscape as p4c
    from .cytoscape_utils import apply_builtin_style

    G = nx.MultiDiGraph()
    for n, d in graph.nodes(data=True):
        G.add_node(n, **_clean_attrs(d))
    for u, v, k, d in graph.edges(keys=True, data=True):
        G.add_edge(u, v, key=k, **_clean_attrs(d))

    print(f"Sending to Cytoscape: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges …")

    suid = p4c.create_network_from_networkx(G, title=network_name)
    apply_builtin_style(suid, base_style)

    current_style = p4c.get_current_style(network=suid)
    pair_style = f"{current_style}-enrich-pair"
    if pair_style in p4c.get_visual_style_names():
        p4c.delete_visual_style(pair_style)
    p4c.copy_visual_style(current_style, pair_style)

    p4c.set_edge_color_mapping(
        "origin", table_column_values=["PSS", "CKN"],
        colors=["#2171B5", "#E6550D"], mapping_type="d", style_name=pair_style,
    )
    p4c.set_edge_line_width_mapping(
        "edge_status",
        table_column_values=["pss_original", "ckn_context", "new_from_ckn"],
        widths=[2.0, 2.0, 6.0],
        mapping_type="d", style_name=pair_style,
    )
    p4c.set_node_border_color_mapping(
        "origin", table_column_values=["PSS", "CKN"],
        colors=["#2171B5", "#E6550D"], mapping_type="d", style_name=pair_style,
    )
    p4c.set_node_border_width_mapping(
        "role", table_column_values=["focus", "intermediate", "context"],
        widths=[10, 6, 2], mapping_type="d", style_name=pair_style,
    )

    p4c.set_visual_style(pair_style, network=suid)
    p4c.layout_network("force-directed", network=suid)
    p4c.fit_content(network=suid)

    print(f"OK: Network '{network_name}' created in Cytoscape (SUID: {suid})")
    return suid


# ─────────────────────────────────────────────────────────────────────────────
#  PART THREE — Node neighbourhood inspection
# ─────────────────────────────────────────────────────────────────────────────
#
# Recommended radius
# -------------------
# radius=1 (default) — direct neighbours only. Answers "did CKN change who
#     this exact node is connected to?" Cheap, always readable, the right
#     default for a random spot-check.
# radius=2 — neighbours-of-neighbours. Useful to see pathway crosstalk
#     opened up by a CKN intermediate, but a handful of CKN hub nodes
#     bridge 100s-1000s of PSS pairs, so radius=2 can blow up fast. Always
#     pair radius=2 with `max_neighbours` (e.g. 15-25) to keep it readable.
# radius>=3 is not recommended for visual inspection of a single node.

def pick_random_node(enriched, node_type=None, seed=None):
    """
    Pick a random node from the enriched network.

    Parameters
    ----------
    enriched : nx.MultiDiGraph
    node_type : "PSS" | "CKN" | None
        Restrict to nodes with this origin. None = any node.
    seed : int, optional
        For a reproducible pick.

    Returns
    -------
    str
    """
    import random
    rng = random.Random(seed)

    if node_type is None:
        candidates = list(enriched.nodes())
    else:
        candidates = [n for n, d in enriched.nodes(data=True)
                     if d.get("origin") == node_type]

    if not candidates:
        raise ValueError(f"No nodes found for node_type={node_type!r}")

    return rng.choice(candidates)


def extract_node_neighbourhood(enriched, node, radius=1, max_neighbours=None):
    """
    Ego-graph around `node` in the enriched network.

    Parameters
    ----------
    enriched : nx.MultiDiGraph
    node : str
    radius : int
        See module-level "Recommended radius" note above. Default 1.
    max_neighbours : int, optional
        Cap on neighbours kept per hop (random sample) — keeps hub-node
        neighbourhoods visualisable. Recommended whenever radius > 1.

    Returns
    -------
    nx.MultiDiGraph
        Subgraph with `role` = "focus" (the node itself) or "context".
    """
    if node not in enriched:
        raise ValueError(f"{node!r} not found in the enriched network.")

    undirected = enriched.to_undirected(as_view=True)

    if max_neighbours is None:
        node_set = set(nx.ego_graph(undirected, node, radius=radius).nodes())
    else:
        import random
        rng = random.Random(0)
        node_set = {node}
        frontier = {node}
        for _ in range(radius):
            next_frontier = set()
            for n in frontier:
                nbrs = set(undirected.neighbors(n)) - node_set
                if len(nbrs) > max_neighbours:
                    nbrs = set(rng.sample(sorted(nbrs), max_neighbours))
                next_frontier |= nbrs
            node_set |= next_frontier
            frontier = next_frontier

    sub = enriched.subgraph(node_set).copy()
    for n in sub.nodes():
        sub.nodes[n]["role"] = "focus" if n == node else "context"

    return sub


def get_node_neighbourhood_enrichment(pss, enriched, node, radius=1, max_neighbours=None):
    """
    Quantify how much CKN enriched the neighbourhood of a single node.

    Returns
    -------
    dict with keys:
      node, radius, pss_neighbours, enriched_neighbours,
      ckn_added_neighbours, pct_new
    """
    sub = extract_node_neighbourhood(enriched, node, radius=radius,
                                     max_neighbours=max_neighbours)
    enriched_neighbours = sub.number_of_nodes() - 1

    if node in pss:
        pss_ego = nx.ego_graph(pss.to_undirected(as_view=True), node, radius=radius)
        pss_neighbours = pss_ego.number_of_nodes() - 1
    else:
        pss_neighbours = 0

    ckn_added = sum(1 for n in sub.nodes() if sub.nodes[n].get("origin") == "CKN")
    pct_new = (ckn_added / enriched_neighbours * 100) if enriched_neighbours else 0.0

    return {
        "node": node,
        "radius": radius,
        "pss_neighbours": pss_neighbours,
        "enriched_neighbours": enriched_neighbours,
        "ckn_added_neighbours": ckn_added,
        "pct_new": pct_new,
    }


def print_node_neighbourhood_summary(stats):
    """Pretty-print the output of `get_node_neighbourhood_enrichment`."""
    print(f"Node: {stats['node']}  (radius={stats['radius']})")
    print(f"  PSS neighbours (before)      : {stats['pss_neighbours']:,}")
    print(f"  Enriched neighbours (after)  : {stats['enriched_neighbours']:,}")
    print(f"  of which new from CKN        : {stats['ckn_added_neighbours']:,} "
          f"({stats['pct_new']:.1f}%)")


def visualize_node_neighbourhood(graph, node, network_name=None, base_style="ckn"):
    """
    Push a node-neighbourhood subgraph (from `extract_node_neighbourhood`)
    to Cytoscape: PSS/CKN colouring (blue/orange) plus a thick border on
    the focal node.

    Requires Cytoscape running and py4cytoscape installed.
    """
    import py4cytoscape as p4c
    from .cytoscape_utils import apply_builtin_style

    if network_name is None:
        network_name = f"neighbourhood_{node}"

    G = nx.MultiDiGraph()
    for n, d in graph.nodes(data=True):
        G.add_node(n, **_clean_attrs(d))
    for u, v, k, d in graph.edges(keys=True, data=True):
        G.add_edge(u, v, key=k, **_clean_attrs(d))

    print(f"Sending to Cytoscape: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges …")

    suid = p4c.create_network_from_networkx(G, title=network_name)
    apply_builtin_style(suid, base_style)

    current_style = p4c.get_current_style(network=suid)
    nb_style = f"{current_style}-neighbourhood"
    if nb_style in p4c.get_visual_style_names():
        p4c.delete_visual_style(nb_style)
    p4c.copy_visual_style(current_style, nb_style)

    p4c.set_edge_color_mapping(
        "origin", table_column_values=["PSS", "CKN"],
        colors=["#2171B5", "#E6550D"], mapping_type="d", style_name=nb_style,
    )
    p4c.set_node_border_color_mapping(
        "origin", table_column_values=["PSS", "CKN"],
        colors=["#2171B5", "#E6550D"], mapping_type="d", style_name=nb_style,
    )
    p4c.set_node_border_width_mapping(
        "role", table_column_values=["focus", "context"],
        widths=[10, 2], mapping_type="d", style_name=nb_style,
    )

    p4c.set_visual_style(nb_style, network=suid)
    p4c.layout_network("force-directed", network=suid)
    p4c.fit_content(network=suid)

    print(f"OK: Network '{network_name}' created in Cytoscape (SUID: {suid})")
    return suid
