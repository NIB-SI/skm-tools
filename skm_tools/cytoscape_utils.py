'''Cytoscape automation: load networks, apply the SKM styles, highlight nodes, edges and
paths, create subnetworks, and export images.

Requires the ``cytoscape`` extra (``pip install skm-tools[cytoscape]``) and a running
Cytoscape (https://cytoscape.org) for py4cytoscape to talk to.

Functions acting on a network take it as ``network`` (a SUID or a name; default: the
current network in Cytoscape). Needs Cytoscape 3.10 or later.

py4cytoscape's console output (e.g. the text of errors that are handled, such as the
retries in :func:`clone_network`) is silenced when this module is imported, see
:func:`silence_py4cytoscape`.

The helpers that don't need Cytoscape (matching image files to nodes, chart definitions)
are in :mod:`skm_tools.node_images`, and also available here.

Nodes and edges are matched to Cytoscape by name: a node's ``name`` is its networkx node
id, and an edge's ``name`` is ``"source (interaction) target"``, where ``interaction`` is
the edge's ``interaction`` attribute (``"interacts with"`` if it has none), as set by
:func:`load_network` / ``py4cytoscape.create_network_from_networkx``. Parallel edges
between the same two nodes are matched together.
'''

import logging
import re
import sys
import time
import uuid
from collections import defaultdict
from pathlib import Path

import pandas as pd
import py4cytoscape as p4c

from . import resources
from .neighbors import neighborhood_nodes
from .node_images import (  # noqa: F401 (also here, for convenience)
    chart_column,
    match_files_to_nodes,
    node_file_key,
    unique_image_copies,
)
from .paths import path_edges
from .utils import as_list, contrast_color

logger = logging.getLogger(__name__)


_P4C_QUIET = True
_P4C_ORIGINAL = {}


def _quiet_wrapper(name, original):
    def wrapper(text, *args, **kwargs):
        if _P4C_QUIET:
            return text  # narrate returns its text
        return original(text, *args, **kwargs)
    wrapper._skm_tools_wrapped = True
    return wrapper


def silence_py4cytoscape(silence=True):
    '''Silence (or restore) py4cytoscape's console output.

    py4cytoscape prints the text of every error it raises, also of errors that are caught
    and handled (e.g. the retries in :func:`clone_network`), and prints progress messages
    in notebooks. Silenced when this module is imported. The errors are still raised,
    usually with the same text. But when Cytoscape's error response isn't JSON,
    py4cytoscape prints the response and raises a plain ``HTTPError`` ("500 Server Error
    ... for url ..."): silenced, Cytoscape's message is then only in py4cytoscape's log file
    (``logs/py4cytoscape.log`` in the working directory, not affected by this). When
    debugging, call ``silence_py4cytoscape(False)``.

    Parameters
    ----------
    silence : bool
        True (default): no console output from py4cytoscape. False: as py4cytoscape does
        by default.
    '''
    global _P4C_QUIET
    _P4C_QUIET = silence

    # show_error (error text) and narrate (progress) are imported into many py4cytoscape
    # modules, so wrap them in each
    for module_name, module in list(sys.modules.items()):
        if not module_name.startswith("py4cytoscape") or module is None:
            continue
        for name in ("show_error", "narrate"):
            f = getattr(module, name, None)
            if f is None or getattr(f, "_skm_tools_wrapped", False):
                continue
            _P4C_ORIGINAL.setdefault(name, f)
            setattr(module, name, _quiet_wrapper(name, _P4C_ORIGINAL[name]))

    summary = logging.getLogger("py4...S")
    if silence:
        _P4C_ORIGINAL.setdefault("summary_level", summary.level)
        summary.setLevel(logging.CRITICAL)
    elif "summary_level" in _P4C_ORIGINAL:
        summary.setLevel(_P4C_ORIGINAL["summary_level"])


silence_py4cytoscape()


# ---------------------------------------------------------------------------
# Matching networkx nodes/edges to Cytoscape SUIDs
# ---------------------------------------------------------------------------

def _node_suids(nodes, network=None, warn=True):
    '''Cytoscape SUID -> name of the nodes named `nodes`; names not in the network are
    skipped, with a warning.'''
    names = {str(n) for n in as_list(nodes) or []}
    table = p4c.tables.get_table_columns(table="node", columns=["name"], network=network)
    found = {int(suid): name for suid, name in table["name"].items() if name in names}
    missing = names - set(found.values())
    if missing and warn:
        shown = ", ".join(sorted(missing)[:10]) + (", ..." if len(missing) > 10 else "")
        logger.warning("%d of %d nodes not in the Cytoscape network: %s",
                       len(missing), len(names), shown)
    return found


def _match_edge_names(edge_names, edge_pairs):
    '''For Cytoscape edge names "u (interaction) v", the ones that link a pair in `edge_pairs`.

    Parameters
    ----------
    edge_names : dict
        SUID -> edge name.
    edge_pairs : iterable of tuple
        (u, v) pairs, or (u, v, key) (the key is ignored).

    Returns
    -------
    dict
        SUID -> (u, v) for every matching edge.
    '''
    wanted = {(str(e[0]), str(e[1])) for e in edge_pairs}
    sources = {u for u, _ in wanted}

    matches = {}
    for suid, name in edge_names.items():
        if not isinstance(name, str):
            continue
        # node names can contain " (" and ") ", so try every split point
        i = name.find(" (")
        while i != -1:
            u = name[:i]
            if u in sources:
                rest = name[i + 2:]
                j = rest.find(") ")
                while j != -1:
                    if (u, rest[j + 2:]) in wanted:
                        matches[suid] = (u, rest[j + 2:])
                        break
                    j = rest.find(") ", j + 1)
            i = name.find(" (", i + 1)
    return matches


def _edge_suids(edge_pairs, network=None):
    '''Cytoscape SUID -> (u, v) of the edges between each (u, v) pair (all parallel edges).'''
    edge_pairs = list(edge_pairs)
    if not edge_pairs:
        return {}
    table = p4c.tables.get_table_columns(table="edge", columns=["name"], network=network)
    return {int(suid): pair for suid, pair in _match_edge_names(table["name"].to_dict(), edge_pairs).items()}


# ---------------------------------------------------------------------------
# Loading and styling
# ---------------------------------------------------------------------------

def set_style(style, network=None):
    '''Apply a visual style to one network.

    Cytoscape also applies the style to the network(s) selected in its *Network* panel,
    so ``py4cytoscape.set_visual_style`` can change other networks too. This makes
    `network` the current (selected) network first.

    Parameters
    ----------
    style : str
        Visual style name.
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).
    '''
    if network is not None:
        p4c.set_current_network(network)
    p4c.set_visual_style(style, network=network)


def load_network(g, title, collection=None, style=None):
    '''Load a networkx graph into Cytoscape.

    Parameters
    ----------
    g : networkx.Graph
        Network to load. Node ids become the Cytoscape ``name``; node and edge attributes
        become table columns.
    title : str
        Network name in Cytoscape.
    collection : str, optional
        Collection to add the network to (default: a new collection named `title`).
    style : {"skm", "skm-reactions"}, optional
        Built-in SKM style to apply (see :func:`apply_builtin_style`); without it,
        Cytoscape's default style.

    Returns
    -------
    int
        SUID of the new Cytoscape network.
    '''
    suid = p4c.networks.create_network_from_networkx(g, title=title, collection=collection)
    if style:
        apply_builtin_style(style, network=suid)
    return suid


def apply_builtin_style(style="skm", network=None):
    '''Apply one of the bundled SKM visual styles to a Cytoscape network.

    The styles are imported into Cytoscape on first use (from ``skm-styles.xml``).

    Parameters
    ----------
    style : {"skm", "skm-reactions"}
        ``"skm"`` (style *SKM*, default): for CKN, the PSS interaction network and the PSS
        gene networks. Nodes by ``node_type`` (the PSS classes, in the PSS Explorer's
        colours), labelled with ``display_label``; edges coloured and with arrows by
        ``interaction``, dashed for mutual (``directed`` False) edges, and thicker for better
        supported ones (``rank``).
        ``"skm-reactions"`` (style *SKM-reactions*): for the PSS reaction graph, with
        reactions as nodes, and edges by ``edge_type`` (activation, inhibition, substrate,
        product, ...) labelled with their ``role``.
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).

    Returns
    -------
    str
        The style name (*SKM* or *SKM-reactions*).
    '''
    key = style.lower()
    if key not in resources.BUILTIN_STYLES:
        raise ValueError(f"apply_builtin_style expects one of {sorted(resources.BUILTIN_STYLES)}, "
                         f"not {style!r}.")
    style_name = resources.BUILTIN_STYLES[key]

    if style_name not in p4c.styles.get_visual_style_names():
        p4c.import_visual_styles(resources.get_style_xml_path())

    set_style(style_name, network)
    return style_name


# ---------------------------------------------------------------------------
# Highlighting (style bypasses)
# ---------------------------------------------------------------------------

_NODE_BYPASSES = ("NODE_FILL_COLOR", "NODE_LABEL_COLOR", "NODE_BORDER_PAINT", "NODE_BORDER_WIDTH",
                  "NODE_HEIGHT", "NODE_WIDTH")
_EDGE_BYPASSES = ("EDGE_STROKE_UNSELECTED_PAINT", "EDGE_UNSELECTED_PAINT", "EDGE_WIDTH")


def highlight_nodes(nodes, color=None, label_color=None, border_color=None, border_width=None,
                    node_height=None, node_width=None, network=None):
    '''Highlight nodes with style bypasses.

    Only the given properties are changed. Bypasses stay when the visual style changes;
    see :func:`clear_highlights` to undo them.

    Parameters
    ----------
    nodes : node or iterable of nodes
        Nodes to highlight (names not in the network are skipped, with a warning).
    color, label_color, border_color : str, optional
        Fill, label and border colours, as hex (``"#FF0000"``).
    border_width, node_height, node_width : float, optional
        Sizes in pixels.
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).

    Returns
    -------
    list of str
        The names of the highlighted nodes.
    '''
    found = _node_suids(nodes, network=network)
    if not found:
        return []

    # to set back, see bug https://github.com/cytoscape/py4cytoscape/issues/114
    og_style = p4c.styles.get_current_style(network)

    bypasses = [
        (color, p4c.style_bypasses.set_node_color_bypass),
        (label_color, p4c.style_bypasses.set_node_label_color_bypass),
        (border_color, p4c.style_bypasses.set_node_border_color_bypass),
        (border_width, p4c.style_bypasses.set_node_border_width_bypass),
        (node_height, p4c.style_bypasses.set_node_height_bypass),
        (node_width, p4c.style_bypasses.set_node_width_bypass),
    ]
    for value, set_bypass in bypasses:
        if value is not None:
            set_bypass(list(found), value, network=network)

    set_style(og_style, network)

    return list(found.values())


def highlight_edges(edges, color, skip_edges=None, edge_line_width=10, directed=True, network=None):
    '''Highlight edges with style bypasses (colour and line width).

    Parameters
    ----------
    edges : iterable of tuple
        Edges as (u, v) or (u, v, key) tuples, e.g. from :func:`skm_tools.paths.path_edges`.
        All parallel edges between u and v are highlighted.
    color : str
        Edge colour, as hex.
    skip_edges : iterable of tuple, optional
        (u, v) pairs not to highlight, e.g. edges already highlighted in another colour.
    edge_line_width : float
        Line width in pixels (default 10).
    directed : bool
        If False, also highlight the edges v -> u (e.g. for paths from
        ``get_paths(..., directed=False)``, which can use an edge against its direction).
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).

    Returns
    -------
    list of tuple
        The (u, v) pairs highlighted: found in the network, and not skipped.
    '''
    skip = {(e[0], e[1]) for e in (skip_edges or [])}
    pairs = [(e[0], e[1]) for e in edges]
    if not directed:
        pairs += [(v, u) for u, v in pairs]
    pairs = [p for p in dict.fromkeys(pairs) if p not in skip]

    found = _edge_suids(pairs, network=network)
    if found:
        p4c.style_bypasses.set_edge_line_width_bypass(list(found), edge_line_width, network=network)
        p4c.style_bypasses.set_edge_color_bypass(list(found), color, network=network)

    # as given (the Cytoscape names are strings)
    by_name = {(str(u), str(v)): (u, v) for u, v in pairs}
    return list(dict.fromkeys(by_name[pair] for pair in found.values()))


def highlight_path(nodes, color, skip_nodes=None, skip_edges=None, label_color="white",
                   border_color="black", border_width=10, edge_line_width=10, directed=True,
                   network=None):
    '''Highlight the nodes and edges of a path.

    To highlight several paths in different colours without painting over earlier ones,
    pass the nodes and edges returned by the earlier calls as `skip_nodes` and `skip_edges`.

    Parameters
    ----------
    nodes : list
        The path, as a list of nodes.
    color : str
        Node fill and edge colour, as hex.
    skip_nodes : iterable, optional
        Nodes not to highlight.
    skip_edges : iterable of tuple, optional
        (u, v) pairs not to highlight.
    label_color, border_color : str
        Node label and border colours (default white and black).
    border_width, edge_line_width : float
        Node border and edge widths in pixels (default 10).
    directed : bool
        If False, also highlight the edge v -> u of each step u, v (for paths from
        ``get_paths(..., directed=False)``).
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).

    Returns
    -------
    nodes : list
        The path's nodes that were highlighted.
    edges : list of tuple
        The path's (u, v) edges that were highlighted.
    '''
    nodes = list(nodes)
    skip_nodes = set(skip_nodes or [])
    to_highlight = [n for n in nodes if n not in skip_nodes]

    highlighted = []
    # the edges are highlighted even if all nodes already are (e.g. a path between two coloured nodes)
    if to_highlight:
        names = set(highlight_nodes(to_highlight, color=color, label_color=label_color,
                                    border_color=border_color, border_width=border_width,
                                    network=network))
        highlighted = [n for n in to_highlight if str(n) in names]

    edges = highlight_edges(list(zip(nodes, nodes[1:])), color, skip_edges=skip_edges,
                            edge_line_width=edge_line_width, directed=directed, network=network)

    return highlighted, edges


def clear_highlights(nodes=None, edges=None, network=None):
    '''Remove the highlights (style bypasses) of :func:`highlight_nodes`,
    :func:`highlight_edges` and :func:`highlight_path`.

    Cytoscape clears bypasses one node and one property at a time, which takes a while for
    many nodes: pass the nodes and edges that were highlighted (as returned by the highlight
    functions) rather than clearing everything.

    Parameters
    ----------
    nodes : iterable, optional
        Nodes to clear. Default: all nodes, if `edges` isn't given either.
    edges : iterable of tuple, optional
        (u, v) edges to clear. Default: all edges, if `nodes` isn't given either.
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).
    '''
    if nodes is None and edges is None:
        node_suids = list(p4c.tables.get_table_columns(table="node", columns=["name"], network=network).index)
        edge_suids = list(p4c.tables.get_table_columns(table="edge", columns=["name"], network=network).index)
    else:
        node_suids = list(_node_suids(nodes, network=network)) if nodes is not None else []
        edge_suids = list(_edge_suids(edges, network=network)) if edges is not None else []

    for prop in _NODE_BYPASSES if node_suids else ():
        p4c.style_bypasses.clear_node_property_bypass(node_suids, prop, network=network)
    for prop in _EDGE_BYPASSES if edge_suids else ():
        p4c.style_bypasses.clear_edge_property_bypass(edge_suids, prop, network=network)


def apply_shortest_paths_style(g, sources, path_lists, edge_colors=None, node_colors=None,
                               style_name=None, network=None):
    '''Colour the results of several path searches (e.g. one per source) in a new visual style.

    Adds node columns ``distance-to-target`` (the number of steps to the end of the path; on
    several paths, the smallest) and ``node-path-source`` (the first search whose paths
    have the node), and the edge column ``edge-priority`` (``"direct path (<source>)"``), then
    copies the network's current style to `style_name`, with colour mappings on them.

    Parameters
    ----------
    g : networkx.Graph
        Graph the paths were found in (to find their edges, see
        :func:`skm_tools.paths.path_edges`).
    sources : list
        The source of each search (its label in ``edge-priority``).
    path_lists : list of list
        Paths found by each search: one list of paths per source, in the same order.
    edge_colors : list of str, optional
        Hex colour for the edges of each source's paths, one per source.
    node_colors : list of str, optional
        Hex colours for a continuous node colour mapping on ``distance-to-target``, from the
        target (distance 0) to the farthest node. With one colour, the second is white.
    style_name : str, optional
        Name of the new style (default: ``"<current style>-shortest-paths"``; replaced if it
        exists).
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).

    Returns
    -------
    str
        The style name.
    '''
    og_style = p4c.styles.get_current_style(network)
    if style_name is None:
        style_name = og_style if og_style.endswith("-shortest-paths") else f"{og_style}-shortest-paths"
    if style_name != og_style:
        copy_style(og_style, style_name, overwrite=True)

    node_attributes = defaultdict(dict)
    edge_priority = {}

    for source, paths in zip(sources, path_lists):
        for p in paths:
            for i, node in enumerate(p):
                distance = len(p) - 1 - i
                d = node_attributes[str(node)]
                d["distance-to-target"] = min(distance, d.get("distance-to-target", distance))
                # the first search that reaches a node claims it
                d.setdefault("node-path-source", str(source))

        for suid in _edge_suids(path_edges(g, paths), network=network):
            edge_priority.setdefault(suid, {"edge-priority": f"direct path ({source})"})

    if node_attributes:
        p4c.tables.load_table_data(
            pd.DataFrame.from_dict(node_attributes, orient="index").rename_axis("name").reset_index(),
            data_key_column="name", table="node", table_key_column="name", network=network,
        )
    if edge_priority:
        p4c.tables.load_table_data(
            pd.DataFrame.from_dict(edge_priority, orient="index").rename_axis("SUID").reset_index(),
            data_key_column="SUID", table="edge", table_key_column="SUID", network=network,
        )

    if node_colors and node_attributes:
        node_colors = list(node_colors) + (["#FFFFFF"] if len(node_colors) == 1 else [])
        max_distance = max(d["distance-to-target"] for d in node_attributes.values()) or 1
        n = len(node_colors)
        values = [max_distance * i / (n - 1) for i in range(n)]

        p4c.style_mappings.set_node_color_mapping(
            "distance-to-target", table_column_values=values, colors=node_colors,
            mapping_type="c", style_name=style_name, network=network,
        )
        p4c.style_mappings.set_node_label_color_mapping(
            "distance-to-target", table_column_values=values,
            colors=[contrast_color(x) for x in node_colors],
            mapping_type="c", style_name=style_name, network=network,
        )

    if edge_colors:
        p4c.style_mappings.set_edge_color_mapping(
            "edge-priority",
            table_column_values=[f"direct path ({source})" for source in sources],
            colors=edge_colors, mapping_type="d", style_name=style_name, network=network,
        )

    set_style(style_name, network)
    return style_name


# ---------------------------------------------------------------------------
# Networks, collections and styles
# ---------------------------------------------------------------------------

_RETRIES = 20
_RETRY_WAIT = 0.5  # seconds


def clone_network(name=None, collection=None, network=None):
    '''Clone a Cytoscape network into a new collection, optionally renaming both.

    py4cytoscape's ``clone_network`` doesn't rename the new collection; this does.

    Parameters
    ----------
    name : str, optional
        Name of the clone.
    collection : str, optional
        Name of the clone's new collection.
    network : int or str, optional
        Cytoscape network to clone, SUID or name (default: the current network).

    Returns
    -------
    int
        SUID of the clone.
    '''
    suid = p4c.networks.clone_network(network=network)

    if collection:
        # py4cytoscape can't rename a collection: set the name in the collection's (root
        # network's) default table through CyREST
        collection_suid = p4c.collections.get_collection_suid(suid)
        p4c.commands.cyrest_put(
            f"collections/{collection_suid}/tables/default",
            body={"key": "SUID", "data": [{"SUID": collection_suid, "name": collection}]},
            require_json=False,
        )

    if name:
        # after renaming the collection (renaming the network first didn't work in the
        # original notebooks). Right after cloning, Cytoscape may not know the new network
        # yet ("unrecognized table entry"), so retry for a few seconds
        for attempt in range(_RETRIES):
            try:
                p4c.rename_network(name, network=suid)
                break
            except p4c.CyError:
                if attempt == _RETRIES - 1:
                    raise
                time.sleep(_RETRY_WAIT)

    return suid


def copy_style(style, new_style, networks=(), overwrite=False):
    '''Copy a visual style, and apply the copy to `networks`.

    Parameters
    ----------
    style : str
        Visual style to copy.
    new_style : str
        Name of the copy.
    networks : iterable of int or str, optional
        Networks to apply the copy to.
    overwrite : bool
        If a style `new_style` exists: replace it (True), or raise ValueError (False,
        default). Cytoscape itself would add a number to the name.

    Returns
    -------
    str
        `new_style`.
    '''
    if new_style in p4c.styles.get_visual_style_names():
        if not overwrite:
            raise ValueError(f"A style {new_style!r} exists already; use overwrite=True to replace it.")
        p4c.styles.delete_visual_style(new_style)
    p4c.copy_visual_style(style, new_style)
    if new_style not in p4c.styles.get_visual_style_names():
        raise RuntimeError(f"Cytoscape didn't name the copy of {style!r} {new_style!r}.")
    for network in networks:
        set_style(new_style, network)
    return new_style


def delete_other_networks(keep):
    '''Delete every Cytoscape network except `keep`, and make it the current network.

    Handy for starting over while debugging a notebook. Can't be undone (save the
    session first if needed).

    Parameters
    ----------
    keep : int
        SUID of the network to keep.

    Returns
    -------
    list of int
        SUIDs of the deleted networks.
    '''
    deleted = [x["suid"] for x in p4c.get_network_list(get_suids=True) if x["suid"] != keep]
    for suid in deleted:
        p4c.delete_network(suid)
    p4c.set_current_network(keep)
    return deleted


# ---------------------------------------------------------------------------
# Subnetworks
# ---------------------------------------------------------------------------

def subnetwork_edge_induced(edges, parent_suid, name="subnetwork (edge induced)"):
    '''New Cytoscape network with only the given edges (and their nodes).

    Parameters
    ----------
    edges : iterable of tuple
        Edges as (u, v) or (u, v, key) tuples; all parallel edges between u and v are included.
    parent_suid : int
        SUID of the Cytoscape network to take them from.
    name : str
        Name of the new network.

    Returns
    -------
    int
        SUID of the new network.
    '''
    edges = list(edges)
    nodes = list(dict.fromkeys(n for e in edges for n in e[:2]))

    return p4c.networks.create_subnetwork(
        nodes=list(_node_suids(nodes, network=parent_suid)),
        edges=list(_edge_suids(edges, network=parent_suid)),
        subnetwork_name=name,
        network=parent_suid,
        exclude_edges=True,
    )


def subnetwork_edge_induced_from_paths(g, paths, parent_suid, name="subnetwork (edge induced)"):
    '''New Cytoscape network with only the edges along `paths`.

    Parameters
    ----------
    g : networkx.Graph
        Graph the paths were found in (to find their edges, see
        :func:`skm_tools.paths.path_edges`).
    paths : list of list
        Paths as lists of nodes, e.g. from :func:`skm_tools.paths.get_paths`.
    parent_suid : int
        SUID of the Cytoscape network to take them from.
    name : str
        Name of the new network.

    Returns
    -------
    int
        SUID of the new network.
    '''
    return subnetwork_edge_induced(path_edges(g, paths), parent_suid, name=name)


def subnetwork_node_induced(nodes, parent_suid, name="subnetwork (node induced)"):
    '''New Cytoscape network with the given nodes and every edge between them.

    Node-induced: the result includes every edge of `parent_suid` between the given
    nodes, not just edges from whatever selection process produced `nodes` (e.g. a
    curated reaction/edge list). To keep an exact edge subset instead, use
    :func:`subnetwork_edge_induced`.

    Parameters
    ----------
    nodes : iterable
        Nodes to include (names not in the network are skipped, with a warning).
    parent_suid : int
        SUID of the Cytoscape network to take them from.
    name : str
        Name of the new network.

    Returns
    -------
    int
        SUID of the new network.
    '''
    return p4c.networks.create_subnetwork(
        nodes=list(_node_suids(nodes, network=parent_suid)),
        subnetwork_name=name,
        network=parent_suid
    )


def get_or_create_subnetwork(nodes, parent_suid, name):
    '''The network called `name` in the collection of `parent_suid`, created with
    :func:`subnetwork_node_induced` if it doesn't exist yet.

    For re-running notebook cells without making duplicate networks.

    Parameters
    ----------
    nodes : iterable
        Nodes to include, if the network is created.
    parent_suid : int
        SUID of the Cytoscape network to take them from.
    name : str
        Name of the network.

    Returns
    -------
    int
        SUID of the (existing or new) network.
    '''
    collection = p4c.get_collection_networks(p4c.get_collection_suid(parent_suid))
    for suid in collection:
        if p4c.get_network_name(suid) == name:
            return suid
    return subnetwork_node_induced(nodes, parent_suid, name=name)


def subnetwork_neighbours(g, nodes, parent_suid, depth=1, direction="both",
                          name="subnetwork (neighbours)"):
    '''New Cytoscape network with `nodes`, their neighbours, and every edge between them.

    The neighbours are found in `g` (see :func:`skm_tools.neighbors.neighborhood_nodes`),
    so the selection in Cytoscape doesn't change.

    Parameters
    ----------
    g : networkx.Graph
        The network loaded as `parent_suid` (or the part of it to search).
    nodes : iterable
        Nodes whose neighbourhood to include.
    parent_suid : int
        SUID of the Cytoscape network to take them from.
    depth : int or None
        Number of steps (default 1: first neighbours).
    direction : {"both", "out", "in"}
        For directed graphs: neighbours in both directions (default), downstream or upstream.
    name : str
        Name of the new network.

    Returns
    -------
    int
        SUID of the new network.
    '''
    neighbours = neighborhood_nodes(g, nodes, depth=depth, direction=direction)
    return subnetwork_node_induced(list(neighbours), parent_suid, name=name)


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

def layout_from_coords(positions, flip_y=True, scale=1.0, network=None):
    '''Place nodes at given coordinates (e.g. from a networkx or graphviz layout).

    Uses a temporary visual style that maps temporary node columns to the node positions,
    and removes both afterwards (also on errors); the network's style and columns are left
    as they were.

    Parameters
    ----------
    positions : dict or pandas.DataFrame
        ``{node: (x, y)}`` (as returned by the networkx layout functions), or a DataFrame
        indexed by node name with columns ``x`` and ``y``.
    flip_y : bool
        Mirror the y coordinates (default True): in networkx and graphviz y grows upwards,
        in Cytoscape downwards.
    scale : float
        Multiply the coordinates by this (default 1). networkx layouts are within -1 and 1,
        so use e.g. ``scale=500`` for them.
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).
    '''
    if isinstance(positions, pd.DataFrame):
        coords = positions[["x", "y"]].astype(float)
    else:
        coords = pd.DataFrame.from_dict({str(n): xy for n, xy in positions.items()},
                                        orient="index", columns=["x", "y"]).astype(float)
    coords.index = coords.index.map(str)
    coords = coords * scale
    if flip_y:
        coords["y"] = -coords["y"]

    tag = uuid.uuid4().hex[:8]
    x_column, y_column, tmp_style = f"skm_tools_x_{tag}", f"skm_tools_y_{tag}", f"skm-tools-layout-{tag}"
    current_style = p4c.styles.get_current_style(network=network)

    p4c.load_table_data(coords.rename(columns={"x": x_column, "y": y_column}),
                        table_key_column="name", network=network)
    try:
        p4c.copy_visual_style(current_style, tmp_style)
        set_style(tmp_style, network)
        p4c.update_style_mapping(tmp_style, p4c.map_visual_property('NODE_X_LOCATION', x_column, 'p'))
        p4c.update_style_mapping(tmp_style, p4c.map_visual_property('NODE_Y_LOCATION', y_column, 'p'))
    finally:
        set_style(current_style, network)
        if tmp_style in p4c.styles.get_visual_style_names():
            p4c.delete_visual_style(tmp_style)
        for column in (x_column, y_column):
            p4c.tables.delete_table_column(column, network=network)


# ---------------------------------------------------------------------------
# Node images and charts
# ---------------------------------------------------------------------------

# NODE_CUSTOMGRAPHICS_POSITION values: "<node anchor>,<graphic anchor>,<justification>,x,y"
IMAGE_POSITIONS = {
    "below": "S,N,c,0.00,0.00",
    "above": "N,S,c,0.00,0.00",
    "right": "E,W,c,0.00,0.00",
    "left": "W,E,c,0.00,0.00",
    "center": "C,C,c,0.00,0.00",
}


def load_node_images(images, column, network=None, unique_dir=None, to_png=False):
    '''Load image file locations into a node table column, for :func:`show_node_images`.

    Parameters
    ----------
    images : dict
        Node name -> image path (e.g. from :func:`skm_tools.node_images.match_files_to_nodes`).
        Nodes not in the network are ignored.
    column : str
        Node table column to load them into, e.g. ``"image_heat_Desiree"``.
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).
    unique_dir : str or pathlib.Path, optional
        Copy the images to this folder first, each to a new, unique file name (see
        :func:`skm_tools.node_images.unique_image_copies`), so that Cytoscape's image cache
        can't show an old or a different image instead.
    to_png : bool
        With `unique_dir`: convert SVG images to PNG.

    Returns
    -------
    dict
        Node -> the path loaded.
    '''
    if unique_dir is not None:
        images = unique_image_copies(images, unique_dir, to_png=to_png)
    images = {str(n): Path(p).absolute() for n, p in images.items()}
    if images:
        table = pd.DataFrame({column: {n: f"file:{p}" for n, p in images.items()}})
        p4c.load_table_data(table, table_key_column="name", network=network)
    return images


def show_node_images(style, column, slot=1, position="below", size=110):
    '''Show the images in a node table column (see :func:`load_node_images`) in a style.

    Maps the column to one of Cytoscape's custom graphics slots
    (``NODE_CUSTOMGRAPHICS_<slot>``), so a node can show several images (e.g. one per data
    type) in different slots and positions.

    Parameters
    ----------
    style : str
        Visual style to change. To show different columns (e.g. conditions) on copies of a
        network, use a copy of the style per column, see :func:`copy_style`.
    column : str
        Node table column with the image locations.
    slot : int
        Custom graphics slot, 1 to 9 (default 1).
    position : str
        ``"below"`` (default), ``"above"``, ``"left"``, ``"right"``, ``"center"``, or a
        Cytoscape custom graphics position (``"S,N,c,0.00,0.00"``: the image's north side
        on the node's south side, centred, no offset).
    size : float
        Image size (default 110). The images aren't resized with the node.
    '''
    p4c.style_dependencies.sync_node_custom_graphics_size(False, style_name=style)
    p4c.style_mappings.update_style_mapping(
        style,
        p4c.style_mappings.map_visual_property(
            visual_prop=f"NODE_CUSTOMGRAPHICS_{slot}", table_column=column, mapping_type="p"
        ),
    )
    p4c.style_defaults.set_visual_property_default(
        {"visualProperty": f"NODE_CUSTOMGRAPHICS_POSITION_{slot}",
         "value": IMAGE_POSITIONS.get(position, position)},
        style_name=style,
    )
    p4c.style_defaults.set_visual_property_default(
        {"visualProperty": f"NODE_CUSTOMGRAPHICS_SIZE_{slot}", "value": size},
        style_name=style,
    )


def add_custom_png(create_png, style=None, column="fig_location", slot=1, position="below",
                   size=110, unique_dir=None, to_png=False, network=None, **kwargs):
    '''Show an image (e.g. a small plot of experimental data) by each node, from a
    function that gives the image of a node.

    Combines :func:`load_node_images` and :func:`show_node_images`.

    Parameters
    ----------
    create_png : callable
        ``create_png(node_name, **kwargs)`` returns the path of a node's image, or None for
        no image. It can make the image (e.g. plot it with matplotlib), or only return the
        path of an existing one.
    style : str, optional
        Visual style to show the images in (default: the network's current style). The
        bundled SKM styles (*SKM*, *SKM-reactions*) are shared by all networks, so they are
        not changed: a copy, ``"<style>-<column>"``, is made (or replaced) and applied to the
        network instead.
    column : str
        Node table column for the image locations (default ``"fig_location"``).
    slot, position, size
        As for :func:`show_node_images`.
    unique_dir, to_png
        As for :func:`load_node_images`.
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).
    **kwargs
        Passed on to `create_png`.

    Returns
    -------
    str
        The style the images are shown in.
    '''
    images = {}
    for node in p4c.get_all_nodes(network=network):
        path = create_png(node, **kwargs)
        if path:
            images[node] = path

    if style is None:
        style = p4c.styles.get_current_style(network=network)
    if style in set(resources.BUILTIN_STYLES.values()):
        style = copy_style(style, f"{style}-{column}", overwrite=True)
        set_style(style, network)

    load_node_images(images, column, network=network, unique_dir=unique_dir, to_png=to_png)
    show_node_images(style, column, slot=slot, position=position, size=size)
    return style


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

def export_network(filename, format="PDF", wait=1.0, network=None, **kwargs):
    '''Export a network view as an image, fitted to the content and with nothing selected.

    Parameters
    ----------
    filename : str or pathlib.Path
        File to write (overwritten if it exists). py4cytoscape adds the extension of the
        format if the name doesn't end with it.
    format : str
        ``"PDF"`` (default), ``"PNG"``, ``"SVG"``, ``"JPEG"``, ...: see
        ``py4cytoscape.network_views.export_image``.
    wait : float
        Seconds to wait before exporting (default 1). Cytoscape applies style changes (e.g.
        a new style, node images or charts) in the background, and an image exported
        right after them can miss some of them.
    network : int or str, optional
        Cytoscape network, SUID or name (default: the current network).
    **kwargs
        Passed on to ``py4cytoscape.network_views.export_image``, e.g. ``zoom=300`` or
        ``transparent_background=True`` for PNG. Needs Cytoscape 3.10 or later.

    Returns
    -------
    pathlib.Path
        The file written.
    '''
    # fit content ignores edges that may extend
    # past the node boundaries
    # Reported bug: CSD-979
    p4c.network_views.fit_content(network=network)
    p4c.network_selection.clear_selection(type='both', network=network)
    time.sleep(wait)

    result = p4c.network_views.export_image(
        filename=str(Path(filename).absolute()),
        type=format,
        network=network,
        overwrite_file=True,
        **{"all_graphics_details": True, **kwargs},
    )
    return Path(result["file"])


def export_collection(folder, format="PNG", crop=False, network=None, **kwargs):
    '''Export every network of a collection as an image (see :func:`export_network`).

    Parameters
    ----------
    folder : str or pathlib.Path
        Folder for the images (created if needed), named ``<network name>_<SUID>``, with the
        extension of the format.
    format : str
        Image format (default ``"PNG"``).
    crop : bool
        PDF only: crop the margins of each file (see :func:`skm_tools.pdf_utils.crop_pdf`;
        needs the ``pdf`` extra). To combine the files into one PDF, see
        :func:`skm_tools.pdf_utils.combine_pdfs`.
    network : int or str, optional
        Any network of the collection, SUID or name (default: the current network).
    **kwargs
        Passed on to :func:`export_network`.

    Returns
    -------
    list of pathlib.Path
        The files written.
    '''
    if crop and format.upper() != "PDF":
        raise ValueError("crop=True is for format='PDF'.")
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    if network is None:
        network = p4c.get_network_suid()

    files = []
    for suid in sorted(p4c.get_collection_networks(p4c.get_collection_suid(network))):
        name = re.sub(r"[^A-Za-z0-9]+", "_", p4c.get_network_name(suid)).strip("_")
        files.append(export_network(folder / f"{name}_{suid}", format=format, network=suid, **kwargs))
    if crop:
        from .pdf_utils import crop_pdf
        for f in files:
            crop_pdf(f)
    return files
