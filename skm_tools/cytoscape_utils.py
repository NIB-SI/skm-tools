'''Cytoscape automation: load networks, apply the SKM styles, highlight nodes, edges and
paths, create subnetworks, and export images.

Requires the ``cytoscape`` extra (``pip install skm-tools[cytoscape]``) and a running
Cytoscape (https://cytoscape.org) for py4cytoscape to talk to.

py4cytoscape's console output (e.g. the text of errors that are handled, such as the
retries in :func:`clone_network`) is silenced when this module is imported, see
:func:`silence_py4cytoscape`.

Nodes and edges are matched to Cytoscape by name: a node's ``name`` is its networkx node
id, and an edge's ``name`` is ``"source (interaction) target"``, where ``interaction`` is
the edge's ``interaction`` attribute (``"interacts with"`` if it has none), as set by
:func:`load_network` / ``py4cytoscape.create_network_from_networkx``. Parallel edges
between the same two nodes are matched together.
'''

import logging
import re
import shutil
import sys
import time
import uuid
from collections import defaultdict
from pathlib import Path

import networkx as nx
import pandas as pd
import py4cytoscape as p4c

from . import resources
from .paths import path_edges
from .utils import to_node_list


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
    and handled (e.g. the retries in :func:`clone_network`), prints progress messages in
    notebooks, and has a console logger. Silenced when this module is imported. The errors
    are still raised, with the same text, so nothing is lost; py4cytoscape's detailed log
    file (``logs/py4cytoscape.log`` in the working directory) is not affected.

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

def _node_suids(nodes, network=None):
    '''Cytoscape SUIDs of the nodes named `nodes` (names not in the network are skipped).'''
    nodes = {str(n) for n in nodes}
    table = p4c.tables.get_table_columns(table="node", columns=["name"], network=network)
    return [int(suid) for suid, name in table["name"].items() if name in nodes]


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
    '''Cytoscape SUIDs of the edges between each (u, v) pair (all parallel edges).'''
    edge_pairs = list(edge_pairs)
    if not edge_pairs:
        return []
    table = p4c.tables.get_table_columns(table="edge", columns=["name"], network=network)
    return [int(suid) for suid in _match_edge_names(table["name"].to_dict(), edge_pairs)]


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
        Cytoscape network (default: the current network).
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
        apply_builtin_style(suid, style)
    return suid


def apply_builtin_style(suid, style="skm"):
    '''Apply one of the bundled SKM visual styles to a Cytoscape network.

    The styles are imported into Cytoscape on first use (from ``skm-styles.xml``).

    Parameters
    ----------
    suid : int
        Cytoscape network SUID.
    style : {"skm", "skm-reactions"}
        ``"skm"`` (style *SKM*, default): for CKN, the PSS interaction network and the PSS
        gene networks. Nodes by ``node_type`` (the PSS classes, in the PSS Explorer's
        colours), labelled with ``display_label``; edges coloured and with arrows by
        ``interaction``, dashed for mutual (``directed`` False) edges, and thicker for better
        supported ones (``rank``).
        ``"skm-reactions"`` (style *SKM-reactions*): for the PSS reaction graph, with
        reactions as nodes, and edges by ``edge_type`` (activation, inhibition, substrate,
        product, ...) labelled with their ``role``.
        The older names ``"pss"`` and ``"ckn"`` apply *SKM*.
    '''
    key = style.lower()
    if key not in resources.BUILTIN_STYLES:
        raise ValueError(f"apply_builtin_style expects one of {sorted(resources.BUILTIN_STYLES)}, "
                         f"not {style!r}.")
    style_name = resources.BUILTIN_STYLES[key]

    if style_name not in p4c.styles.get_visual_style_names():
        p4c.import_visual_styles(resources.get_style_xml_path())

    set_style(style_name, suid)
    print(f"Applied {style_name} to {suid}")


# ---------------------------------------------------------------------------
# Highlighting (style bypasses)
# ---------------------------------------------------------------------------

def highlight_nodes(node_names, colour=None, label_color=None, border_color=None, border_width=None,
                    node_height=None, node_width=None, network=None):
    '''Highlight nodes with style bypasses.

    Only the given properties are changed. Bypasses stay when the visual style changes;
    clear them in Cytoscape (or with py4cytoscape) to undo.

    Parameters
    ----------
    node_names : node or iterable of nodes
        Nodes to highlight (names not in the network are skipped).
    colour, label_color, border_color : str, optional
        Fill, label and border colours, as hex (``"#FF0000"``).
    border_width, node_height, node_width : float, optional
        Sizes in pixels.
    network : int or str, optional
        Cytoscape network (default: the current network).

    Returns
    -------
    list of int
        SUIDs of the highlighted nodes.
    '''
    nodes_by_suid = _node_suids(to_node_list(node_names), network=network)
    if not nodes_by_suid:
        return nodes_by_suid

    # to set back, see bug https://github.com/cytoscape/py4cytoscape/issues/114
    og_style = p4c.styles.get_current_style(network)

    bypasses = [
        (colour, p4c.style_bypasses.set_node_color_bypass),
        (label_color, p4c.style_bypasses.set_node_label_color_bypass),
        (border_color, p4c.style_bypasses.set_node_border_color_bypass),
        (border_width, p4c.style_bypasses.set_node_border_width_bypass),
        (node_height, p4c.style_bypasses.set_node_height_bypass),
        (node_width, p4c.style_bypasses.set_node_width_bypass),
    ]
    for value, set_bypass in bypasses:
        if value:
            set_bypass(nodes_by_suid, value, network=network)

    set_style(og_style, network)

    return nodes_by_suid


def highlight_edges(edge_pairs, colour, skip_edges=None, edge_line_width=10, network=None):
    '''Highlight edges with style bypasses (colour and line width).

    Parameters
    ----------
    edge_pairs : iterable of tuple
        Edges as (u, v) or (u, v, key) tuples, e.g. from :func:`skm_tools.paths.path_edges`.
        All parallel edges between u and v are highlighted.
    colour : str
        Edge colour, as hex.
    skip_edges : iterable of tuple, optional
        (u, v) pairs not to highlight, e.g. edges already highlighted in another colour.
    edge_line_width : float
        Line width in pixels (default 10).
    network : int or str, optional
        Cytoscape network (default: the current network).

    Returns
    -------
    list of tuple
        The (u, v) pairs highlighted (after skipping).
    '''
    skip = {(e[0], e[1]) for e in (skip_edges or [])}
    edges = list(dict.fromkeys((e[0], e[1]) for e in edge_pairs if (e[0], e[1]) not in skip))

    edges_by_suid = _edge_suids(edges, network=network)
    if edges_by_suid:
        p4c.style_bypasses.set_edge_line_width_bypass(edges_by_suid, edge_line_width, network=network)
        p4c.style_bypasses.set_edge_color_bypass(edges_by_suid, colour, network=network)

    return edges


def highlight_path(node_names, colour, skip_nodes=None, skip_edges=None, label_color="white",
                   border_color="black", border_width=10, edge_line_width=10, network=None):
    '''Highlight the nodes and edges of a path.

    To highlight several paths in different colours without overwriting, pass the nodes
    and edges returned by earlier calls as `skip_nodes` and `skip_edges`.

    Parameters
    ----------
    node_names : list
        The path, as a list of nodes.
    colour : str
        Node fill and edge colour, as hex.
    skip_nodes : iterable, optional
        Nodes not to highlight.
    skip_edges : iterable of tuple, optional
        (u, v) pairs not to highlight.
    label_color, border_color : str
        Node label and border colours (default white and black).
    border_width, edge_line_width : float
        Node border and edge widths in pixels (default 10).
    network : int or str, optional
        Cytoscape network (default: the current network).

    Returns
    -------
    nodes : list
        The path's nodes.
    edges : list of tuple
        The path's (u, v) edges that were highlighted.
    '''
    node_names = list(node_names)
    skip_nodes = set(skip_nodes or [])
    nodes_for_highlight = [n for n in node_names if n not in skip_nodes]

    # the edges are highlighted even if all nodes already are (e.g. a path between two coloured nodes)
    if nodes_for_highlight:
        highlight_nodes(
            nodes_for_highlight,
            colour=colour,
            label_color=label_color,
            border_color=border_color,
            border_width=border_width,
            network=network
        )

    edge_pairs = list(zip(node_names, node_names[1:]))
    edges = highlight_edges(edge_pairs, colour, skip_edges=skip_edges, edge_line_width=edge_line_width,
                            network=network)

    return node_names, edges


def contrast_colour(colour):
    '''The complementary colour, e.g. for a label on a node of colour `colour`.

    Parameters
    ----------
    colour : str
        Hex colour (``"#RRGGBB"``).

    Returns
    -------
    str
        Hex colour.
    '''
    rgb = int(colour.lstrip('#'), 16)
    complementary_colour = 0xffffff-rgb
    return f'#{complementary_colour:06X}'


def apply_shortest_paths_style(sources, path_lists, target, g, edge_colors=None, node_colors=None, network=None):
    '''Colour shortest-path results from several searches into a new visual style.

    Adds node columns ``distance-to-target`` and ``node-path-source``, and the edge column
    ``edge-priority`` (``"direct path (<source>)"``), then copies the current style to
    ``<style>-shortest-paths-query`` with colour mappings on them.

    Parameters
    ----------
    sources : list
        Source node of each search.
    path_lists : list of list
        Paths found by each search: one list of paths per source, in the same order.
    target : node
        The common target of the searches.
    g : networkx.Graph
        Graph the paths were found in (for distances to the target).
    edge_colors : list of str, optional
        Hex colours for the edges of each source's paths, one per source.
    node_colors : list of str, optional
        Hex colours for a continuous node colour mapping on distance to the target
        (closest first); at least two.
    network : int or str, optional
        Cytoscape network (default: the current network).
    '''
    og_style = p4c.styles.get_current_style(network)
    new_style = f'{og_style}-shortest-paths-query'
    p4c.styles.copy_visual_style(og_style, new_style)

    path_searches_attributes = defaultdict(dict)
    edge_priority = {}

    for source, paths in zip(sources, path_lists):
        nodes = {n for p in paths for n in p}
        for node in nodes:
            path_searches_attributes[node]['distance-to-target'] = nx.shortest_path_length(g, source=node, target=target)
            # the first search that reaches a node claims it
            path_searches_attributes[node].setdefault('node-path-source', source)

        for suid in _edge_suids(path_edges(paths, g), network=network):
            edge_priority.setdefault(suid, {'edge-priority': f"direct path ({source})"})

    p4c.tables.load_table_data(
        pd.DataFrame.from_dict(path_searches_attributes, orient='index').rename_axis('name').reset_index(),
        data_key_column='name',
        table='node',
        table_key_column="name",
        network=network
    )

    p4c.tables.load_table_data(
        pd.DataFrame.from_dict(edge_priority, orient='index').rename_axis('SUID').reset_index(),
        data_key_column='SUID',
        table='edge',
        table_key_column="SUID",
        network=network
    )

    if node_colors:
        max_len = max(x['distance-to-target'] for x in path_searches_attributes.values())
        n = len(node_colors)
        node_color_mapping_range = [max_len * i / (n - 1) for i in range(n)]

        p4c.style_mappings.set_node_color_mapping(
            'distance-to-target',
            table_column_values=node_color_mapping_range,
            colors=node_colors,
            mapping_type='c',
            style_name=new_style,
            network=network
        )

        p4c.style_mappings.set_node_label_color_mapping(
            'distance-to-target',
            table_column_values=node_color_mapping_range,
            colors=[contrast_colour(x) for x in node_colors],
            mapping_type='c',
            style_name=new_style,
            network=network
        )

    if edge_colors:
        p4c.style_mappings.set_edge_color_mapping(
            'edge-priority',
            table_column_values=[f'direct path ({source})' for source in sources],
            colors=edge_colors,
            mapping_type='d',
            style_name=new_style,
            network=network
        )

    set_style(new_style, network)


# ---------------------------------------------------------------------------
# Networks, collections and styles
# ---------------------------------------------------------------------------

_RETRIES = 20
_RETRY_WAIT = 0.5  # seconds

def clone_network(network, name=None, collection=None):
    '''Clone a Cytoscape network into a new collection, optionally renaming both.

    py4cytoscape's ``clone_network`` doesn't rename the new collection; this does.

    Parameters
    ----------
    network : int or str
        Cytoscape network to clone.
    name : str, optional
        Name of the clone.
    collection : str, optional
        Name of the clone's new collection.

    Returns
    -------
    int
        SUID of the clone.
    '''
    suid = p4c.networks.clone_network(network=network)

    if collection:
        collection_suid = p4c.collections.get_collection_suid(suid)
        p4c.commands.cyrest_put(
            f"collections/{collection_suid}/tables/default",
            body={"key": "SUID", "data": [{"SUID": collection_suid, "name": collection}]},
            require_json=False,
        )

    if name:
        # right after cloning, Cytoscape may not know the new network yet ("unrecognized
        # table entry"), so retry for a few seconds
        for attempt in range(_RETRIES):
            try:
                p4c.rename_network(name, network=suid)
                break
            except p4c.CyError:
                if attempt == _RETRIES - 1:
                    raise
                time.sleep(_RETRY_WAIT)

    return suid


def copy_style(style, new_style, networks=()):
    '''Copy a visual style, and apply the copy to `networks`.

    Parameters
    ----------
    style : str
        Visual style to copy.
    new_style : str
        Name of the copy.
    networks : iterable of int or str, optional
        Networks to apply the copy to.

    Returns
    -------
    str
        `new_style`.
    '''
    p4c.copy_visual_style(style, new_style)
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

def subnetwork_edge_induced(edge_pairs, parent_suid, name="subnetwork (edge induced)"):
    '''New Cytoscape network with only the given edges (and their nodes).

    Parameters
    ----------
    edge_pairs : iterable of tuple
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
    edge_pairs = list(edge_pairs)
    nodes = {n for e in edge_pairs for n in e[:2]}

    return p4c.networks.create_subnetwork(
        nodes=_node_suids(nodes, network=parent_suid),
        edges=_edge_suids(edge_pairs, network=parent_suid),
        subnetwork_name=name,
        network=parent_suid,
        exclude_edges=True,
    )


def subnetwork_edge_induced_from_paths(paths, g, parent_suid, name="subnetwork (edge induced)"):
    '''New Cytoscape network with only the edges along `paths`.

    Parameters
    ----------
    paths : list of list
        Paths as lists of nodes, e.g. from :func:`skm_tools.paths.get_paths`.
    g : networkx.Graph
        Graph the paths were found in (to find the edges, see :func:`skm_tools.paths.path_edges`).
    parent_suid : int
        SUID of the Cytoscape network to take them from.
    name : str
        Name of the new network.

    Returns
    -------
    int
        SUID of the new network.
    '''
    return subnetwork_edge_induced(path_edges(paths, g), parent_suid, name=name)


def subnetwork_node_induced(nodes, parent_suid, name="subnetwork (node induced)"):
    '''New Cytoscape network with the given nodes and every edge between them.

    Node-induced: the result includes every edge of `parent_suid` between the given
    nodes, not just edges from whatever selection process produced `nodes` (e.g. a
    curated reaction/edge list). To keep an exact edge subset instead, use
    :func:`subnetwork_edge_induced`.

    Parameters
    ----------
    nodes : iterable
        Nodes to include.
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
        nodes=_node_suids(to_node_list(nodes), network=parent_suid),
        subnetwork_name=name,
        network=parent_suid
    )


def get_or_create_subnetwork(nodes, parent_suid, name):
    '''The network called `name`, created with :func:`subnetwork_node_induced` if it
    doesn't exist yet.

    For re-running notebook cells without making duplicate networks. Network names are
    matched across all collections.

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
    if name in p4c.get_network_list():
        return p4c.get_network_suid(name)
    return subnetwork_node_induced(nodes, parent_suid, name=name)


def subnetwork_neighbours(nodes, parent_suid, name="subnetwork (1st neighbours)"):
    '''New Cytoscape network with `nodes`, their first neighbours, and every edge between them.

    For deeper or directed neighbourhoods, use :func:`skm_tools.neighbors.neighborhood_nodes`
    and pass the result to :func:`subnetwork_node_induced`.

    Parameters
    ----------
    nodes : iterable
        Nodes whose neighbourhood to include.
    parent_suid : int
        SUID of the Cytoscape network to take them from.
    name : str
        Name of the new network.

    Returns
    -------
    int
        SUID of the new network.
    '''
    p4c.select_nodes(
        _node_suids(to_node_list(nodes), network=parent_suid),
        by_col="SUID",
        preserve_current_selection=False,
        network=parent_suid
    )
    neighbours = p4c.select_first_neighbors(network=parent_suid)

    return p4c.networks.create_subnetwork(
        nodes=neighbours['nodes'],
        subnetwork_name=name,
        network=parent_suid
    )


# ---------------------------------------------------------------------------
# Layout, custom graphics and export
# ---------------------------------------------------------------------------

def layout_from_coords(network, table):
    '''Place nodes at given coordinates (e.g. from a graphviz or networkx layout).

    Parameters
    ----------
    network : int
        Cytoscape network SUID.
    table : pandas.DataFrame
        Indexed by node name, with columns ``x`` and ``y``.
    '''
    p4c.load_table_data(table, network=network)

    current_style = p4c.styles.get_current_style(network=network)

    tmp_style = 'tmp-layout'
    p4c.copy_visual_style(current_style, tmp_style)
    set_style(tmp_style, network)

    p4c.update_style_mapping(tmp_style, p4c.map_visual_property('NODE_X_LOCATION', 'x', 'p'))
    p4c.update_style_mapping(tmp_style, p4c.map_visual_property('NODE_Y_LOCATION', 'y', 'p'))

    p4c.delete_visual_style(tmp_style)

    p4c.tables.delete_table_column('x', network=network)
    p4c.tables.delete_table_column('y', network=network)

    set_style(current_style, network)


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


def node_file_key(name):
    '''Default key for matching file names to node names (see :func:`match_files_to_nodes`).

    Replaces the characters that can't (or shouldn't) be in file names with ``_``, e.g.
    ``"WRKY33[fc00166]"`` -> ``"WRKY33_fc00166_"`` and ``"11-/12-OH-JA"`` ->
    ``"11-_12-OH-JA"``; gene ids stay as they are.

    Parameters
    ----------
    name : str
        Node name, or file name without extension.

    Returns
    -------
    str
    '''
    return re.sub(r"[^A-Za-z0-9_-]", "_", str(name))


def match_files_to_nodes(folder, nodes, key=node_file_key, aliases=None, extensions=("png", "svg")):
    '''Match image files (e.g. one plot per gene or metabolite) to nodes by name.

    Parameters
    ----------
    folder : str or pathlib.Path
        Folder with the files.
    nodes : iterable
        Node names, e.g. ``g.nodes()`` or ``py4cytoscape.get_all_nodes()``.
    key : callable
        Applied to node names and file names (without extension); they match if the keys
        are equal. Default :func:`node_file_key`.
    aliases : dict, optional
        File name (without extension) -> node name, for files not named after their node
        (e.g. ``{"Pro": "Proline accumulation"}``).
    extensions : iterable of str
        File extensions to use (default png and svg).

    Returns
    -------
    dict
        Node -> file path. Nodes with the same key get the same file. Files that match no
        node are printed.
    '''
    aliases = aliases or {}
    nodes_by_key = defaultdict(list)
    for n in nodes:
        nodes_by_key[key(n)].append(n)

    matched, unmatched = {}, []
    for ext in extensions:
        for path in sorted(Path(folder).glob(f"*.{ext}")):
            k = key(aliases.get(path.stem, path.stem))
            if k in nodes_by_key:
                matched.update({n: path for n in nodes_by_key[k]})
            else:
                unmatched.append(path.name)

    if unmatched:
        print(f"{len(unmatched)} files in {folder} match no node: {unmatched}")
    return matched


def unique_image_copies(images, folder, to_png=False):
    '''Copy images to new, unique file names, optionally converting SVG to PNG.

    Cytoscape caches images, and doesn't always notice that a file is a different one (e.g.
    the same file name in another folder) or has changed: it may show an old image instead.
    Every call makes new copies, with names that no other file has had (``<name>_<random
    id>``), so Cytoscape has to load them again. Old copies are not deleted.

    Parameters
    ----------
    images : dict
        Node -> image path.
    folder : str or pathlib.Path
        Folder for the copies (created if needed).
    to_png : bool
        Convert SVG images to PNG. Needs ``cairosvg`` (``pip install cairosvg``).

    Returns
    -------
    dict
        Node -> path of the copy. Nodes with the same image share one copy.
    '''
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)

    copy_of = {}
    for path in {Path(p) for p in images.values()}:
        convert = to_png and path.suffix.lower() == ".svg"
        copy = folder / f"{path.stem}_{uuid.uuid4().hex[:12]}{'.png' if convert else path.suffix}"
        if convert:
            import cairosvg
            cairosvg.svg2png(url=str(path), write_to=str(copy))
        else:
            shutil.copyfile(path, copy)
        copy_of[path] = copy
    return {node: copy_of[Path(p)] for node, p in images.items()}


def load_node_images(images, column, network=None, unique_dir=None, to_png=False):
    '''Load image file locations into a node table column, for :func:`show_node_images`.

    Parameters
    ----------
    images : dict
        Node name -> image path (e.g. from :func:`match_files_to_nodes`). Nodes not in the
        network are ignored.
    column : str
        Node table column to load them into, e.g. ``"image_heat_Desiree"``.
    network : int or str, optional
        Cytoscape network (default: the current network).
    unique_dir : str or pathlib.Path, optional
        Copy the images to this folder first, each to a new, unique file name (see
        :func:`unique_image_copies`), so that Cytoscape's image cache can't show an old or
        a different image instead.
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


def add_custom_png(network, create_png, style=None, column="fig_location", **kwargs):
    '''Show an image (e.g. a small plot of experimental data) below each node, from a
    function that gives the image of a node.

    Combines :func:`load_node_images` and :func:`show_node_images`.

    Parameters
    ----------
    network : int
        Cytoscape network SUID.
    create_png : callable
        ``create_png(node_name, **kwargs)`` returns the path of a node's image, or None for
        no image. It can make the image (e.g. plot it with matplotlib), or only return the
        path of an existing one.
    style : str, optional
        Visual style to show the images in (default: the network's current style).
    column : str
        Node table column for the image locations (default ``"fig_location"``).
    **kwargs
        Passed on to `create_png`.
    '''
    images = {}
    for node in p4c.get_all_nodes(network=network):
        path = create_png(node, **kwargs)
        if path:
            images[node] = path

    if style is None:
        style = p4c.styles.get_current_style(network=network)

    load_node_images(images, column, network=network)
    show_node_images(style, column)


def chart_column(df, columns, colours, chart="barchart", value_range=None, labels=None,
                 na_value=0, **options):
    '''Chart definitions for each row of `df`, to show as node charts in Cytoscape.

    Builds `enhancedGraphics <https://apps.cytoscape.org/apps/enhancedgraphics>`_ chart
    strings (e.g. ``barchart: colorlist="..." valuelist="..."``); load them as a node table
    column and show them like images, with :func:`show_node_images`. Needs the
    enhancedGraphics app in Cytoscape.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per node.
    columns : list of str
        Columns of `df` to chart, in order.
    colours : str or list of str
        Hex colour per column, or one for all.
    chart : str
        enhancedGraphics chart type, e.g. ``"barchart"`` (default), ``"linechart"``,
        ``"heatstripchart"``, ``"piechart"``.
    value_range : tuple of float, optional
        (min, max) of the value axis (default: per chart, from its values).
    labels : list of str, optional
        Label per column (default: the column names).
    na_value : float
        Value for missing data (default 0).
    **options
        Further enhancedGraphics options, e.g. ``separation=2``, ``ybase=1``.

    Returns
    -------
    pandas.Series
        Chart string per row of `df`.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({"t1": [1.5], "t2": [-0.5]}, index=["AT2G38470"])
    >>> chart_column(df, ["t1", "t2"], "#E41A1C", value_range=(-2, 2)).iloc[0]
    'barchart: colorlist="#E41A1C,#E41A1C" valuelist="1.5,-0.5" labellist="t1,t2" range="-2,2"'
    '''
    if isinstance(colours, str):
        colours = [colours] * len(columns)
    labels = labels or columns

    fixed = f'colorlist="{",".join(colours)}"'
    rest = f'labellist="{",".join(map(str, labels))}"'
    if value_range is not None:
        rest += f' range="{value_range[0]:g},{value_range[1]:g}"'
    rest += "".join(f" {k}={v}" for k, v in options.items())

    values = df[columns].astype(float).fillna(na_value)
    return values.apply(
        lambda row: f'{chart}: {fixed} valuelist="{",".join(f"{v:g}" for v in row)}" {rest}',
        axis=1,
    )


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

def export_network(network, filename, format="PDF", **kwargs):
    '''Export a network view as an image, fitted to the content and with nothing selected.

    Parameters
    ----------
    network : int
        Cytoscape network SUID.
    filename : str or pathlib.Path
        File to write (overwritten if it exists).
    format : str
        ``"PDF"`` (default), ``"PNG"``, ``"SVG"``, ...: see
        ``py4cytoscape.network_views.export_image``.
    **kwargs
        Passed on to ``py4cytoscape.network_views.export_image``, e.g. ``zoom=300`` or
        ``transparent_background=True`` for PNG.
    '''
    # fit content ignores edges that may extend
    # past the node boundaries
    # Reported bug: CSD-979
    p4c.network_views.fit_content(network=network)
    p4c.network_selection.clear_selection(type='both', network=network)

    p4c.network_views.export_image(
        filename=str(Path(filename).absolute()),
        type=format,
        network=network,
        overwrite_file=True,
        **{"all_graphics_details": True, **kwargs},
    )


def export_collection(network, folder, format="PNG", **kwargs):
    '''Export every network of a collection as an image (see :func:`export_network`).

    Parameters
    ----------
    network : int
        Any network of the collection.
    folder : str or pathlib.Path
        Folder for the images (created if needed), named ``<network name>_<SUID>``.
    format : str
        Image format (default ``"PNG"``).
    **kwargs
        Passed on to :func:`export_network`.

    Returns
    -------
    list of pathlib.Path
        The files written.
    '''
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)

    files = []
    for suid in sorted(p4c.get_collection_networks(p4c.get_collection_suid(network))):
        name = re.sub(r"[^A-Za-z0-9]+", "_", p4c.get_network_name(suid)).strip("_")
        filename = folder / f"{name}_{suid}.{format.lower()}"
        export_network(suid, filename, format=format, **kwargs)
        files.append(filename)
    return files
