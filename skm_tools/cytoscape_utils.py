'''Cytoscape automation: load networks, apply the SKM styles, highlight nodes, edges and
paths, create subnetworks, and export images.

Requires the ``cytoscape`` extra (``pip install skm-tools[cytoscape]``) and a running
Cytoscape (https://cytoscape.org) for py4cytoscape to talk to.

Nodes and edges are matched to Cytoscape by name: a node's ``name`` is its networkx node
id, and an edge's ``name`` is ``"source (interaction) target"``, where ``interaction`` is
the edge's ``interaction`` attribute (``"interacts with"`` if it has none), as set by
:func:`load_network` / ``py4cytoscape.create_network_from_networkx``. Parallel edges
between the same two nodes are matched together.
'''

from collections import defaultdict
from pathlib import Path

import networkx as nx
import pandas as pd
import py4cytoscape as p4c

from . import resources
from .paths import path_edges
from .utils import to_node_list


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
    style : {"pss", "ckn"}, optional
        Built-in SKM style to apply (see :func:`apply_builtin_style`).

    Returns
    -------
    int
        SUID of the new Cytoscape network.
    '''
    suid = p4c.networks.create_network_from_networkx(g, title=title, collection=collection)
    if style:
        apply_builtin_style(suid, style)
    return suid


def apply_builtin_style(suid, style):
    '''Apply one of the bundled SKM visual styles to a Cytoscape network.

    The styles are imported into Cytoscape on first use.

    Parameters
    ----------
    suid : int
        Cytoscape network SUID.
    style : {"pss", "ckn"}
        ``"pss"`` (PSS-default: shapes by ``node_type``, fill by ``pathway``, arrows by
        ``interaction`` influence, dashed for mutual (``directed`` False) edges) or
        ``"ckn"`` (CKN-default).
    '''
    style = style.lower()

    if not (style in resources.BUILTIN_STYLES):
        raise ValueError(f"apply_builtin_style expects a value in {resources.BUILTIN_STYLES}.")

    style_name = {
        'ckn':resources.CKN_DEFAULT_STYLE,
        'pss':resources.PSS_DEFAULT_STYLE
    }[style]

    if not style_name in p4c.styles.get_visual_style_names():
        p4c.import_visual_styles(resources.get_style_xml_path())

    p4c.styles.set_visual_style(style_name, network=suid)
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

    p4c.styles.set_visual_style(og_style, network=network)

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

    p4c.styles.set_visual_style(new_style, network=network)


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
    p4c.set_visual_style(tmp_style, network=network)

    p4c.update_style_mapping(tmp_style, p4c.map_visual_property('NODE_X_LOCATION', 'x', 'p'))
    p4c.update_style_mapping(tmp_style, p4c.map_visual_property('NODE_Y_LOCATION', 'y', 'p'))

    p4c.delete_visual_style(tmp_style)

    p4c.tables.delete_table_column('x', network=network)
    p4c.tables.delete_table_column('y', network=network)

    p4c.set_visual_style(current_style, network=network)


def add_custom_png(network, create_png, style=None, **kwargs):
    '''Show an image (e.g. a small plot of experimental data) under each node.

    Parameters
    ----------
    network : int
        Cytoscape network SUID.
    create_png : callable
        ``create_png(node_name, **kwargs)`` makes an image for a node and returns its
        path (pathlib.Path), or None for no image.
    style : str, optional
        Visual style to add the image mapping to (default: the network's current style).
    **kwargs
        Passed on to `create_png`.
    '''
    nodes = p4c.get_all_nodes(network=network)

    node_pngs = {}
    for node in nodes:
        node_png_fname = create_png(node, **kwargs)
        if node_png_fname:
            node_pngs[node] = {'fig_location':f"file:{str(Path(node_png_fname).absolute())}"}

    table = pd.DataFrame.from_dict(node_pngs, orient='index')

    if style is None:
        style = p4c.styles.get_current_style(network=network)

    p4c.load_table_data(table, network=network)
    p4c.style_dependencies.sync_node_custom_graphics_size(False, style_name=style)

    style_mapping = p4c.style_mappings.map_visual_property(
        visual_prop="NODE_CUSTOMGRAPHICS_1",
        table_column="fig_location",
        mapping_type="p",
    )
    p4c.style_mappings.update_style_mapping(style, style_mapping)

    p4c.style_defaults.set_visual_property_default({
        'visualProperty':'NODE_CUSTOMGRAPHICS_POSITION_1',
        'value':'S,N,c,0.00,10.00'},
        style_name=style
    )

    p4c.style_defaults.set_visual_property_default({
        'visualProperty':'NODE_CUSTOMGRAPHICS_SIZE_1',
        'value':110},
        style_name=style
    )


def export_network(network, filename, format="PDF"):
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
        all_graphics_details=True,
    )
