'''Functions to import, extract subnetworks and style in Cytoscape'''

from collections.abc import Sequence
from collections import defaultdict
import pandas as pd
import networkx as nx
import py4cytoscape as p4c
from . import resources

import io, re

def _cytoscape_safe_names(names):
    '''Escape names for Cystoscape '''
    return [x.replace(",", r"\,") for x in list(set(names))]

def apply_builtin_style(suid, style):
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

def highlight_nodes(node_names, colour=None, label_color=None, border_color=None, border_width=None, node_height=None, node_width=None, network=None):
    '''
    Highlight nodes in a Cytoscape networks using style bypasses.
    '''

    # to set back, see bug https://github.com/cytoscape/py4cytoscape/issues/114
    og_style = p4c.styles.get_current_style(network)

    escaped_names = _cytoscape_safe_names(node_names)
    nodes_by_suid = p4c.select_nodes(escaped_names, by_col="name", network=network)['nodes']

    if colour:
        p4c.style_bypasses.set_node_color_bypass(
            node_names=nodes_by_suid,
            new_colors=colour,
            network=network
        )

    if label_color:
        p4c.style_bypasses.set_node_label_color_bypass(
            nodes_by_suid,
            label_color,
            network=network
        )

    if border_color:
        p4c.style_bypasses.set_node_border_color_bypass(
            nodes_by_suid,
            border_color,
            network=network
        )

    if border_width:
        p4c.style_bypasses.set_node_border_width_bypass(
            nodes_by_suid,
            border_width,
            network=network
        )

    if node_height:
        p4c.style_bypasses.set_node_height_bypass(
            nodes_by_suid,
            node_height,
            network=network
        )

    if node_width:
        p4c.style_bypasses.set_node_width_bypass(
            nodes_by_suid,
            node_width,
            network=network
        )

    p4c.styles.set_visual_style(og_style)

def highlight_path(node_names, colour, skip_nodes=None, skip_edges=None, label_color="white", border_color="black", border_width=10, edge_line_width=10, network=None):
    '''Bypasses
    '''

    if isinstance(skip_nodes, Sequence) and not isinstance(skip_nodes, str):
        node_names = [n for n in node_names if not n in skip_nodes]

    if len(node_names) == 0:
        print("No more nodes to colour")
        return [], []

    if not isinstance(skip_edges, Sequence) and not isinstance(skip_edges, str):
        skip_edges = []

    highlight_nodes(
        node_names,
        colour=colour,
        label_color=label_color,
        border_color=border_color,
        border_width=border_width,
        network=network
    )

    edge_pairs = []
    for s, t in zip(node_names, node_names[1:]):
        edge_pairs.append((s, t))
    edges = highlight_edges(edge_pairs, colour, skip_edges=skip_edges, edge_line_width=edge_line_width, network=network)

    return node_names, edges

def highlight_edges(edge_pairs, colour, skip_edges=None, edge_line_width=10, network=None):

    if skip_edges is None:
        skip_edges = []

    edges = []
    for s, t in edge_pairs:
        e = f"{s} (interacts with) {t}"
        if not (e in skip_edges):
            edges.append(e)

    escaped_names = _cytoscape_safe_names(set(edges))
    edges_by_suid = p4c.select_edges(escaped_names, by_col="name", network=network)['edges']

    p4c.style_bypasses.set_edge_line_width_bypass(edges_by_suid, edge_line_width, network=network)
    p4c.style_bypasses.set_edge_color_bypass(edges_by_suid, colour, network=network)

    return edges

def get_path_edges(paths, g):
    '''path [a, b, c] --> edges along path [(a, b), (b, c)]
    '''
    edges = []

    if g.is_multigraph():
        for path in paths:
            for s, t in zip(path, path[1:]):
                edges_here = [(u, v, k) for (u, v, k) in g.edges(keys=True) if (u==s) and (v==t)]
                edges += edges_here
        cytoscape_edges = _cytoscape_safe_names([f"{u} (interacts with) {v}" for u, v, k in edges])

    else:
        for path in paths:
            for s, t in zip(path, path[1:]):
                edges_here = [(u, v) for (u, v) in g.edges() if (u==s) and (v==t)]
                edges += edges_here
        cytoscape_edges = _cytoscape_safe_names([f"{u} (interacts with) {v}" for u, v in edges])

    return edges, cytoscape_edges

def subnetwork_edge_induced_from_paths(paths, g, parent_suid, name="subnetwork (edge induced)"):

    '''Subnetwork from existing Cytosape network.

    Notes
    -----

    First we find the Cytoscape identifiers for the nodes we are interested in
    by "selecting" the nodes in Cytoscape. Due to a py4cytoscape bug, we need
    to escape the commas in the node names.
    '''

    all_node_names = {y for x in paths for y in x}
    escaped_names = _cytoscape_safe_names(set(all_node_names))
    select_result = p4c.select_nodes(
        escaped_names,
        by_col="name",
        preserve_current_selection=False,
        network=parent_suid
    )

    all_path_nodes = select_result["nodes"]

    _, cytoscape_edges = get_path_edges(paths, g)

    network_suid = p4c.networks.create_subnetwork(
        nodes=all_path_nodes,
        edges=cytoscape_edges,
        subnetwork_name=name,
        network=parent_suid,
        exclude_edges=True,
        edges_by_col="shared name"
    )

    return network_suid

def subnetwork_node_induced(nodes, parent_suid, name="subnetwork (node induced)"):
    '''Subnetwork from existing Cystoscape network.
    '''

    escaped_names = _cytoscape_safe_names(set(nodes))
    select_result = p4c.select_nodes(
        escaped_names,
        by_col="name",
        preserve_current_selection=False,
        network=parent_suid
    )

    network_suid = p4c.networks.create_subnetwork(
        nodes=select_result['nodes'],
        subnetwork_name=name,
        network=parent_suid
    )

    return network_suid

def subnetwork_neighbours(nodes, parent_suid, name="subnetwork (1st neighbours)"):
    ''' Will be node-induced
    '''

    escaped_names = _cytoscape_safe_names(set(nodes))
    p4c.select_nodes(
        escaped_names,
        by_col="name",
        preserve_current_selection=False,
        network=parent_suid
    )

    neighbours = p4c.select_first_neighbors(network=parent_suid)

    network_suid = p4c.networks.create_subnetwork(
        nodes=neighbours['nodes'],
        subnetwork_name=name,
        network=parent_suid
    )

    return network_suid

def contrast_colour(colour):
    '''To contrast lable to node colours'''
    rgb = int(colour.lstrip('#'), 16)
    complementary_colour = 0xffffff-rgb
    return f'#{complementary_colour:06X}'

def apply_shortest_paths_style(sources, path_lists, target, g, edge_colors=None, node_colors=None, network=None ):
    '''Generate a style showing shortest paths, from multiple queries

    Parameters
    ----------

    sources : str
        List of names of sources of the shorstest paths
    path_lists : list of lists
        List of lists of paths, same number of first level lists as in `sources`
    target: str
        Name of target node in shortest paths results
    g: nx.Graph
        Graph in which searches were performed
    edge_colors: list
        Hex values for colouring paths based on sources of the search, same number as sources
    node_colors: list
        Hex values for colouring node based on distance from source in g, any number of colours
    network : int
        Cytoscape suid, where to apply the style
    '''

    og_style = p4c.styles.get_current_style(network)
    new_style = f'{og_style}-shortest-paths-query'
    p4c.styles.copy_visual_style(og_style, new_style)

    path_searches_attributes = defaultdict(dict)
    edge_priority = {}

    skip_source = set()

    for source, paths in zip(sources, path_lists):
        nodes = {n for p in paths for n in p}
        for node in nodes:
            d = len(nx.shortest_path(g, source=node, target=target)) -1
            path_searches_attributes[node]['distance-to-target'] = d
            if not (node in skip_source):
                path_searches_attributes[node]['node-path-source'] = source
                skip_source.update([node])

        _, cytoscape_edges = get_path_edges(paths, g)
        for e in cytoscape_edges:
            if not e in edge_priority:
                edge_priority[e] = {'edge-priority': f"direct path ({source})"}

    p4c.tables.load_table_data(
        pd.DataFrame.from_dict(path_searches_attributes, orient='index').reset_index(),
        data_key_column='index',
        table='node',
        table_key_column="shared name",
        network=network
    )

    p4c.tables.load_table_data(
        pd.DataFrame.from_dict(edge_priority, orient='index').reset_index(),
        data_key_column='index',
        table='edge',
        table_key_column="name",
        network=network
    )

    if node_colors:
        max_len = max([path_searches_attributes[x]['distance-to-target'] \
                       for x in path_searches_attributes])
        node_color_mapping_range = range(0, max_len+1, int(max_len / (len(node_colors)-1)))

        p4c.style_mappings.set_node_color_mapping(
            'distance-to-target',
            table_column_values = node_color_mapping_range,
            colors=node_colors,
            mapping_type='c',
            style_name=new_style,
            network=network
        )

        p4c.style_mappings.set_node_label_color_mapping(
            'distance-to-target',
            table_column_values = node_color_mapping_range,
            colors=[contrast_colour(x) for x in node_colors],
            mapping_type='c',
            style_name=new_style,
            network=network
        )

    if edge_colors:
        p4c.style_mappings.set_edge_color_mapping(
            'edge-priority',
            table_column_values = [f'direct path ({source})' for source in sources],
            colors=edge_colors,
            mapping_type='d',
            style_name=new_style,
            network=network
        )

    p4c.styles.set_visual_style(new_style)


def layout_from_coords(network, table):
    '''
    network : int
        Cytoscape suid, where to apply the layout
    table: DataFrame
        Pandas DataFrame, index as node names, and columns 'x' and 'y' with coordinates

    '''

    p4c.load_table_data(table, network=network)

    current_style = p4c.styles.get_current_style(network=network)

    tmp_style = 'tmp-layout'
    p4c.copy_visual_style(current_style, tmp_style)
    p4c.set_visual_style(tmp_style, network=network)

    p4c.update_style_mapping(tmp_style, p4c.map_visual_property('NODE_X_LOCATION', 'x', 'p'))
    p4c.update_style_mapping(tmp_style, p4c.map_visual_property('NODE_Y_LOCATION', 'y', 'p'))

    p4c.delete_visual_style(tmp_style)

    p4c.tables.delete_table_column('x')
    p4c.tables.delete_table_column('y')

    p4c.set_visual_style(current_style, network=network)

def add_custom_png(network, create_png, style=None):
    '''
    network : int
        Cytoscape suid, where to apply the layout
    create_png: function
        function that creates a png per node

    '''

    nodes = p4c.get_all_nodes(network=network)

    node_pngs = {}
    for node in nodes:
        node_png_fname = create_png(node)
        if node_png_fname:
            node_pngs[node] = {'fig_location':f"file:{str(node_png_fname.absolute())}"}

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

    # fit content ignores edges that may extend
    # past the node boundaries
    # Reported bug: CSD-979
    p4c.network_views.fit_content(network=network)
    p4c.network_selection.clear_selection(type='both', network=network)

    p4c.network_views.export_image(
        filename=str(filename.absolute()),
        type=format,
        network=network,
        overwrite_file=True,
        resolution=600
        # all_graphics_details=True,
        # hide_labels=False
    )

def export_collection_to_pdf(collection_suid, filename, font_size=20, font='Helvetica'):
    '''requires pdf libraries'''

    from pypdf import PdfWriter, PdfReader
    from pdfCropMargins import crop

    from reportlab.pdfgen import canvas

    def create_text_pdf(text, mb, font, font_size):

        packet = io.BytesIO()

        can = canvas.Canvas(packet, pagesize=(mb.width, mb.height))
        try:
            can.setFont(font, font_size)
        except KeyError as e:
            print(f"Font {font} not available. Using {can._fontname}. ")
            font = can._fontname
            can.setFont(font, font_size)

        text_width = canvas.pdfmetrics._fonts[font].stringWidth(text, size=font_size)
        # rtl_fonts['Helvetica'].stringWidth(text, size=font_size)

        center = (mb.right - mb.left)/2 + mb.left
        text_left = center - text_width/2
        bottom = mb.bottom + font_size/4


        can.drawString(text_left, bottom, text)
        can.save()

        packet.seek(0)
        new_pdf = PdfReader(packet)

        return new_pdf

    networks = sorted(p4c.collections.get_collection_networks(collection_suid))

    i = 0
    while (filename.parent / f"tmp{i}").exists():
        i += 1
    tmp_folder = filename.parent / f"tmp{i}"
    tmp_folder.mkdir()


    writer = PdfWriter()
    for network in networks:
        network_name = p4c.get_network_name(network)
        print(network, network_name)
        pdf = tmp_folder / f"{re.sub('[^a-zA-Z0-9]+', '_', network_name).strip('_')}.pdf"
        export_network(network, pdf, format="pdf")

        cropped_pdf1 = tmp_folder / (pdf.stem + "_cropped" + pdf.suffix)
        cropped_pdf2 = tmp_folder / (pdf.stem + "_2nd_cropped" + pdf.suffix)

        crop(['--noundosave', '-o', str(cropped_pdf1), str(pdf)])
        crop(['--noundosave', '-a4', '0', f'-{font_size}', '0', '0', '-o', str(cropped_pdf2), str(cropped_pdf1)])

        reader = PdfReader(cropped_pdf2)
        page = reader.pages[0]
        mb = page.cropbox
        caption_pdf = create_text_pdf(network_name, mb, font, font_size)

        new_page = caption_pdf.pages[0]
        new_page.cropbox = page.cropbox
        page.merge_page(caption_pdf.pages[0])

        writer.add_page(page)

        pdf.unlink()
        cropped_pdf1.unlink()
        cropped_pdf2.unlink()


    with open(filename, "wb") as output_stream:
        writer.write(output_stream)
    writer.close()

    tmp_folder.rmdir()

    print(f'Collection save to {filename}')