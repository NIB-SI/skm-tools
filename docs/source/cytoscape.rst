=====================
Cytoscape automation
=====================

:mod:`skm_tools.cytoscape_utils` sends networks and analysis results to
`Cytoscape <https://cytoscape.org>`_, to look at, style and export them. It uses
`py4cytoscape <https://py4cytoscape.readthedocs.io>`_, which controls a running Cytoscape
through its REST interface (CyREST).

Getting started
===============

You need:

- the ``cytoscape`` extra: ``pip install "skm-tools[cytoscape]"`` (see :doc:`installation`),
- Cytoscape, open on the same computer. Python talks to it while it runs, so leave it open.

Check the connection:

.. code-block:: python

   import py4cytoscape as p4c
   from skm_tools import cytoscape_utils as cu

   p4c.cytoscape_ping()   # "You are connected to Cytoscape!"

A few Cytoscape terms used below:

- **SUID**: Cytoscape's id for a network (and for each node and edge). Functions here
  return the SUID of the networks they create, and take it as ``network``. SUIDs are not
  kept: they change when a session is saved and opened again. To find a network again,
  use its name, e.g. ``p4c.get_network_suid("PSS")``.
- **Collection**: a group of networks, e.g. a network and its subnetworks; shown as a
  folder in Cytoscape's *Network* panel. The networks of a collection share their node and
  edge tables: a column added in one network is there in all of them.
- **Visual style**: how a network is drawn (shapes, colours, labels), from the values in
  the node and edge tables. Several networks can share a style.
- **Node and edge tables**: the node and edge attributes, one column per attribute.
  Loading data into a new column, and using that column in a style, is how data is shown
  on the network.

Loading a network
=================

:func:`~skm_tools.cytoscape_utils.load_network` sends a networkx graph to Cytoscape. The
node and edge attributes become table columns. ``style`` is optional: it applies one of
the SKM styles; without it, Cytoscape uses its default style.

.. code-block:: python

   suid = cu.load_network(pss, title="PSS", style="skm")

The SKM styles are bundled with skm-tools, and imported into Cytoscape the first time
they are used (also with :func:`~skm_tools.cytoscape_utils.apply_builtin_style`). CKN and
PSS use the same node classes (``node_type``) and edge attributes, so they share a style:

- ``"skm"`` (*SKM*): for CKN, the PSS interaction network and the PSS gene networks.
  Nodes are coloured and shaped by ``node_type``, in the colours of the PSS Explorer, and
  labelled with ``display_label``. Edges are coloured, and get their arrow, by
  ``interaction`` (positive: green, negative: red, unknown: grey); mutual (``directed``
  False) edges are dashed, and better supported edges (lower ``rank``) are thicker.
- ``"skm-reactions"`` (*SKM-reactions*): for the PSS reaction graph. The same nodes, plus
  the reactions (blue); edges as in the PSS Explorer, by ``edge_type`` (activation: green,
  inhibition: red, substrate and product: grey), labelled with the participant's
  ``role``.

To change a style for one analysis, copy it first (see
:func:`~skm_tools.cytoscape_utils.copy_style`), so the original stays as it is.

Highlighting
============

Highlighting changes the look of single nodes and edges, on top of the style (Cytoscape
calls these *bypasses*: they stay when the style changes, until cleared in Cytoscape).

- :func:`~skm_tools.cytoscape_utils.highlight_nodes`: fill, label and border colour,
  border width, size.
- :func:`~skm_tools.cytoscape_utils.highlight_edges`: colour and line width.
- :func:`~skm_tools.cytoscape_utils.highlight_path`: the nodes and edges of a path.

To highlight several paths in different colours, pass what was already highlighted as
``skip_nodes`` and ``skip_edges``, so later paths don't paint over earlier ones. The
example below highlights the shortest paths from JA to SA in red, and then those from ABA
to SA in blue; nodes and edges on both (e.g. SA itself) stay red:

.. code-block:: python

   from skm_tools.paths import get_paths

   done_nodes, done_edges = [], []
   for paths, colour in [(get_paths(pss, "JA", "SA"), "#E41A1C"),
                         (get_paths(pss, "ABA", "SA"), "#377EB8")]:
       for p in paths:
           nodes, edges = cu.highlight_path(p, colour, skip_nodes=done_nodes,
                                            skip_edges=done_edges, network=suid)
           done_nodes += nodes
           done_edges += edges

:func:`~skm_tools.cytoscape_utils.apply_shortest_paths_style` shows the results of several
path searches to the same target in one go, through a style instead of highlights. It
adds three columns to the tables: for the nodes on the paths, their distance to the
target (``distance-to-target``) and the source whose paths they are on
(``node-path-source``; the first source, if on several), and for the edges on the paths,
``edge-priority`` (``"direct path (<source>)"``). It then makes a copy of the network's
style, ``<style>-shortest-paths-query``, in which node colour shows the distance to the
target (from the first of ``node_colors``, at the target, to the last), and edge colour
the source (one of ``edge_colors`` per source):

.. code-block:: python

   sources = ["JA", "ABA"]
   path_lists = [get_paths(pss, source, "SA") for source in sources]
   cu.apply_shortest_paths_style(sources, path_lists, "SA", pss,
                                 edge_colors=["#E41A1C", "#377EB8"],
                                 node_colors=["#FFFFFF", "#FDAE61", "#D7191C"], network=suid)

Subnetworks
===========

Subnetworks are new networks in the same collection, with part of a network; they keep
its layout and style.

- :func:`~skm_tools.cytoscape_utils.subnetwork_node_induced`: the given nodes and *every*
  edge between them.
- :func:`~skm_tools.cytoscape_utils.subnetwork_edge_induced`: only the given edges (and
  their nodes).
- :func:`~skm_tools.cytoscape_utils.subnetwork_edge_induced_from_paths`: only the edges
  along paths, e.g. the paths found by :func:`~skm_tools.paths.get_paths`.
- :func:`~skm_tools.cytoscape_utils.subnetwork_neighbours`: the given nodes and their
  first neighbours.

.. code-block:: python

   paths = get_paths(pss, "JA", "SA")
   paths_suid = cu.subnetwork_edge_induced_from_paths(paths, pss, suid, name="JA to SA")

:func:`~skm_tools.cytoscape_utils.get_or_create_subnetwork` returns the network of that
name if it already exists, so notebook cells can be run again without making duplicates.

For layouts, use Cytoscape's own (``p4c.layout_network``), or place the nodes at
coordinates computed in Python (e.g. a networkx or graphviz layout) with
:func:`~skm_tools.cytoscape_utils.layout_from_coords`.

Showing data as images on nodes
===============================

Small plots, e.g. a heatmap of a gene's logFC per time point, can be shown next to each
node. This takes three steps:

1. **Find the image of each node.** :func:`~skm_tools.cytoscape_utils.match_files_to_nodes`
   matches the image files in a folder to the nodes by name.
2. **Load the image locations into a node table column.**
   :func:`~skm_tools.cytoscape_utils.load_node_images`.
3. **Show that column in a style.** :func:`~skm_tools.cytoscape_utils.show_node_images`.

.. code-block:: python

   images = cu.match_files_to_nodes("heatmaps", pss.nodes())
   cu.load_node_images(images, "heatmap", network=suid)
   cu.show_node_images("SKM", "heatmap")

Matching files to nodes
-----------------------

A file matches a node if the file name (without the extension) is the node name, with the
characters that can't be in file names replaced by ``_``. For example:

- ``WRKY33_fc00166_.png`` is the image of ``WRKY33[fc00166]``,
- ``11-_12-OH-JA.png`` is the image of ``11-/12-OH-JA``,
- ``AT2G38470.png`` is the image of the gene ``AT2G38470`` (gene ids stay as they are).

:func:`~skm_tools.cytoscape_utils.node_file_key` gives the file name for a node name, to
use when making the images. For files that aren't named after their node, give the node
name in ``aliases``, e.g. ``aliases={"Pro": "Proline accumulation"}``. Files that match no
node are printed, so they can be checked.

Several images per node
-----------------------

A node can show up to nine images, each in its own *slot*, e.g. transcriptomics above the
node and metabolomics below it. Load each kind into its own column, and show each column
in its own slot and position:

.. code-block:: python

   for omics, slot, position in [("transcriptomics", 1, "above"), ("metabolomics", 2, "below")]:
       images = cu.match_files_to_nodes(f"heatmaps/{omics}", pss.nodes())
       cu.load_node_images(images, f"image_{omics}", network=suid)
       cu.show_node_images("SKM", f"image_{omics}", slot=slot, position=position)

The position puts a side of the image against a side of the node, centred:

- ``"below"``: the top of the image against the bottom of the node,
- ``"above"``: the bottom of the image against the top of the node,
- ``"left"`` and ``"right"``: the image's right or left side against the node's left or
  right side,
- ``"center"``: the image on the middle of the node.

Other positions can be given in Cytoscape's own form, ``"<node point>,<image
point>,c,<x offset>,<y offset>"``, with the points as compass directions (``N``, ``NE``,
``E``, ..., ``C`` for the centre): ``"below"`` is ``"S,N,c,0.00,0.00"``, the image's north
(top) point on the node's south (bottom) point. ``size`` sets the size of the images.

Comparing conditions
--------------------

To compare conditions (e.g. heat and drought), show each condition on its own copy of the
network. Each copy also needs its own style: a style says which column's images to show
(e.g. "show the images in column ``image_heat`` below the nodes"), so all networks with
the same style show the same images. The example loads each condition's images, makes a
copy of the network and a copy of the style for it, and sets that style to show the
condition's column:

.. code-block:: python

   for condition in ["heat", "drought"]:
       images = cu.match_files_to_nodes(f"heatmaps/{condition}", pss.nodes())
       cu.load_node_images(images, f"image_{condition}", network=suid, unique_dir="images")

       copy = cu.clone_network(suid, name=condition, collection=f"PSS - {condition}")
       style = cu.copy_style("SKM", f"SKM-{condition}", networks=[copy])
       cu.show_node_images(style, f"image_{condition}")

:func:`~skm_tools.cytoscape_utils.clone_network` copies the network into a new collection,
which has its own tables: load the images before cloning (as here), or into the clone.
:func:`~skm_tools.cytoscape_utils.get_or_create_subnetwork` makes a copy in the same
collection instead, which shares the tables.

Image caching
-------------

Cytoscape keeps the images it has shown in a cache, and how it uses the cache is hard to
predict: it can show an image it loaded before instead of the right one, e.g. an image of
another file with the same name, or the old version of a file that has since been
changed. Two things help:

- Keep the number of images down, e.g. only for the nodes and conditions you need.
- Give ``unique_dir`` to :func:`~skm_tools.cytoscape_utils.load_node_images`. The images
  are then first copied to that folder, each under a new name that no file had before
  (the file name plus a random id), every time, so Cytoscape has to load them again. The
  copies are not deleted afterwards; empty the folder from time to time.

SVG images
----------

Cytoscape can show SVG images, but sizes them differently from PNG images: in our tests,
SVG images were drawn larger than PNG images of the same ``size``. PNG images are
therefore easier to place. To convert SVG images (e.g. plots made in R) to PNG, use
``to_png=True`` together with ``unique_dir``: the SVG images are converted while they are
copied. The conversion needs the ``cairosvg`` package, which is not installed with
skm-tools: install it with ``pip install cairosvg``.

Images made on the fly
----------------------

:func:`~skm_tools.cytoscape_utils.add_custom_png` does the three steps in one call. Instead
of a folder of images, it takes a function you write, which gets a node name and returns
the file name of that node's image (or None for no image). The function can make the
image (e.g. plot the node's data with matplotlib and save it), or only return the name of
an existing file. The images are shown below the nodes.

Charts instead of images
------------------------

Cytoscape can also draw simple charts itself, with the
`enhancedGraphics <https://apps.cytoscape.org/apps/enhancedgraphics>`_ app.
:func:`~skm_tools.cytoscape_utils.chart_column` makes a chart for each node from a table
with a row per node: by default a bar chart, with one bar per given column, in the given
colours, and the value axis from ``value_range``. Other chart types are set with
``chart``, e.g. ``chart="linechart"`` or ``chart="heatstripchart"`` (see the
enhancedGraphics documentation). Load the charts as a node table column, and show them
like images:

.. code-block:: python

   # df: logFC per gene (rows) and time point (columns)
   charts = cu.chart_column(df, ["10min", "30min", "1h"], "#E41A1C", value_range=(-2, 2))
   p4c.load_table_data(charts.to_frame("chart"), network=suid)
   cu.show_node_images("SKM", "chart")

Exporting
=========

- :func:`~skm_tools.cytoscape_utils.export_network`: a network as an image (PDF, PNG, SVG,
  ...), zoomed to fit and with nothing selected.
- :func:`~skm_tools.cytoscape_utils.export_collection`: every network of a collection, one
  image each.
- :mod:`skm_tools.cytoscape_pdf_utils`: every network of a collection to PDF, one file per
  network or one document with captions (needs the ``pdf`` extra).

Tips
====

Running scripts in the background
---------------------------------

Cytoscape brings its window to the front for many operations, which makes it hard to do
anything else while a script runs. On Linux, run Cytoscape in its own window with
`Xephyr <https://freedesktop.org/wiki/Software/Xephyr/>`_, a display inside a window:
Cytoscape can only take the focus inside it.

.. code-block:: bash

   Xephyr :2 -screen 1600x1000 -resizeable &
   cd /path/to/Cytoscape && DISPLAY=:2 bash -c 'sleep infinity | ./cytoscape.sh' &

(``sleep infinity |`` keeps Cytoscape's console input open; without it, Cytoscape shuts
down right away when started in the background.)

Messages from py4cytoscape
--------------------------

py4cytoscape prints the text of every error, also of errors that are expected and
handled (e.g. while waiting for Cytoscape to finish creating a network). This output is
switched off when :mod:`skm_tools.cytoscape_utils` is imported; errors are still raised as
usual. :func:`~skm_tools.cytoscape_utils.silence_py4cytoscape` with ``False`` switches it
back on.

Housekeeping
------------

- :func:`~skm_tools.cytoscape_utils.delete_other_networks`: delete every network but one,
  e.g. to start over while trying things out.

How nodes and edges are found
-----------------------------

The functions here find nodes in Cytoscape by their ``name`` (the networkx node id), and
edges by their ``name``, ``"source (interaction) target"``, both set by
:func:`~skm_tools.cytoscape_utils.load_network`. Networks loaded in other ways need the
same names. All parallel edges between two nodes are found (and highlighted) together.

Useful links
============

- `py4cytoscape documentation <https://py4cytoscape.readthedocs.io/en/latest/>`_
- `Cytoscape Automation wiki <https://github.com/cytoscape/cytoscape-automation/wiki>`_
- `Network layout notebook <https://github.com/cytoscape/cytoscape-automation/blob/master/for-scripters/Python/network-layout.ipynb>`_.
  The Copycat layout app can be installed from the Cytoscape App store; yFiles layouts
  can't be used from scripts.
