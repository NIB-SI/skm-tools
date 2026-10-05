========
Analysis
========

The analysis functions work on any networkx graph, so on PSS and CKN alike, and return
networkx graphs or plain Python data.

Shortest paths
==============

:func:`~skm_tools.paths.get_paths` finds the shortest paths from one or more sources to
one or more targets. Nodes not in the graph, and pairs without a path, are skipped.

.. code-block:: python

   from skm_tools.paths import get_paths, path_subgraph

   paths = get_paths(pss, sources=["ABA", "JA", "SA"], targets="RD29[fc00192]")
   # undirected search, keeping only the overall shortest paths to each target
   paths = get_paths(ckn, sources, targets, directed=False, shortest_overall=True)

   # the network made of the paths' edges
   g_paths = path_subgraph(paths, pss)

Neighbourhoods
==============

:func:`~skm_tools.neighbors.get_neighborhood` returns the subgraph around some nodes, to
any depth, following edges in both directions (default), downstream (``"out"``) or
upstream (``"in"``). Every node gets a ``distance`` attribute.

.. code-block:: python

   from skm_tools.neighbors import get_neighborhood

   g1 = get_neighborhood(ckn, "AT2G38470")                     # first neighbours
   g2 = get_neighborhood(ckn, "AT2G38470", depth=2, direction="out", induced=False)

With ``induced=False`` only the edges followed outwards are kept, not every edge between
the nodes found.

Experimental data
=================

:func:`~skm_tools.experimental_data.overlay_experimental_data` adds logFC and p-values as
node attributes. It returns a copy, so several experiments can be overlaid on separate
copies of the same network.

.. code-block:: python

   import pandas as pd
   from skm_tools.experimental_data import overlay_experimental_data

   df = pd.read_csv("deg.tsv", sep="\t", index_col=0)    # indexed by gene id

   # CKN nodes are genes: match on the node id
   ckn_h2o2 = overlay_experimental_data(ckn, df, "logFC", "padj", prefix="H2O2 ")

   # PSS nodes are functional clusters: match on their genes; the most significant gene is used
   pss_h2o2 = overlay_experimental_data(pss, df, "logFC", "padj", match_attribute="ath_homologues")

Minimum cuts
============

:func:`~skm_tools.cuts.get_cutset` finds the smallest set of edges that disconnects a set
of sources from a set of targets. Set a ``capacity`` on every edge first.

Saving and serialising
======================

:mod:`skm_tools.persistence` saves and loads networks:

.. code-block:: python

   from skm_tools.persistence import save_graph, load_graph

   save_graph(ckn, "ckn-filtered.json", format="json")
   ckn = load_graph("ckn-filtered.json", format="json")

.. list-table:: Full CKN (26k nodes, 899k edges)
   :header-rows: 1

   * - Format
     - Save
     - Load
     - Size
     - Notes
   * - ``pickle`` (default)
     - 2 s
     - 2 s
     - 101 MB
     - Python only; keeps everything; never load untrusted files
   * - ``json``
     - 15 s
     - 6 s
     - 193 MB
     - Portable and safe; keeps lists
   * - ``graphml``
     - 43 s
     - 84 s
     - 262 MB
     - For other network tools; lists are joined with ``;``, empty values left out

:mod:`skm_tools.serialize` converts graphs and paths to JSON-safe dicts and lists
(:func:`~skm_tools.serialize.graph_to_dict`, :func:`~skm_tools.serialize.path_to_dict`).
