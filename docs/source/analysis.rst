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

:func:`~skm_tools.experimental_data.add_experimental_data` adds logFC and p-values as
node attributes. It returns a copy, so several experiments can be added to separate
copies of the same network.

.. code-block:: python

   import pandas as pd
   from skm_tools.experimental_data import add_experimental_data

   df = pd.read_csv("deg.tsv", sep="\t", index_col=0)    # indexed by gene id

   # CKN nodes are genes: match on the node id
   ckn_h2o2 = add_experimental_data(ckn, df, "logFC", "padj", prefix="H2O2 ")

   # PSS nodes are functional clusters: match on their genes; the most significant gene is used
   pss_h2o2 = add_experimental_data(pss, df, "logFC", "padj", match_attribute="ath_homologues")

.. _mapman:

MapMan annotations
==================

PSS and CKN nodes have a ``mapman`` attribute: the MapMan bins (GoMapMan 2, MapMan4) of their
Arabidopsis genes, each as ``<bin code>_<bin name>``. :func:`~skm_tools.annotations.get_nodes_by_mapman`
finds the nodes in any of the given bins, by default including their sub-bins; the same bins
work for both networks.

.. code-block:: python

   from skm_tools.annotations import get_nodes_by_mapman

   # 26.11: External stimuli response.pathogen
   pss_pathogen = get_nodes_by_mapman(pss, "26.11")
   ckn_pathogen = get_nodes_by_mapman(ckn, "26.11")

Minimum cuts
============

:func:`~skm_tools.cuts.get_cutset` finds the smallest set of edges that disconnects a set
of sources from a set of targets, e.g. the bottlenecks between a signal and a response. It
needs a :class:`networkx.DiGraph` (for PSS, run :func:`~skm_tools.pss.simplify_pss` first)
with a ``capacity`` on every edge; a capacity of 1 counts edges. It prints the maximum
flow, and returns the cut edges:

.. code-block:: python

   import networkx as nx
   from skm_tools.cuts import get_cutset
   from skm_tools.pss import simplify_pss

   simple = simplify_pss(pss)
   nx.set_edge_attributes(simple, 1, "capacity")

   get_cutset(["flg22"], ["WRKY33[fc00166]"], simple)
   # max_flow = 1
   # [('MPK3,6[fc00308]', 'WRKY33[fc00166]')]

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
