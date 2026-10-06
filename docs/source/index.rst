==========
skm-tools
==========

Python tools for systems biology analysis of plant stress responses with the
`Stress Knowledge Map <https://skm.nib.si>`_ (SKM).

| **Homepage:** https://github.com/NIB-SI/skm-tools
| **Documentation:** https://nib-si.github.io/skm-tools
| **License:** GPL-3.0

Overview
--------

SKM contains knowledge on molecular interactions in plants, integrated from a wide
diversity of sources into a single, freely available entrypoint. It has two knowledge graphs:

- **PSS** (Plant Stress Signalling): a highly curated, detailed mechanistic model of plant
  stress signalling, mostly compiled from targeted biochemical studies in the literature,
  with over 800 reactions.
- **CKN** (Comprehensive Knowledge Network): molecular interactions in the plant cell,
  mostly from high-throughput experiments, with over 26,000 molecules and 480,000 interactions.

skm-tools loads both as `networkx <https://networkx.org>`_ graphs, so the whole networkx
toolbox is available, and adds functions for the analyses SKM is typically used for:

- **Loading:** PSS (reaction graph, interaction network, per-species gene networks) and CKN,
  with node and edge annotations
- **Filtering:** by node type, species, tissue, MapMan bin, edge rank or type
- **Shortest paths** between sources and targets of interest, directed or undirected
- **Neighbourhoods** of nodes of interest, to any depth
- **Minimum cuts** between sets of nodes
- **Translation** of Arabidopsis networks to other species
- **Saving** networks (pickle, JSON, GraphML) and JSON-safe conversion of graphs and paths
- **Cytoscape automation:** load, style, highlight, extract subnetworks and export images,
  through `py4cytoscape <https://py4cytoscape.readthedocs.io>`_. Includes - **experimental data overlay** (e.g. heatmaps with expression logFC and p-values as node annotations)


.. toctree::
   :maxdepth: 1
   :caption: Manual and guides

   installation
   pss
   ckn
   translations
   analysis
   cytoscape

Citation
========

If you use skm-tools (or SKM) in your work, please cite:

   Bleker C, Ramšak Ž, Bittner A, Podpečan V, Zagorščak M, Wurzinger B, Baebler Š, Petek M,
   Križnik M, van Dieren A, Gruber J, Afjehi-Sadat L, Weckwerth W, Županič A, Teige M,
   Vothknecht UC, Gruden K (2024).
   *Stress Knowledge Map: A knowledge graph resource for systems biology analysis of plant
   stress responses*. Plant Communications.
   https://doi.org/10.1016/j.xplc.2024.100920


Indices
=======

* :ref:`genindex`
* :ref:`modindex`
