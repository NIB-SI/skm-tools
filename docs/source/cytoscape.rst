=====================
Cytoscape automation
=====================

:mod:`skm_tools.cytoscape_utils` shows networks and analysis results in
`Cytoscape <https://cytoscape.org>`_ through `py4cytoscape <https://py4cytoscape.readthedocs.io>`_.
It needs the ``cytoscape`` extra and a running Cytoscape.

.. code-block:: python

   from skm_tools import cytoscape_utils as cu
   from skm_tools.paths import get_paths

   suid = cu.load_network(pss, title="PSS", style="pss")

   paths = get_paths(pss, "JA", "SA")
   done_nodes, done_edges = [], []
   for p in paths:
       nodes, edges = cu.highlight_path(p, "#E41A1C", skip_nodes=done_nodes, skip_edges=done_edges, network=suid)
       done_nodes += nodes
       done_edges += edges

   paths_suid = cu.subnetwork_edge_induced_from_paths(paths, pss, suid, name="JA to SA")

Bundled styles
==============

:func:`~skm_tools.cytoscape_utils.apply_builtin_style` applies one of the SKM styles:

- ``"pss"`` (PSS-default): node shape by ``node_type``, fill colour by ``pathway``, arrow
  by ``interaction`` (positive: arrow, negative: T, unknown: diamond), and dashed lines for
  mutual (``directed`` False) edges.
- ``"ckn"`` (CKN-default).

How nodes and edges are found
=============================

Nodes are found by their Cytoscape ``name`` (the networkx node id) and edges by their
``name``, ``"source (interaction) target"``, as set when loading the network. All
parallel edges between two nodes are matched together.

PDF export
==========

:mod:`skm_tools.cytoscape_pdf_utils` exports every network of a Cytoscape collection to
PDF, either one file per network or a single captioned document. It needs the ``pdf`` extra.

Useful links
============

- `Cytoscape Automation wiki <https://github.com/cytoscape/cytoscape-automation/wiki>`_
- `py4cytoscape documentation <https://py4cytoscape.readthedocs.io/en/latest/>`_
- `Network layout notebook <https://github.com/cytoscape/cytoscape-automation/blob/master/for-scripters/Python/network-layout.ipynb>`_
  (Cytoscape Automation wiki). The Copycat layout app can be installed from the Cytoscape
  App store; yFiles layouts don't support automation.
