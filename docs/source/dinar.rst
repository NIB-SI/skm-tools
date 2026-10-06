=====
DiNAR
=====

`DiNAR <https://github.com/NIB-SI/DiNAR>`_ (Differential Network Analysis in R) shows
differential expression (or other omics data) of several conditions or time points on a
prior knowledge network, one cluster of nodes at a time. A PSS gene network, CKN, or a
subnetwork of either can be used as the knowledge network (as a DiNAR *Custom Network* in
the form of a nodes table and an edges table).

From Python
===========

:func:`~skm_tools.dinar.write_dinar` writes a network as DiNAR's nodes and edges tables
(:func:`~skm_tools.dinar.to_dinar` returns them as :class:`pandas.DataFrame`):

.. code-block:: python

   from skm_tools.dinar import write_dinar
   from skm_tools.pss import pss_gene_network_to_networkx

   pss_ath = pss_gene_network_to_networkx(
       "pss-gene-network-ath-edges-public.tsv",
       "pss-gene-network-ath-nodes-public.tsv",
   )
   write_dinar(pss_ath, "pss-ath-dinar-nodes.txt", "pss-ath-dinar-edges.txt", clusters="pathway")

Nodes get ``shortName`` from ``display_label``, ``shortDescription`` from ``description``
and ``MapManBin`` from ``mapman`` (see :ref:`mapman`); edges get ``reactionType`` from
``interaction`` (``positive-influence``, ...). Missing values are ``-``, as DiNAR expects.

Clusters
--------

DiNAR shows one cluster at a time, so choose clusters that make sense to look at together.
The argument ``clusters`` can be a node attribute (e.g. ``"pathway"``; for a gene in
several functional clusters, its first one is used) or a ``{node: cluster}`` mapping, e.g.
from community detection:

.. code-block:: python

   import networkx as nx

   communities = nx.community.louvain_communities(nx.Graph(ckn_subnetwork), seed=0)
   clusters = {n: i for i, community in enumerate(communities) for n in community}
   write_dinar(ckn_subnetwork, "nodes.txt", "edges.txt", clusters=clusters)

Without ``clusters``, all nodes will be in one cluster. Nodes without a cluster (e.g.
without a ``pathway``) get ``clusterID`` 0, which DiNAR doesn't show.

Coordinates
-----------

``positions`` is a ``{node: (x, y)}`` mapping, as returned by the networkx layout functions,
or a node attribute holding them. By default, the ``pos`` node attribute is used if every
node has one (the networkx convention, ``nx.set_node_attributes(g, pos, "pos")``);
otherwise each cluster gets a spring layout, with the clusters on a grid. To use a layout
made in Cytoscape (see :doc:`cytoscape`):

.. code-block:: python

   import py4cytoscape as p4c

   xy = p4c.get_node_position(network=suid)
   # Cytoscape's y axis points down
   positions = {n: (x, -y) for n, x, y in zip(xy.index, xy["x"], xy["y"])}
   write_dinar(pss_ath, "nodes.txt", "edges.txt", clusters="pathway", positions=positions)

In DiNAR
--------

Run DiNAR (`online <https://nib-si.shinyapps.io/DiNAR>`__ or locally, see its
`README <https://github.com/NIB-SI/DiNAR>`_), choose *Custom network*, upload the nodes and
edges tables, and then the expression tables of your conditions.

Without code
============

Networks from the SKM explorers can be turned into DiNAR tables with DiNAR's
GMM-SKM-KnetMiner app:

1. Get the (sub)network of interest: download the result of a query in the
   `PSS <https://skm.nib.si/biomine/>`_ or `CKN <https://skm.nib.si/ckn/>`_ Explorer,
   exporting the nodes and edges as ``.csv``. In the case of PSS, the entire network is also
   visualisable (download from `skm.nib.si/downloads <https://skm.nib.si/downloads>`_).
2. Run the `GMM-SKM-KnetMiner <https://github.com/NIB-SI/DiNAR/tree/master/subApps/GMM-SKM-KnetMiner>`_
   app, `online <https://nib-si.shinyapps.io/GMM-SKM-KnetMiner/>`__ or locally
   (`script <https://github.com/NIB-SI/DiNAR/blob/master/subApps/GMM-SKM-KnetMiner/scripts/CustomNetwork_KnetMiner-SKM.R>`_).
3. Load the `files <https://github.com/NIB-SI/DiNAR/tree/master/subApps/GMM-SKM-KnetMiner/input>`_.
4. Follow the instructions, and save the Custom Network tables (nodes and edges).
5. Load the `Custom Network <https://github.com/NIB-SI/DiNAR/tree/master/subApps/GMM-SKM-KnetMiner/output-for-DiNAR>`_
   and the expression tables into DiNAR.
6. Export the results as needed
   (`examples <https://github.com/NIB-SI/DiNAR/tree/master/subApps/GMM-SKM-KnetMiner/output-from-DiNAR>`_).

In tables made by hand, replace empty values with ``-``, and keep ``shortName`` short.
