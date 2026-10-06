===
PSS
===

The Plant Stress Signalling model (PSS) is a curated, mechanistic model of plant stress
signalling; see the `PSS documentation <https://skm.nib.si/documentation/pss-explore>`_ on
the SKM website for what it contains and how it is built. Everything here is in
:mod:`skm_tools.pss`.

Loading
=======

PSS is exported in three network forms, each as an edge file and a node file (tab-separated,
with a header). The exports are made by `skm-pss-export <https://github.com/NIB-SI/skm-pss-export>`_,
which describes the formats and every column in detail.

.. list-table::
   :header-rows: 1
   :widths: 25 45 30

   * - Export
     - Content
     - Loader
   * - Reaction graph
     - Entities **and reactions** as nodes, one edge per reaction participant
       (participant → reaction for inputs and modifiers, reaction → participant for products).
       The lossless form of PSS. Edges are keyed by ``role``.
     - :func:`~skm_tools.pss.pss_reaction_graph_to_networkx`
   * - Interaction network
     - Entities as nodes; edges are entity → entity influences through the reactions
       (``interaction``: ``positive-influence``, ``negative-influence`` or ``unknown-influence``).
       Edges are keyed by ``reaction_id``.
     - :func:`~skm_tools.pss.pss_interaction_network_to_networkx`
   * - Gene network (per species)
     - The interaction network with functional clusters expanded into the genes of one
       species (``node_type`` ``gene``), for each species in :doc:`translations`.
     - :func:`~skm_tools.pss.pss_gene_network_to_networkx`

All three are loaded as :class:`networkx.MultiDiGraph`, as several reactions can link the same two
nodes. Mutual influences (``directed`` False, e.g. between binding partners) are listed in
both directions. Every edge has ``rank`` 0, CKN's rank for curated PSS interactions, so PSS
and CKN can be combined and filtered alike. The networks are downloaded from `skm.nib.si <https://skm.nib.si>`_ if they don't exist locally:

.. code-block:: python

   from skm_tools.pss import pss_interaction_network_to_networkx

   pss = pss_interaction_network_to_networkx(
       "pss-interaction-network-edges-public.tsv",
       "pss-interaction-network-nodes-public.tsv",
   )
   pss.nodes["WRKY33[fc00166]"]["display_label"]   # 'WRKY33'

To download a gene network, pass the species code too, e.g.
``pss_gene_network_to_networkx(edge_path, node_path, species="ath")``.

Node attributes include ``node_type`` (the PSS class: ``PlantCoding``, ``Metabolite``,
``Complex``, ..., see the `PSS database schema <https://skm.nib.si/documentation/pss-db#node-labels>`_),
``display_label``, ``short_name``, ``synonyms``, ``pathway``, ``mapman`` (see
:ref:`mapman`), ``components`` and ``component_cluster_ids`` (for complexes), and
``<species>_homologues``. In the files, lists are joined with ``|`` (names can contain ``,``,
e.g. ``AHK2,3,4``, and gene symbols ``;``, e.g. ``PIP1;3``), and the loaders return them as
Python lists. Empty values are ``None``, and booleans
are ``True``/``False``.

Species
=======

PSS is species-independent: its plant nodes are `functional clusters
<https://skm.nib.si/documentation/pss-explore#FC>`_, whose genes are listed per species in the
``<species>_homologues`` node attributes of the reaction graph and the interaction network. A
gene network expands the clusters into the genes of one species, and only has the reactions
whose functional clusters all have genes in that species. For the species and their codes,
see :doc:`translations`.

Filtering and simplifying
=========================

The filtering functions change the graph in place (make a copy first, ``g.copy()``, to keep
the original) and return the reasons nodes were removed.

:func:`~skm_tools.pss.filter_pss_nodes` keeps only some node types, and/or removes plant
nodes without genes in the given species. Complexes that lose a component are removed too,
and so are nodes left without edges:

.. code-block:: python

   from skm_tools.pss import filter_pss_nodes

   reasons = filter_pss_nodes(pss, species=["stu"])
   reasons = filter_pss_nodes(pss, node_types=["PlantCoding", "Complex", "Metabolite"])

:func:`~skm_tools.pss.remove_deadend_complexes` removes complexes that influence nothing
(no outgoing edges), which are dead ends in directed analyses.

:func:`~skm_tools.pss.simplify_pss` merges the parallel edges (one per reaction) between two
nodes into one edge, and returns a new :class:`networkx.DiGraph`. The merged edge keeps all
``reaction_id`` values, and one value of every other attribute (with a printed warning if
they differ, e.g. a positive and a negative influence). To keep such edges apart, pass
``split_on_attrs=["interaction"]``; the result is then a :class:`networkx.MultiDiGraph`.

:func:`~skm_tools.pss.remove_and_rewire` removes nodes while connecting each of their
upstream nodes to each of their downstream nodes, e.g. to hide intermediate steps. It needs
a simple :class:`networkx.DiGraph`, so run :func:`~skm_tools.pss.simplify_pss` first.

:func:`~skm_tools.pss.remove_duplicated_binding_edges` keeps one direction of each mutual
binding edge, for a less cluttered drawing. Don't use it before path or neighbourhood analysis.

.. code-block:: python

   from skm_tools.pss import simplify_pss, remove_and_rewire

   simple = simplify_pss(pss)
   remove_and_rewire(simple, ["JA-Ile"])
