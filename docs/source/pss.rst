===
PSS
===

The Plant Stress Signalling model (PSS) is a curated, mechanistic model of plant stress
signalling; see the `PSS documentation <https://skm.nib.si/documentation/pss-explore>`_ on
the SKM website for what it contains and how it is built. Everything here is in
:mod:`skm_tools.pss`.

Networks
========

PSS is exported in three network forms, each as an edge file and a node file (tab-separated,
with a header). The exports are made by `skm-pss-export <https://github.com/NIB-SI/skm-pss-export>`_,
which describes the formats and every column in detail.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Network
     - Content
     - Loader
   * - Reaction graph
     - Entities **and reactions** as nodes, one edge per reaction participant
       (participant → reaction for inputs and modifiers, reaction → participant for
       products). The lossless form of PSS, as in the PSS Explorer.
     - :func:`~skm_tools.pss.pss_reaction_graph_to_networkx`
   * - Interaction network
     - Entities as nodes; edges are entity → entity influences through the reactions.
       Species-independent: plant nodes are functional clusters.
     - :func:`~skm_tools.pss.pss_interaction_network_to_networkx`
   * - Gene network (per species)
     - The interaction network for one species, with the functional clusters replaced by
       their genes in that species.
     - :func:`~skm_tools.pss.pss_gene_network_to_networkx`

Use the **interaction network** for species-independent analyses, a **gene network** to
combine PSS with CKN or with experimental data of one species, and the **reaction graph**
when the reactions themselves matter (their participants and roles).

Loading
=======

All three are loaded as :class:`networkx.MultiDiGraph`, as several reactions (or roles)
can link the same two nodes. Missing files are downloaded from
`skm.nib.si <https://skm.nib.si>`_: to the given paths, or, without paths, into ``data_dir``
(default: the current folder) under the names the browser gives them, without the export
date (e.g. ``pss-public-interaction-network-edges.tsv``):

.. code-block:: python

   from skm_tools.pss import pss_interaction_network_to_networkx, pss_gene_network_to_networkx

   pss = pss_interaction_network_to_networkx(data_dir="data")
   pss.nodes["WRKY33[fc00166]"]["display_label"]   # 'WRKY33'

   pss_stu = pss_gene_network_to_networkx(species="stu", data_dir="data")   # default species "ath"

   # or with your own file names
   pss = pss_interaction_network_to_networkx("pss-edges.tsv", "pss-nodes.tsv")

The loaders record which network a graph is in ``g.graph["pss_network"]``
(``"reaction_graph"``, ``"interaction_network"`` or ``"gene_network"``, and the species of a
gene network in ``g.graph["species"]``), so the functions below can check that they are
used on a network they work on.

In the files, lists are joined with ``|`` (names can contain ``,``, e.g. ``AHK2,3,4``, and
gene symbols ``;``, e.g. ``PIP1;3``), and the loaders return them as Python lists. Empty
values are ``None``, and booleans are ``True``/``False``.

Nodes
-----

Every node has ``node_type`` (the PSS class: ``PlantCoding``, ``Metabolite``, ``Complex``,
..., see the `PSS database schema <https://skm.nib.si/documentation/pss-db#node-labels>`_),
``display_label`` (the name to show), ``short_name``, ``synonyms``, ``pathway`` and
``mapman`` (see :ref:`mapman`). Complexes have ``components`` and ``component_cluster_ids``.

**Interaction network and reaction graph:** the plant nodes (``PlantCoding``,
``PlantNonCoding``) are `functional clusters <https://skm.nib.si/documentation/pss-explore#FC>`_,
named ``short_name[functional_cluster_id]`` (e.g. ``WRKY33[fc00166]``), with their genes per
species in ``<species>_homologues`` attributes (species codes: :doc:`translations`).

**Gene network:** there are no functional clusters; their genes are nodes instead (gene ids,
with ``species`` set, and the ``node_type`` of their cluster: ``PlantCoding`` or
``PlantNonCoding``), and only the reactions whose functional clusters all have genes in
the species are included. A gene keeps the clusters it comes from: ``short_name``,
``pathway`` and ``functional_cluster_id`` are lists (one entry per cluster, in the same
order, as a gene can be in several clusters); ``display_label`` is a single string. Other
nodes (metabolites, complexes, ...) are as in the interaction network.

**Reaction graph:** reactions are nodes too (``node_type`` ``reaction``, the reaction id as
node id), with ``reaction_type``, ``reaction_effect`` and ``reaction_mechanism``.

Edges
-----

**Interaction network and gene network:** edges are keyed by ``reaction_id``, with
``interaction`` (``positive-influence``, ``negative-influence`` or ``unknown-influence``),
the reaction's ``reaction_type`` and ``reaction_effect``, and the participants' roles
(``source_role``, ``target_role``). In a gene network, ``source_entity`` and
``target_entity`` give the functional cluster a gene acted through. Mutual influences
(``directed`` False, e.g. between binding partners) are listed in both directions. Every
edge has ``rank`` 0, CKN's rank for curated PSS interactions, so PSS and CKN can be combined
and filtered alike.

**Reaction graph:** edges are keyed by ``role`` (e.g. ``substrate``, ``catalyst``,
``product``), as an entity can take part in the same reaction twice (e.g. as template and
stimulator).

Filtering
=========

The filtering functions change the graph in place (make a copy first, ``g.copy()``, to
keep the original), and return the reasons nodes were removed (``{node: reason}``). They
work on whole reactions where it matters, so the networks stay consistent: in the reaction
graph a reaction is a node, and in the interaction network and the gene networks it is the
edges with its ``reaction_id``. Afterwards, nodes without edges are removed too (also those
that had none before; ``remove_isolates=False`` keeps them). Arguments taking a list also
take a single value, e.g. ``species="stu"``.

:func:`~skm_tools.pss.remove_reactions` removes reactions by id, in all three networks.

:func:`~skm_tools.pss.filter_pss_nodes` keeps only the reactions of some species (interaction
network and reaction graph): it removes the functional clusters without genes in the
species, the complexes they are components of, and their reactions (e.g. for a catalysis
whose enzyme has no genes in the species, also the substrate → product edge). Abstract
clusters (``PlantAbstract``, without genes) are kept, as metabolites. The result has the
same reactions as the species' gene network, with functional clusters instead of genes.
It can also keep only some node types (interaction network and gene networks; in the
reaction graph, this would leave reactions with missing participants); complexes that lose
a component are removed too (in a gene network, once all genes of a component cluster are
gone):

.. code-block:: python

   from skm_tools.pss import filter_pss_nodes, remove_reactions

   reasons = filter_pss_nodes(pss, species=["stu"])
   reasons = filter_pss_nodes(pss, node_types=["PlantCoding", "Complex", "Metabolite"])
   reasons = filter_pss_nodes(pss_stu, node_types=["PlantCoding", "Complex", "Metabolite"])
   reasons = remove_reactions(pss, ["rx00001"])

:func:`~skm_tools.pss.remove_deadend_complexes` removes complexes that influence nothing
(not an input or modifier of any reaction), which are dead ends in directed analyses. In
the reaction graph, the reactions forming them are removed too; in the interaction network
and the gene networks, the binding partners' mutual edges from those reactions are kept.
The two can give different results: in the interaction network, the substrate of a
degradation has no outgoing edge, so a complex that is only degraded is a dead end there,
but not in the reaction graph.

Simplifying
===========

These work on the interaction network and the gene networks; they raise a ``ValueError``
for the reaction graph.

:func:`~skm_tools.pss.simplify_pss` merges the parallel edges (one per reaction) between two
nodes into one edge, and returns a new :class:`networkx.DiGraph`. The merged edge keeps all
``reaction_id`` values (a sorted list; every edge of the result has a list), and one value of
every other attribute, by a fixed rule: ``interaction`` becomes ``unknown-influence`` if the
merged edges disagree (e.g. a positive and a negative influence), ``directed`` is True if
any edge is directed, ``rank`` the lowest, and anything else the value of the first
reaction. How many merged edges had differing values is logged (``verbose=True`` lists
them). To keep such edges apart, pass ``split_on_attrs=["interaction"]``; the result is then
a :class:`networkx.MultiDiGraph`.

:func:`~skm_tools.pss.remove_and_rewire` removes nodes while connecting each of their
upstream nodes to each of their downstream nodes, e.g. to hide intermediate steps. It needs
a simple :class:`networkx.DiGraph`, so run :func:`~skm_tools.pss.simplify_pss` first. The
nodes are removed one at a time, so chains of removed nodes are bridged. A new edge A → C
(replacing A → B → C) has the sign of the chain (e.g. an activation followed by an
inhibition: ``negative-influence``) and the reaction ids of both edges; if A → C exists
already, it is merged into it. Rewiring doesn't go through the mutual edges between binding
partners, only through complex formation (partner → complex). ``dry_run=True`` returns the
edges that would be added, without changing the graph.

:func:`~skm_tools.pss.remove_duplicated_binding_edges` keeps one direction of each mutual
binding edge, for a less cluttered drawing. Don't use it before path or neighbourhood analysis.

.. code-block:: python

   from skm_tools.pss import simplify_pss, remove_and_rewire

   simple = simplify_pss(pss)
   planned = remove_and_rewire(simple, ["JA-Ile"], dry_run=True)   # [(u, v, data), ...]
   reasons = remove_and_rewire(simple, ["JA-Ile"])
