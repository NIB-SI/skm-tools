=================
Loading networks
=================

All loaders return networkx graphs (:mod:`skm_tools.load_networks`).

PSS
===

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
     - :func:`~skm_tools.load_networks.pss_reaction_graph_to_networkx`
   * - Interaction network
     - Entities as nodes; edges are entity → entity influences through the reactions
       (``interaction``: ``positive-influence``, ``negative-influence`` or ``unknown-influence``).
       Edges are keyed by ``reaction_id``.
     - :func:`~skm_tools.load_networks.pss_interaction_network_to_networkx`
   * - Gene network (per species)
     - The interaction network with functional clusters expanded into the genes of one
       species (``node_type`` ``gene``), e.g. ``ath``, ``stu`` or ``mdo`` (see `Species`_).
     - :func:`~skm_tools.load_networks.pss_gene_network_to_networkx`

All three are :class:`networkx.MultiDiGraph`, as several reactions can link the same two
nodes. Mutual influences (``directed`` False, e.g. between binding partners) are listed in
both directions.

.. code-block:: python

   from skm_tools.load_networks import pss_interaction_network_to_networkx

   pss = pss_interaction_network_to_networkx(
       "pss-interaction-network-edges-public.tsv",
       "pss-interaction-network-nodes-public.tsv",
   )
   pss.nodes["WRKY33[fc00166]"]["display_label"]   # 'WRKY33'

Node attributes include ``node_type`` (the PSS class: ``PlantCoding``, ``Metabolite``,
``Complex``, ..., see the `PSS database schema <https://skm.nib.si/documentation/pss-db#node-labels>`_),
``display_label``, ``short_name``, ``synonyms``, ``pathway``, ``components`` (for complexes)
and ``<species>_homologues``. In the files, lists are joined with ``;`` (names can contain
commas); the loaders return them as Python lists. Empty values are ``None``, and booleans
are ``True``/``False``.

PSS-specific functions are in :mod:`skm_tools.pss_utils`: filtering by node type or
species (:func:`~skm_tools.pss_utils.filter_pss_nodes`), merging parallel edges
(:func:`~skm_tools.pss_utils.simplify_pss`), and removing nodes while keeping their
upstream and downstream nodes connected (:func:`~skm_tools.pss_utils.remove_and_rewire`).

Species
-------

PSS is species-independent: its plant nodes are functional clusters, whose genes are listed
per species in the ``<species>_homologues`` node attributes of the reaction graph and the
interaction network. A gene network expands the clusters into the genes of one species, and
only has the reactions whose functional clusters all have genes in that species.

.. list-table::
   :header-rows: 1

   * - Code
     - Species
   * - ``ath``
     - *Arabidopsis thaliana*
   * - ``stu``
     - potato (*Solanum tuberosum*)
   * - ``sly``
     - tomato (*Solanum lycopersicum*)
   * - ``mdo``
     - apple (*Malus domestica*)
   * - ``vvi``
     - grapevine (*Vitis vinifera*)
   * - ``ppe``
     - peach (*Prunus persica*)
   * - ``pavi``
     - sweet cherry (*Prunus avium*)
   * - ``pcer``
     - sour cherry (*Prunus cerasus*)
   * - ``pdul``
     - almond (*Prunus dulcis*)
   * - ``parm``
     - apricot (*Prunus armeniaca*)
   * - ``pcox``
     - pear (*Pyrus communis*)
   * - ``psib``
     - Siberian apricot (*Prunus sibirica*)

The genes of the crop species are translated from Arabidopsis by
`skm-translate <https://github.com/NIB-SI/skm-translate>`_, which combines several orthology
methods. The same translations can be used to translate CKN with :mod:`skm_tools.translate`.

CKN
===

:func:`~skm_tools.load_networks.ckn_to_networkx` loads CKN as a :class:`networkx.DiGraph`,
downloading the files from skm.nib.si if they don't exist yet:

.. code-block:: python

   from skm_tools.load_networks import ckn_to_networkx

   ckn = ckn_to_networkx("ckn-edges.tsv.gz", "ckn-nodes.tsv.gz")

By default the reverse of every undirected edge is added, so directed path searches can
use undirected edges in both directions (``add_reciprocal_edges=True``); or keep only the
directed edges with ``directed=True``.

CKN-specific functions are in :mod:`skm_tools.ckn_utils`: filtering edges by rank or
type (:func:`~skm_tools.ckn_utils.filter_ckn_edges`), nodes by type, species or tissue
(:func:`~skm_tools.ckn_utils.filter_ckn_nodes`), and finding nodes by MapMan (GMM)
annotation (:func:`~skm_tools.ckn_utils.get_nodes_by_annotation`). To translate CKN to
another species' genes, see :mod:`skm_tools.translate`.

Filter functions change the graph in place and return the reasons nodes were removed;
make a copy first (``g.copy()``) to keep the original.
