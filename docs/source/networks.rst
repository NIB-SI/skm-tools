=================
Loading networks
=================

All loaders return networkx graphs (:mod:`skm_tools.load_networks`).

PSS
===

PSS is exported in three forms, each as an edge file and a node file (tab-separated,
with a header; lists are joined with ``;``):

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
       species (``node_type`` ``gene``).
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
``Complex``, ...), ``display_label``, ``short_name``, ``synonyms``, ``pathway``,
``components`` (for complexes) and ``<species>_homologues``. Empty values are ``None``,
lists are Python lists, and booleans are ``True``/``False``.

.. note::

   The download URLs of the PSS exports on skm.nib.si are not published yet: the PSS
   loaders read local files only.

PSS-specific functions are in :mod:`skm_tools.pss_utils`: filtering by node type or
species (:func:`~skm_tools.pss_utils.filter_pss_nodes`), merging parallel edges
(:func:`~skm_tools.pss_utils.simplify_pss`), and removing nodes while keeping their
upstream and downstream nodes connected (:func:`~skm_tools.pss_utils.remove_and_rewire`).

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
