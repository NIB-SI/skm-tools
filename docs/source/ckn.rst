===
CKN
===

The Comprehensive Knowledge Network (CKN) is a network of molecular interactions in
*Arabidopsis thaliana*, mostly from high-throughput experiments. Everything here is in
:mod:`skm_tools.ckn`.

Loading
=======

:func:`~skm_tools.ckn.ckn_to_networkx` loads CKN as a :class:`networkx.DiGraph`,
downloading the files from `skm.nib.si <https://skm.nib.si>`_ if they don't exist yet:

.. code-block:: python

   from skm_tools.ckn import ckn_to_networkx

   ckn = ckn_to_networkx("ckn-edges.tsv.gz", "ckn-nodes.tsv.gz")

By default the reverse of every undirected edge is added, so directed path searches can
use undirected edges in both directions (``add_reciprocal_edges=True``); or keep only the
directed edges with ``directed=True``.

Node attributes include ``node_type`` (``protein_coding``, ``mirna``,
``transposable_element_gene``, ``metabolite``, ``complex``, ...), ``species`` (``ath``,
``foreign`` or ``ath/foreign``; missing for metabolites), ``TAIR``, ``short_name``,
``synonyms``, ``full_name``, ``GMM`` (MapMan bins) and ``tissue`` (``leaf``, ``flower``,
``stem``, ``root``, ``seed`` or ``not assigned``); ``synonyms``, ``GMM`` and ``tissue`` are
lists. Edge attributes include ``type`` (``binding``, ``small RNA interactions``,
``transcription factor regulation``, ``post-translational modification``, ``other``),
``effect`` and ``rank``.

Filtering
=========

The filtering functions change the graph in place (make a copy first, ``g.copy()``, to keep
the original), and remove nodes left without edges.

:func:`~skm_tools.ckn.filter_ckn_edges` keeps only edges of some ranks (0, the best
supported, to 4) and/or types, or those for which your own function returns True.
:func:`~skm_tools.ckn.rank_counts` counts the edges of each rank.

.. code-block:: python

   from skm_tools.ckn import filter_ckn_edges, filter_ckn_nodes

   filter_ckn_edges(ckn, keep_edge_ranks=[0, 1, 2])
   filter_ckn_edges(ckn, keep_edge_types=["binding", "transcription factor regulation"])
   filter_ckn_edges(ckn, filter_function=lambda d: d["isDirected"] == 1)

:func:`~skm_tools.ckn.filter_ckn_nodes` keeps only nodes of some types, species or tissues,
and returns the reasons nodes were removed. Complexes that lose a component are removed too.
Nodes without a species (metabolites) are kept by the species filter.

.. code-block:: python

   reasons = filter_ckn_nodes(ckn, node_types=["protein_coding", "metabolite"], tissues=["leaf"])

Annotations
===========

:func:`~skm_tools.ckn.get_nodes_by_annotation` finds the nodes with any of the given MapMan
(GMM) bins, by default including their sub-bins; :func:`~skm_tools.ckn.get_all_annotations`
lists all values of an annotation.

.. code-block:: python

   from skm_tools.ckn import get_nodes_by_annotation

   tfs = get_nodes_by_annotation(ckn, gmm=["27.3"])   # 27.3: regulation of transcription

Translation to other species
============================

:mod:`skm_tools.translate` translates the Arabidopsis genes of CKN to another species,
using the `skm-translate <https://github.com/NIB-SI/skm-translate>`_ translation files from the
`SKM downloads page <https://skm.nib.si/downloads>`_ (see :doc:`pss` for the species codes):

.. code-block:: python

   from skm_tools.translate import load_translation_file, integrate_translation_ckn

   translations = load_translation_file("parm", "translations-parm.tsv")
   ckn_apricot = integrate_translation_ckn(ckn, translations, t_target_col="apricot")
