===
CKN
===

The Comprehensive Knowledge Network (CKN) is a network of molecular interactions in
*Arabidopsis thaliana*, mostly from high-throughput experiments; see the
`CKN documentation <https://skm.nib.si/documentation/ckn-explore>`_ on the SKM website for its
sources and interaction types. Everything here is in :mod:`skm_tools.ckn`.

Loading
=======

:func:`~skm_tools.ckn.ckn_to_networkx` loads CKN as a :class:`networkx.DiGraph`,
downloading the files from `skm.nib.si <https://skm.nib.si>`_ if they don't exist locally:

.. code-block:: python

   from skm_tools.ckn import ckn_to_networkx

   ckn = ckn_to_networkx("ckn-edges.tsv.gz", "ckn-nodes.tsv.gz")

By default the reverse of every undirected edge is added, so directed path searches can
use undirected edges in both directions (``add_reciprocal_edges=True``); or keep only the
directed edges with ``directed=True``.

Node attributes include ``node_type``, ``species``, ``TAIR``, ``display_label``,
``short_name``, ``synonyms``, ``description``, ``mapman`` (see :ref:`mapman`) and ``tissue``.
Edge attributes include ``interaction`` (``positive-influence``, ``negative-influence`` or
``unknown-influence``, as in PSS), ``directed``, ``rank``, ``effect``, ``type`` and
``interactionSources``. Lists (e.g. ``synonyms``, ``tissue``) are joined with ``|`` in the
files and loaded as Python lists; empty values are ``None``. Files in the older CKN v2 format
are converted to these attributes when loaded.

Filtering
=========

The filtering functions change the graph in place (make a copy first, ``g.copy()``, to keep
the original), and remove nodes left without edges.

:func:`~skm_tools.ckn.filter_ckn_edges` keeps only edges of some
`ranks <https://skm.nib.si/documentation/ckn-explore#ranks>`_ (from 0, curated in PSS, to 4,
purely predicted) and/or types, or those for which your own function returns True.
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

MapMan annotations
==================

To find nodes by MapMan bin, in CKN or PSS, see :ref:`mapman`.

Translation to other species
============================

To translate CKN's Arabidopsis genes to another species, see :doc:`translations`.
