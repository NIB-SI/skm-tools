===
CKN
===

The Comprehensive Knowledge Network (CKN) is a network of molecular interactions in
*Arabidopsis thaliana*, mostly from high-throughput experiments; see the
`CKN documentation <https://skm.nib.si/documentation/ckn-explore>`_ on the SKM website for its
sources and interaction types. Everything here is in :mod:`skm_tools.ckn`.

Loading
=======

:func:`~skm_tools.ckn.ckn_to_networkx` loads CKN (v2.0.1) as a :class:`networkx.DiGraph`,
downloading the files from `skm.nib.si <https://skm.nib.si>`_ if they don't exist locally:
to the given paths, or, without paths, into ``data_dir`` (default: the current folder) under
the names the browser gives them (``AtCKN-v2.0.1-2026.10.tsv.gz`` and
``AtCKN-v2.0.1-2026.10_node-annot.tsv.gz``). A file name ending in ``.gz`` is saved gzipped,
any other name uncompressed.

.. code-block:: python

   from skm_tools.ckn import ckn_to_networkx

   ckn = ckn_to_networkx(data_dir="data")

By default, undirected edges (e.g. binding) are added in both directions, so directed path
searches can use them either way (``as_directed=True``); with ``as_directed=False``, every
edge is added once, as listed in the file. Nodes of the node file without edges are left
out. Older CKN versions (before v2.0.1) are not supported.

Node attributes include ``node_type``, ``locus_type``, ``species``, ``TAIR``,
``display_label``, ``short_name``, ``synonyms``, ``description``, ``mapman`` (see
:ref:`mapman`) and ``tissue``. ``node_type`` is the PSS class, as in PSS (``PlantCoding``,
``PlantNonCoding``, ``Metabolite``, ``Complex``, ``ForeignCoding``, ..., see the
`PSS database schema <https://skm.nib.si/documentation/pss-db#node-labels>`_), plus
``PlantPseudogene`` for pseudogenes. For genes and RNAs, ``locus_type`` has the TAIR locus
type (e.g. ``protein_coding``, ``mirna``, ``transposable_element_gene``).
Edge attributes include ``interaction`` (``positive-influence``, ``negative-influence`` or
``unknown-influence``, as in PSS), ``directed``, ``rank``, ``effect``, ``type`` and
``interactionSources``. Lists (e.g. ``synonyms``, ``tissue``) are joined with ``|`` in the
files and loaded as Python lists; empty values are ``None``. The files are read as written,
without quoting (a ``"`` is part of a value), as the PSS exports.

Filtering
=========

The filtering functions change the graph in place (make a copy first, ``g.copy()``, to keep
the original). Afterwards, nodes without edges are removed (also those that had none before;
``remove_isolates=False`` keeps them). Arguments taking a list also take a single value or
any other collection, e.g. ``keep_edge_ranks=0`` or ``species="ath"``.

:func:`~skm_tools.ckn.filter_ckn_edges` keeps only edges of some
`ranks <https://skm.nib.si/documentation/ckn-explore#ranks>`_ (from 0, curated in PSS, to 4,
purely predicted) and/or types, or those for which your own function returns True. It
returns the removed edges and the removed nodes (``{node: "isolate"}``).
:func:`~skm_tools.ckn.rank_counts` counts the edges of each rank.

.. code-block:: python

   from skm_tools.ckn import filter_ckn_edges, filter_ckn_nodes

   removed_edges, removed_nodes = filter_ckn_edges(ckn, keep_edge_ranks=[0, 1, 2])
   filter_ckn_edges(ckn, keep_edge_types=["binding", "transcription factor regulation"])
   filter_ckn_edges(ckn, filter_function=lambda d: d["directed"])   # only directed edges

:func:`~skm_tools.ckn.filter_ckn_nodes` keeps only nodes of some types, species or tissues,
and returns the reasons nodes were removed (``{node: reason}``). Complexes that lose a
component are removed too: a complex can't exist without its components. So complexes of a
plant protein and a pathogen protein (species ``ath/foreign``, e.g. ``RISC|virus_vsiRNA``)
are removed by ``species=["ath"]``. Nodes without a species (metabolites) are kept by the
species filter.

.. code-block:: python

   reasons = filter_ckn_nodes(ckn, node_types=["PlantCoding", "Metabolite"], tissues=["leaf"])

MapMan annotations
==================

To find nodes by MapMan bin, in CKN or PSS, see :ref:`mapman`.

Translation to other species
============================

To translate CKN's Arabidopsis genes to another species, see :doc:`translations`.
