============
Translations
============

SKM's networks are in majority built on *Arabidopsis thaliana*. Their genes are translated to crop
species by `skm-translate <https://github.com/NIB-SI/skm-translate>`_, a pipeline that
combines several orthology methods. The translation files can be downloaded from the
`SKM downloads page <https://skm.nib.si/downloads>`_, which has the current list of species
and their references.

Species
=======

.. list-table::
   :header-rows: 1

   * - Code
     - Species
     - Common name
   * - ``ath``
     - *Arabidopsis thaliana*
     - thale cress
   * - ``mdo``
     - *Malus domestica*
     - apple
   * - ``pdul``
     - *Prunus amygdalus* syn. *Prunus dulcis*
     - almond
   * - ``parm``
     - *Prunus armeniaca*
     - apricot
   * - ``pavi``
     - *Prunus avium*
     - wild cherry
   * - ``pcer``
     - *Prunus cerasifera*
     - cherry plum
   * - ``ppe``
     - *Prunus persica*
     - peach
   * - ``psib``
     - *Prunus sibirica*
     - siberian apricot
   * - ``pcox``
     - *Pyrus communis*
     - pear
   * - ``sly``
     - *Solanum lycopersicum*
     - tomato
   * - ``stu``
     - *Solanum tuberosum*
     - potato
   * - ``vvi``
     - *Vitis vinifera*
     - grapevine

PSS
===

The genes of PSS functional clusters are listed per species in the ``<species>_homologues``
node attributes, and each species has its own gene network: see :doc:`pss`.

CKN
===

:mod:`skm_tools.translate` translates the Arabidopsis genes of CKN to another species. Each
Arabidopsis gene node is replaced by a node for each of its translations, keeping its edges:

.. code-block:: python

   from skm_tools.translate import load_translation_file, integrate_translation_ckn

   translations = load_translation_file("parm", "translations-parm.tsv")
   ckn_apricot = integrate_translation_ckn(ckn, translations, t_target_col="apricot")
