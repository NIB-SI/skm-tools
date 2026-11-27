# Changelog

Versions are tagged `v<version>`. The `publication` branch is skm-tools as used for the SKM publication
(doi:10.1016/j.xplc.2024.100920).

## 0.2.0 (2026-11-27)

A rewrite for the network exports of the SKM web app (skm.nib.si, from 2026-11-27): the PSS downloads made by
[pss-export](https://github.com/NIB-SI/skm-pss-export) and CKN v2.1. Not compatible with 0.1 (the old modules and
the old PSS and CKN files are gone); the case studies of the publication still need the `publication` branch.

- Loading: the PSS reaction graph, interaction network and gene networks (one per species) and CKN v2.1 as
  networkx graphs, downloaded from skm.nib.si when missing (`skm_tools.pss`, `skm_tools.ckn`); lists, booleans and
  integers typed, empty values `None`.
- PSS: filtering by node type and species and removing reactions, by whole reactions (`filter_pss_nodes`,
  `remove_reactions`, `remove_deadend_complexes`); `simplify_pss`, `remove_and_rewire`,
  `remove_duplicated_binding_edges`; `metabolite_enzymes`, `pathway_gene_sets`.
- CKN: `filter_ckn_edges`, `filter_ckn_nodes`, `rank_counts`, `transcription_factors`, `tf_target_enrichment`.
- Analysis: shortest paths, neighbourhoods, minimum cuts, connectivity of a node set against degree-matched random
  nodes (`induced_connectivity`), hypergeometric test and BH FDR without scipy (`skm_tools.stats`).
- Annotations: MapMan lookups, MapMan gene sets from Mercator files, PlantGSEA gene sets.
- Experimental data and any table columns as node attributes (`add_experimental_data`, `add_node_data`).
- Translation of CKN's Arabidopsis genes to other species with the SKM translation files.
- Cytoscape automation (py4cytoscape, `skm-tools[cytoscape]`): the SKM styles, flag mappings, highlighting,
  subnetworks, node images and charts (e.g. heatmaps), export; heatmap images without Cytoscape.
- DiNAR custom networks (`skm_tools.dinar`); saving graphs (pickle, JSON) and JSON-safe dicts of graphs and paths;
  PDF cropping and combining (`skm-tools[pdf]`).
- Documentation (https://nib-si.github.io/skm-tools) with tutorials (PSS, CKN, heatmaps) that run in CI.
- Packaging: `pyproject.toml` (setuptools); `networkx` and `pandas` the only required dependencies; Python 3.10+.

## 0.1 (2024-10-16)

The scripts and notebooks of the SKM publication.
