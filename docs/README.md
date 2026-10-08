# skm-tools

[![Python >=3.10](https://img.shields.io/badge/python-%3E%3D3.10-blue)](https://github.com/NIB-SI/skm-tools)
[![License: GPL-3.0](https://img.shields.io/github/license/NIB-SI/skm-tools)](https://github.com/NIB-SI/skm-tools/blob/main/LICENSE)
[![Docs](https://img.shields.io/github/actions/workflow/status/NIB-SI/skm-tools/deploy-docs.yml?label=docs)](https://nib-si.github.io/skm-tools)
[![Tests](https://img.shields.io/github/actions/workflow/status/NIB-SI/skm-tools/tests.yml?label=tests)](https://github.com/NIB-SI/skm-tools/actions/workflows/tests.yml)

Python tools for systems biology analysis of plant stress responses with the
[Stress Knowledge Map](https://skm.nib.si) (SKM).

SKM contains knowledge on molecular interactions in plants, integrated from a wide diversity of sources into a single,
freely available entrypoint. It has two knowledge graphs:

* **PSS** (Plant Stress Signalling): a highly curated, detailed mechanistic model of plant stress signalling, mostly
  compiled from targeted biochemical studies in the literature, with over 800 reactions.
* **CKN** (Comprehensive Knowledge Network): molecular interactions in the plant cell, mostly from high-throughput
  experiments, with over 26,000 molecules and 480,000 interactions.

skm-tools loads both as [networkx](https://networkx.org) graphs, and adds functions for the analyses SKM is typically
used for: filtering, shortest paths, neighbourhoods, minimum cuts, experimental data, translation to other species,
Cytoscape automation and DiNAR export.

## CI status

| Tests | Python 3.10 | Python 3.12 | Python 3.13 |
|---|---|---|---|
| Ubuntu | [![3.10](https://img.shields.io/github/check-runs/NIB-SI/skm-tools/main?nameFilter=test%20%283.10%29)](https://github.com/NIB-SI/skm-tools/actions/workflows/tests.yml) | [![3.12](https://img.shields.io/github/check-runs/NIB-SI/skm-tools/main?nameFilter=test%20%283.12%29)](https://github.com/NIB-SI/skm-tools/actions/workflows/tests.yml) | [![3.13](https://img.shields.io/github/check-runs/NIB-SI/skm-tools/main?nameFilter=test%20%283.13%29)](https://github.com/NIB-SI/skm-tools/actions/workflows/tests.yml) |

| Tutorials | PSS | CKN | Heatmaps |
|---|---|---|---|
| Ubuntu, Python 3.12 | [![PSS](https://img.shields.io/github/check-runs/NIB-SI/skm-tools/main?nameFilter=tutorials%20%28tutorial-pss%29)](https://github.com/NIB-SI/skm-tools/actions/workflows/tests.yml) | [![CKN](https://img.shields.io/github/check-runs/NIB-SI/skm-tools/main?nameFilter=tutorials%20%28tutorial-ckn%29)](https://github.com/NIB-SI/skm-tools/actions/workflows/tests.yml) | [![Heatmaps](https://img.shields.io/github/check-runs/NIB-SI/skm-tools/main?nameFilter=tutorials%20%28tutorial-heatmaps%29)](https://github.com/NIB-SI/skm-tools/actions/workflows/tests.yml) |

## Documentation

[nib-si.github.io/skm-tools](https://nib-si.github.io/skm-tools): installation, the PSS and CKN networks, the analysis
functions, Cytoscape automation, and tutorials. For using [DiNAR](https://github.com/NIB-SI/DiNAR/) (Differential
Network Analysis in R) with SKM, also without code, see [DiNAR](https://nib-si.github.io/skm-tools/dinar.html).

## Installation

```bash
pip install "skm-tools @ git+https://github.com/NIB-SI/skm-tools.git"
```

With the optional extras (Cytoscape automation, PDF export), see the
[installation page](https://nib-si.github.io/skm-tools/installation.html).

## Case studies from the publication

The [`publication`](https://github.com/NIB-SI/skm-tools/tree/main/publication) folder has the case studies of the
[SKM publication](https://doi.org/10.1016/j.xplc.2024.100920), and the `publication` branch contains skm-tools as it
was when the analyses for the publication were made. The [tutorials](https://nib-si.github.io/skm-tools) follow the
case studies with the current skm-tools:

* Case study 1: interaction of ABA, JA and SA in the activation of RD29 transcription (PSS),
* Case study 2: the impact of the Ca<sup>2+</sup> channel inhibitor LaCl<sub>3</sub> on proteome-wide peroxide
  responses (CKN).

## Citation

If you use skm-tools (or SKM) in your work, please cite:

> Carissa Bleker, Živa Ramšak, Andras Bittner, Vid Podpečan, Maja Zagorščak, Bernhard Wurzinger, Špela Baebler, Marko
> Petek, Maja Križnik, Annelotte van Dieren, Juliane Gruber, Leila Afjehi-Sadat, Wolfram Weckwerth, Anže Županič,
> Markus Teige, Ute C. Vothknecht, Kristina Gruden. (2024). Stress Knowledge Map: A knowledge graph resource for
> systems biology analysis of plant stress responses. Plant Communications. doi:10.1016/j.xplc.2024.100920

## License

[GPL-3.0](https://github.com/NIB-SI/skm-tools/blob/main/LICENSE)
