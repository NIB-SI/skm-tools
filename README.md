# skm-tools
Scripts and utilities for leveraging SKM resources

Stress Knowledge Map (SKM) is a knowledge graph resource for systems biology analysis of plant stress responses. SKM contains knowledge on molecular interactions in plants which have been integrated from a wide diversity of sources into a single, freely available entrypoint. The two knowledge graphs in SKM are: 
* The Plant Stress Signalling model (PSS) -- a highly curated and detailed mechanistic plant stress signalling model, in majority compiled from targeted biochemical studies available in literature and containing over 500 interactions,
* The Comprehensive Knowledge Network (CKN) -- a compilation of molecular interactions in the plant cell, in majority is based on high-throughput experiments, and contains over 26.000 molecules and 480.000 interactions.

SKM is avalible at https://skm.nib.si, and provides tools for exploring both PSS and CKN. 

SKM-tools (https://github.com/NIB-SI/skm-tools) is a collection of Python scripts and notebooks, incorporating network analysis and visualisation tools, that facilitates interrogation of CKN and PSS with targeted questions beyond the scope of the web application. Included are demonstrations of:

1) Creating a networkX graph object – allowing access to the multitude of graph analysis and graph operation functions available in the networkX library,
2) Network filtering – removing edges or nodes from the network that are not relevant to the biological question at hand,
3) Neighbourhood extraction – identifying the immediate interactors of target entities,
4) Shortest path analysis – identifying directed (or undirected) paths between source and target nodes of interest,
5) Subnetwork extraction – extracting and exporting a subnetwork (that may be the result of previous steps),
6) Inclusion of experimental data – importing experimental data (logFC and p-values values), and adding as context to the networks,
7) Cytoscape Automation – loading, styling, highlighting, and exporting results of network analysis performed in Python, in the Cytoscape interface using py4cytoscape, and
8) Interoperability with the DiNAR application, allowing integration and visualisation of multiple condition high-throughput data in a knowledge graph context.

## DiNAR & SKM
Instructions for using [DiNAR](https://github.com/NIB-SI/DiNAR/) (Differential Network Analysis in R) with SKM are here: [DiNAR & SKM](https://github.com/NIB-SI/skm-tools/tree/main/DiNAR.md)

## Case studies from the publication

See: https://doi.org/10.1016/j.xplc.2024.100920.

The `publication` branch contains skm-tools as it existed when performing the analysis in the publication.

### [Case study 1: Interaction of ABA, JA, and SA in the activation of RD29 transcription](https://github.com/NIB-SI/skm-tools/tree/main/publication/case-study-1)
Experimental results showed that ABA is able to induce expression of RD29 in both arabidopsis and potato plants. However, the addition of SA or JA showed attenuation of this induction (while alone they had no effect). Here, we interrogate the PSS network of stress signalling to identify potential points of intersection between the ABA activation of RD29 and SA or JA pathways, to identify potential mechanistic explanations of the observed experimental results.

Includes demonstration of: 1, 2, 3, 4, 5, 7

### [Case study 2: The impact of Ca<sup>2+</sup> channel inhibitor LaCl<sub>3</sub> on proteome-wide peroxide responses](https://github.com/NIB-SI/skm-tools/tree/main/publication/case-study-2)
In this case study, we analysed the results of a proteimics study of arabidopsis rosettes treated with either H<sub>2</sub>O<sub>2</sub> or a combination of H<sub>2</sub>O<sub>2</sub> and LaCl<sub>3</sub>. We found 119 proteins that showed significantly changed abundances in response to H<sub>2</sub>O<sub>2</sub> compared to mock treatment after 10 or 30 min of treatment. Out of these, 49 proteins did not significantly change in abundance upon pretreatment with LaCl<sub>3</sub>, indicating that a significant number of H<sub>2</sub>O<sub>2</sub> induced changes in protein abundance required a Ca<sup>2+</sup> signal (Ca<sup>2+</sup>-dependent redox-responsive proteins). Here we analysis CKN to identify potential explanations behind the observed data.

Includes demonstration of: 1, 2, 3, 4, 5, 6, 7

## Install and requirements

The only required non-default library is networkX (https://networkx.org/).

Additional __optional__ libraries are
* py4cytoscape (https://py4cytoscape.readthedocs.io/) – necessary if you wish to view the results interactively in Cytoscape.
* graphviz and pygraphviz – necessary for using graphviz layouts instead of networkX or Cytoscape layouts.
* pypdf, pdfCropMargins, reportlab – necessary if you wish to batch export multiple network views from Cytoscape to a single file.

### Install

To install from GitHub:

    pip install "skm-tools @ git+ssh://git@github.com/NIB-SI/skm-tools.git"

To install with Cytoscape and PDF features enabled:

    pip install "skm-tools[cytoscape,pdf] @ git+ssh://git@github.com/NIB-SI/skm-tools.git"

## Citation

If you used skm-tools (or SKM) in your work, please cite: 

> Carissa Bleker, Živa Ramšak, Andras Bittner, Vid Podpečan, Maja Zagorščak, Bernhard Wurzinger, Špela Baebler, Marko Petek, Maja Križnik, Annelotte van Dieren, Juliane Gruber, Leila Afjehi-Sadat, Wolfram Weckwert, Anže Županič, Markus Teige, Ute C. Vothknecht, Kristina Gruden. (2024). Stress Knowledge Map: A knowledge graph resource for systems biology analysis of plant stress responses. Plant Communications. doi:10.1016/j.xplc.2024.100920 .

## Notes and additional resources

### Handy links:
- [Cytosape Automation Wiki](https://github.com/cytoscape/cytoscape-automation/wiki)
- [py4cytoscape.readthedocs.io](https://py4cytoscape.readthedocs.io/en/latest/)

### Network layouts

More information on layout automation:
- [Network layout notebook from Cytoscape-automation wiki](https://github.com/cytoscape/cytoscape-automation/blob/master/for-scripters/Python/network-layout.ipynb)

#### Copycat
If the Copycat Layout app is not installed in Cytoscape, you can install it from the App store.

#### yFiles
yFiles does not support Cytoscape Automation.
For small networks, often the yFiles Organic layout is the most optimal. You can set it as you default layout in the layout settings,
and then not apply a layout in the automation.

### Pygraphviz on windows

If you have issues installing pygraphviz on windows, the following link may have helpful suggestions: [stackoverflow|howto install pygraphviz on windows 10 64bit](https://stackoverflow.com/questions/40809758/howto-install-pygraphviz-on-windows-10-64bit).

Please note that any issues installing dependencies should be directed to the developer of the dependency in question.
