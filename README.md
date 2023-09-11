# skm-tools
Scripts and utilities for leveraging SKM resources

SKM-tools ( https://github.com/NIB-SI/skm-tools) is a collection of Python scripts and notebooks, incorporating network analysis and visualisation tools, that facilitates interrogation of CKN and PSS with targeted questions beyond the scope of the web application. Included are demonstrations of:

1) Creating a networkX graph object – allowing access to the multitude of graph analysis and graph operation functions available in the networkX library,
2) Network filtering – removing edges or nodes from the network that are not relevant to the biological question at hand,
3) Neighbourhood extraction – identifying the immediate interactors of target entities,
4) Shortest path analysis – identifying directed (or undirected) paths between source and target nodes of interest,
5) Subnetwork extraction – extracting and exporting a subnetwork (that may be the result of previous steps),
6) Inclusion of experimental data – importing experimental data (logFC and p-values values), and adding as context to the networks,
7) Cytoscape Automation18 – loading, styling, highlighting, and exporting results of network analysis performed in Python, in the Cytoscape interface using py4cytoscape1, and
8) Interoperability with the DiNAR application, allowing integration and visualisation of multiple condition high-throughput data in a knowledge graph context.

## Case studies

### [Case study 1: Interaction of ABA, JA, and SA in the activation of RD29 transcription](https://github.com/NIB-SI/skm-tools/tree/main/case-study-1)
Experimental results showed that ABA is able to induce expression of RD29 in both arabidopsis and potato plants. However, the addition of SA or JA showed attenuation of this induction (while alone they had no effect). Here, we interogate the PSS network of stress signalling to identify potential points of intersection between the ABA activation of RD29 and SA or JA pathways, to identify potential mechanistic explainations of the observed experimental results. 

Includes demonstration of: 1, 2, 3, 4, 5, 7

### [Case study 2: The impact of Ca<sup>2+</sup> channel inhibitor LaCl<sub>3</sub> on proteome-wide peroxide responses](https://github.com/NIB-SI/skm-tools/tree/main/case-study-2)
In this case study, we analysed the results of a proteimics study of arabidopsis rosettes treated with either H<sub>2</sub>O<sub>2</sub> or a combination of H<sub>2</sub>O<sub>2</sub> and LaCl<sub>3</sub>. We found 119 proteins that showed significantly changed abundances in response to H<sub>2</sub>O<sub>2</sub> compared to mock treatment after 10 or 30 min of treatment. Out of these, 49 proteins did not significantly change in abundance upon pretreatment with LaCl<sub>3</sub>, indicating that a significant number of H<sub>2</sub>O<sub>2</sub> induced changes in protein abundance required a Ca<sup>2+</sup> signal (Ca<sup>2+</sup>-dependent redox-responsive proteins). Here we analysis CKN to identify potential explainations behind the observed data. 

Includes demonstration of: 1, 2, 3, 4, 5, 6, 7

## Requirements

The only required non-default library is networkX (https://networkx.org/).

Additional optional libraries are
* py4cytoscape (https://py4cytoscape.readthedocs.io/) – necessary if you wish to view the results in Cytoscape.
* graphviz and pygraphviz – necessary for using graphviz layouts instead of builtin networkX or Cytoscape layouts. 
* pypdf, pdfCropMargins, reportlab – necessary if you wish to export multiple network views from Cytoscape to a single file.


## DiNAR
[DiNAR](https://github.com/NIB-SI/DiNAR/) (Differential Network Analysis in R) allows for the examination and visulasation of differential (gene) expression within prior knowledge graphs. 

To use CKN or PSS as the prior knowledge networks in DiNAR, you can follow the steps here:

1) Obtain the (sub)network of interest. This can be done by using extracting a subnetwrk yourself (e.g. see the use cases) or downloading the result of a query in the [PSS](https://skm.nib.si/biomine/) or [CKN](https://skm.nib.si/ckn/) Explorers. Export nodes and edges as `.csv` In the case of PSS, the entire network is also visualisable (download from https://skm.nib.si/downloads/pss/public/dinar-edges). 
<p align="center" width="100%">
  <img src="https://github.com/NIB-SI/DiNAR/blob/master/subApps/GMM-SKM-KnetMiner/figs/SKM-PSS.png" width=45%>
</p>

2) [Download](https://github.com/NIB-SI/DiNAR/blob/master/subApps/GMM-SKM-KnetMiner_customNet/scripts/CustomNetwork_KnetMiner-SKM.R) and run _Custom Network_ [GMM-SKM-KnetMiner](https://github.com/NIB-SI/DiNAR/tree/master/subApps/GMM-SKM-KnetMiner/) app locally or from [shinyapps](https://nib-si.shinyapps.io/GMM-SKM-KnetMiner/)
3) Load [files](https://github.com/NIB-SI/DiNAR/tree/master/subApps/GMM-SKM-KnetMiner/input)
4) Follow instructions and save Custom Network Tables (Nodes/Vertices and Reactions/Edges)
5) Load [Custom Network](https://github.com/NIB-SI/DiNAR/tree/master/subApps/GMM-SKM-KnetMiner/output-for-DiNAR) and Expression Tables into [DiNAR](https://github.com/NIB-SI/DiNAR)
6) Export results as wished ([e.g.](https://github.com/NIB-SI/DiNAR/tree/master/subApps/GMM-SKM-KnetMiner/output-from-DiNAR))

</br>

(*) Take care of empty strings in tables in general, replace them with `-`

(**) Try to keep shortName short


## Publication
### TODO


## Notes and additional resources

### Handy links:
- [Cytosape Automation Wiki](https://github.com/cytoscape/cytoscape-automation/wiki)
- [py4cytoscape.readthedocs.io](https://py4cytoscape.readthedocs.io/en/latest/)

### Network layouts

More information on layout automation:
- [Network layout notebook from cytoscape-automation wiki](https://github.com/cytoscape/cytoscape-automation/blob/master/for-scripters/Python/network-layout.ipynb)

#### Copycat
If the Copycat Layout app is not installed in Cytoscape, you can install it from the App store.

#### yFiles
yFiles does not support Cytoscape Automation.
For small networks, often the yFiles Organic layout is the most optimal. You can set it as you default layout in the layout settings,
and then not apply a layout in the automation.




