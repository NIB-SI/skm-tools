# skm-tools
Scripts and utilities for leveraging SKM resources

SKM-tools ( https://github.com/NIB-SI/skm-tools) is a collection of Python scripts and notebooks, incorporating network analysis and visualisation tools, that facilitates interrogation of  CKN and PSS with targeted questions beyond the scope of the web application. Included are demonstrations of:

* Creating a networkX1 graph object – allowing access to the multitude of graph analysis and graph operation functions available in the networkX library,
* Network edge filtering – removing edges or nodes from the network that are not relevant to the biological question at hand,
* Neighbourhood extraction – identifying the immediate interactors of target entities,
* Shortest path analysis – identifying directed (or undirected) paths between source and target nodes of interest,
* Subnetwork extraction – extracting and exporting a subnetwork (that may be the result of previous steps),
* Inclusion of experimental data – importing experimental data (logFC and p-values values), and adding as context to the networks,
* Cytoscape Automation18 – loading, styling, highlighting, and exporting results of network analysis performed in Python, in the Cytoscape interface using py4cytoscape1, and
* Interoperability with the DiNAR application, allowing integration and visualisation of multiple condition high-throughput data in a knowledge graph context.



## Use cases
### Extracting neighbourhoods for genes of interest

### Finding paths between genes of interest

### Interactive visulaisation using Cytoscape automation

## Additional resources

SKM-tools ( https://github.com/NIB-SI/skm-tools) is a collection of Python scripts and notebooks, incorporating network analysis and visualisation tools, that facilitates interrogation of  CKN and PSS with targeted questions beyond the scope of the web application. Included are demonstrations of:

* Creating a networkX1 graph object – allowing access to the multitude of graph analysis and graph operation functions available in the networkX library,
* Network edge filtering – removing edges or nodes from the network that are not relevant to the biological question at hand,
* Neighbourhood extraction – identifying the immediate interactors of target entities,
* Shortest path analysis – identifying directed (or undirected) paths between source and target nodes of interest,
* Subnetwork extraction – extracting and exporting a subnetwork (that may be the result of previous steps),
* Inclusion of experimental data – importing experimental data (logFC and p-values values), and adding as context to the networks,
* Cytoscape Automation18 – loading, styling, highlighting, and exporting results of network analysis performed in Python, in the Cytoscape interface using py4cytoscape1, and
* Interoperability with the DiNAR19 application, allowing integration and visualisation of multiple condition high-throughput data in a knowledge graph context.

## Requirements

The required non-default library is networkX (https://networkx.org/).

Additional optional libraries are
* py4cytoscape (https://py4cytoscape.readthedocs.io/), necessary if you wish to view the results in Cytoscape.
* pdf () necessary if you wish to export multiple network views from Cytoscape to a single file.

## Publication

### Case Study 1



### Case Study 2










## Examples



## DiNAR




## Notes (and warnings)

### Network layouts

More information on layout automation:
- [https://github.com/cytoscape/cytoscape-automation/blob/master/for-scripters/Python/network-layout.ipynb](Network layout notebook from cytoscape-automation wiki)

#### Copycat
If the Copycat Layout app is not installed in Cytoscape, you can install it from the App store.

#### yFiles
yFiles does not support Cytoscape Automation.
For small networks, often the yFiles Organic layout is the most optimal. You can set it as you default layout in the layout settings,
and then not apply a layout in the automation.



















## Additional Resources:

- [https://github.com/cytoscape/cytoscape-automation/wiki](Cytosape Automation Wiki)
- [https://py4cytoscape.readthedocs.io/en/latest/](py4cytoscape.readthedocs.io)
>>>>>>> Stashed changes
