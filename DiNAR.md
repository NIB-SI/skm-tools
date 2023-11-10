# DiNAR & SKM
[DiNAR](https://github.com/NIB-SI/DiNAR/) (Differential Network Analysis in R) allows for the examination and visualisation of differential (gene) expression within prior knowledge graphs.

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


