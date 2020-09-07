# Project for Network Analysis

This project seeks to replicate the “immunization robustness” described by Rosenblatt et al. (2020) on the Colorado Springs network and compares it with the “traditional robustness” described by Martin and Niemeyer (2019) with respect to the same network and four different network centrality measures.

Summary of the files: 

* NetworkAnalysis_FinalProject.py: Main script which runs the outbreak simulations within randomly generated configuration networks based on the Colorado Springs network
* NetworkAnalysis_FinalProject.R: Script to plot the final results of the simulations 
* data_ ...: centrality-specific data results of the simulation outbreaks (betweenness centrality, degree centrality, eigenvector centrality and pagerank were examined)
* edgelist.truecolsprings.csv: edgelist of original Colorado Springs network
* EstimatingCentralityRobustness.py: Functions created by Martin and Niemeyer (2019) to calculate "traditional robustness"


Sources: 

* Rosenblatt, S., Smith, J., Gauthier, G., & Hébert-Dufresne, L. (2020). Immunization strategies in networks with missing data. _PLOS Computational Biology_.
* Martin, C., & Niemeyer, P. (2019). Influence of measurement errors on networks: Estimating the robustness of centrality measures. _Network Science, 7(2)_, 180–195. https://doi.org/10. 1017/nws.2019.12

