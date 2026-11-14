
# PL4246-Final-Report

This project removes top-ranked connector hubs from a functional brain network, one at a time, and measures how each removal changes the network's modularity and average shortest path length.

**Files**

BiswalNetworkLesioning.R	
The full analysis script to simulate lesions.

AveragedMatrix.txt	Group-averaged functional connectivity matrix. 

region_names_abbrev_file.txt	177 unique short region labels, one per line (for example LLOCsd2).

region_names_full_file.txt	177 full region names, one per line (for example Left Lateral Occipital Cortex superior division).
hub_lesion_results.csv	Output created by the script of lesion outcomes.

**Requirements**

R (2026)
R packages: igraph, igraphdata

**How to run**

Download all four input files and the script into one folder.
Run the script on BiswalNetworkLesioning.R from top to bottom.
The results are written to hub_lesion_results.csv in the same folder.

**Outputs**

hub_lesion_results.csv has the following columns for inspection: short_label, degree, betweenness, participation, degree_rank, betweenness_rank, participation_rank, hub_score, long_label, delta_Q, delta_ASPL, PC_group, abs_delta_Q.
