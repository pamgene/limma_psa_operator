# LIMMA Phosphosite Analysis operator for Tercen

##### Description

Differential `Phosphosite Analysis` using the Bioconductor package LIMMA (linear models for micro arrays):
https://bioconductor.org/packages/release/bioc/html/limma.html


##### Usage



Input projection|.
---|---
`y-axis` | values to perform analysis on
`x-axis` | x-axis to indicate the grouping. all n(n-1)/2 group-group comparissons, aka `contrasts`,  will be analyzed (`n` is the number of groups).
`rows` | spot ID
`labels`| labels to indicate observations in the analysis.
`color`| optional single color to indicate a pairing or blocking factor. This factor will be implemented as a random factor in the LIMMA analysis.
`columns`| Optional supergrouping. The Analysis will be performed for each column of the cross tab view separately.
Input parameters|.
---|---


Output relations|.
---|---
.|.

##### Details





 
 
