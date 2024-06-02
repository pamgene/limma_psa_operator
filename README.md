# LIMMA Phosphosite Analysis operator for Tercen

##### Description

Differential `Phosphosite Analysis` using the Bioconductor package LIMMA (linear models for micro arrays):
https://bioconductor.org/packages/release/bioc/html/limma.html


##### Usage



Input projection|.
---|---
`y-axis` | values to perform analysis on. (VSN) normalized or log-transformed data are expected.
`x-axis` | x-axis to indicate the grouping. all n(n-1)/2 group-group comparissons, aka `contrasts`,  will be analyzed (`n` is the number of groups).
`rows` | spot ID
`labels`| labels to indicate observations in the analysis.
`color`| optional single color to indicate a pairing or blocking factor. This factor will be implemented as a random factor in the LIMMA analysis.
`columns`| Optional supergrouping. The Analysis will be performed for each column of the cross tab view separately.


Input parameters|.
---|---


Output relations|.
---|---
`contrast`|Character,The two-group comparison analyzed
`logFC`|Numeric, the log fold change for the contrast analyzed, per peptide
`AveExpr`| Numeric, the average expression per peptide
`t`|Numeric, (moderated)  t-value for the contrast analyzed, per peptide
`p.value`|Numeric, p value for the contrast analyzed, per peptide. The p values are not corrected for multiple comparisons between conditions or peptides, It is left to the user to handle this.
`logp`| Numeric, -log10(p.value)
p.rank| Ranks p-values from low to high, per contrast
...|optional factors indicating supergrouping

##### Details

LIMMA is used to apply linear modeling to PamChip experiments over all included test conditions.

Subsequently, all pairwise comparissons (contrasts) are analyzed.

LIMMA's empirical Bayes method is applied to moderate the t-statistic based on information of all analyzed peptides.

Hence, the statistical power of analyzing per peptide pairwise comparisons is increased by using information of all conditions and peptides included in the analysis.






 
 
