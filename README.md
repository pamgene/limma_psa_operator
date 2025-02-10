# LIMMA Phosphosite Analysis operator for Tercen

##### Description

Differential `Phosphosite Analysis` using the Bioconductor package LIMMA (linear models for micro arrays):
https://bioconductor.org/packages/release/bioc/html/limma.html


##### Usage



Input projection|.
---|---
`y-axis` | Values to perform analysis on. (VSN) normalized or log-transformed data are expected.
`x-axis` | x-axis to indicate the grouping. all n(n-1)/2 group-group comparisons, aka `contrasts`,  are analyzed (`n` is the number of groups).
`rows` | Peptide (spot) ID
`labels`| Labels to indicate observations in the analysis (e.g. Barcode, Row)
`color`| Optional single color to indicate a pairing or blocking factor. This factor is implemented as a random factor in the LIMMA analysis. 
`columns`| Optional supergrouping. The analysis is performed for each column of the cross tab view separately.


Input parameters|.
---|---
`ReverseContrast`| The direction of the contrast (Test vs Control or Control vs Test). True by default, meaning "Test vs Control" type comparison: the condition starting with lower letter in the alphabet is used as the control condition.


Output relations|.
---|---
`contrast`|Character, The two-group comparison analyzed (e.g. Test vs Control)
`logFC`|Numeric, the log fold change for the contrast analyzed, per peptide
`AveExpr`| Numeric, the average expression per peptide
`t`|Numeric, (moderated)  t-value for the contrast analyzed, per peptide
`pvalue`|Numeric, p value for the contrast analyzed, per peptide. 
`FDR`|adjusted p value with False Discovery Rate method (FDR, or BH). The p values are corrected for multiple peptides per comparison, but not for multiple comparisons between conditions. It is left to the user to handle this.
`logp`| Numeric, -log10(p.value)
`rankp`| Ranks p-values from low to high, per contrast
...|factors indicating supergrouping, if applicable

##### Details

LIMMA is used to apply linear modeling to PamChip experiments over all included test conditions.

Subsequently, all pairwise comparisons (contrasts) are analyzed.

LIMMA, instead of analyzing each peptide in isolation, uses information from all the peptides and conditions in the experiment to improve its analysis. This is like asking, "What can we learn from the whole group to help us understand each individual better?" This is done by applying the empirical Bayes method to calculate differences (the t-statistic), borrowing information from all peptides. 

Applying LIMMA leads to results with increased statistical power (i.e. the ability to detect real differences) over individual testing per peptide (with for example t-tests). Therefore, LIMMA can detect smaller differences that might otherwise be missed. This makes the analysis more robust, especially when only limited data is available.







 
 
