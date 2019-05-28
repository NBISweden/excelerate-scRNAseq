Differential expression
================

*Created by Olga Dethlefsen & Ståle Nygård. Parts of the code were adapted from <https://nbisweden.github.io/workshop-scRNAseq/>, courtesy of Åsa Björklund*

Content overview
----------------

-   [Introduction](session-de.md) (~15 min)
-   [Common methods](session-de-methods.md) (~30 min)
-   [Methods performance and evaluation](session-de-methods-evaluation.md) (~30 min)

------------------------------------------------------------------------

Learning objectives
-------------------

-   define differential expression DE
-   identify properties of scRNA-seq data influencing DE
-   assess data distributions

------------------------------------------------------------------------

#### What is differential expression testing

-   taking read count data &
-   performing statistical analysis to discover quantitative changes in expression levels between experimental groups
-   i.e. to decide whether, for a given gene, an observed difference in read counts is significant (greater than what would be expected just due to natural random variation)

#### DE is an "old problem"

-   known from bulk RNA-seq and microarray studies
-   in fact building on one of the most common statistical problems, i.e comparing groups for statistical differences

#### Single-cell vs bulk RNA-seq count matrices

<figure>
<!-- <img src="session-de-files/images/single-cell-vs-bulk.pdf" width="400" height="400"> -->
<img src="session-de-files/images/methods-stats.png">
<figcaption>
*O**u**t**c**o**m**e*<sub>*i*</sub> = *M**o**d**e**l*<sub>*i*</sub> + *e**r**r**o**r*<sub>*i*</sub>
</figcaption>
</figure>
Characteristics of scRNA-seq data
---------------------------------

-   high noise levels (technical and biological factors)
-   low library sizes
-   low amount of available mRNAs results in amplification biases and "dropout events"
-   3' bias, partial coverage and uneven depth (technical)
-   stochastic nature of transcription (biological)
-   multimodality in gene expression; presence of multiple possible cell states within a cell population (biological)

##### Example distributions

![menti](session-de-files/images/intro-distributions.png)

<!-- ## Live-coding -->
<!-- _Idea: add a small example of looking at the data, distributions, number of zeros, boxplots comparisons of selected genes_ -->

------------------------------------------------------------------------

[Back to Schedule](../schedule.md)
==================================

[Next to Methods](session-de-methods.html)
==========================================

------------------------------------------------------------------------
