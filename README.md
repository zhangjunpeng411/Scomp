# Scomp
A method for portraying ncRNA synergistic competition in Pancancer

## Background
Generally, the crosstalk relationships between ncRNAs (acting as ceRNAs) and mRNAs are not one-to-one but many-to-many, indicating synergistic competition of ncRNAs. Portraying ncRNA synergistic competition can greatly help to understand the synergistic competition mechanism of ncRNAs in human complex diseases.

## Description of each file
Pancancer.RData: Matched lncRNA, pseudogene and mRNA expression data, survival data, putative ceRNA network from ENCORI and SPONGEdb, and tumor type information.

Scomp.R: Functions for identifying and analyzing ncRNA synergistic competition networks.

Case_studies.R: Scripts of two case studies for identifying and analyzing ncRNA synergistic competition networks.

## Case studies
Paste all files including scripts and datasets into a single folder (set the folder as the directory of R environment), the scripts of two case studies using Scomp is implemented in Case_studies.R. The users can simply run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("Case_studies.R")
```

## The usage of Scomp
For identifying ncRNA synergistic competition networks, users should prepare datasets including lncRNA, pseudogene and mRNA expression data, survival data, putative ceRNA network from ENCORI and SPONGEdb, and tumor type information.. Paste the datasets and our source file (Scomp.R) into a single folder (set the folder as the directory of R environment), users can use the following scripts to identify ncRNA synergistic competition network. For convenience, the datasets prepared for users are from our datasets (Pancancer.RData). 

```{r echo=FALSE, results='hide', message=FALSE}
## Load required R packages and utility functions
library(corpcor)
library(igraph)
library(miRspongeR)
library(survival)
library(mldr)
library(utiml)
library(e1071)
library(plyr)
source("Scomp.R")

## Load data source
load("Pancancer.RData")

## Identifying ncRNA synergistic competition network, the putative ceRNA network is from ENCORI
Scomp_ENCORI_lncRpseudo2lncRpseudo <- Scomp_sc(rbind(ENCORI_filter_lncRmR, ENCORI_filter_pseudomR), cbind(Pancancer_lncRNA_Exp, Pancancer_pseudogene_Exp), Pancancer_mRNA_Exp, senscorcutoff = 0.1)
```
