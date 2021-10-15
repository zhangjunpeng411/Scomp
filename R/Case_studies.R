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

## Identifying hub ncRNAs in ncRNA synergistic competition network, the putative ceRNA network is from ENCORI
Scomp_ENCORI_lncRpseudo2lncRpseudo_graph <- make_graph(t(Scomp_ENCORI_lncRpseudo2lncRpseudo[, 1:2]), directed = FALSE)
Scomp_ENCORI_lncRpseudo2lncRpseudo_outdegree <- degree(Scomp_ENCORI_lncRpseudo2lncRpseudo_graph, mode="out")
Scomp_ENCORI_lncRpseudo2lncRpseudo_hub <- names(sort(Scomp_ENCORI_lncRpseudo2lncRpseudo_outdegree[which(Scomp_ENCORI_lncRpseudo2lncRpseudo_outdegree!=0)], decreasing=TRUE))[1:ceiling(0.1*length(which(Scomp_ENCORI_lncRpseudo2lncRpseudo_outdegree!=0)))]

## Network module identification, the putative ceRNA network is from ENCORI
Scomp_ENCORI_lncRpseudo2lncRpseudo_module <- netModule(Scomp_ENCORI_lncRpseudo2lncRpseudo[, 1:2], modulesize = 10)

## Suivival analysis of network modules, the putative ceRNA network is from ENCORI
Scomp_ENCORI_lncRpseudo2lncRpseudo_module_survival <- moduleSurvival(Scomp_ENCORI_lncRpseudo2lncRpseudo_module, cbind(Pancancer_lncRNA_Exp, Pancancer_pseudogene_Exp), Pancancer_survival)

## Classification analysis of network modules, the putative ceRNA network is from ENCORI
Scomp_ENCORI_lncRpseudo2lncRpseudo_classify_baseline <- module.classify(cbind(Pancancer_lncRNA_Exp, Pancancer_pseudogene_Exp), tumor_type, Scomp_ENCORI_lncRpseudo2lncRpseudo_module, method = "baseline")
Scomp_ENCORI_lncRpseudo2lncRpseudo_classify_br <- module.classify(cbind(Pancancer_lncRNA_Exp, Pancancer_pseudogene_Exp), tumor_type, Scomp_ENCORI_lncRpseudo2lncRpseudo_module, method = "br")
Scomp_ENCORI_lncRpseudo2lncRpseudo_hub_classify_br <- module.classify(cbind(Pancancer_lncRNA_Exp, Pancancer_pseudogene_Exp), tumor_type, list(Scomp_ENCORI_lncRpseudo2lncRpseudo_hub), method = "br")

## Identifying ncRNA synergistic competition network, the putative ceRNA network is from SPONGEdb
Scomp_SPONGEdb_lncRpseudo2lncRpseudo <- Scomp_sc(rbind(SPONGEdb_filter_lncRmR, SPONGEdb_filter_pseudomR), cbind(Pancancer_lncRNA_Exp, Pancancer_pseudogene_Exp), Pancancer_mRNA_Exp, senscorcutoff = 0.1)

## Identifying hub ncRNAs in ncRNA synergistic competition network, the putative ceRNA network is from SPONGEdb
Scomp_SPONGEdb_lncRpseudo2lncRpseudo_graph <- make_graph(t(Scomp_SPONGEdb_lncRpseudo2lncRpseudo[, 1:2]), directed = FALSE)
Scomp_SPONGEdb_lncRpseudo2lncRpseudo_outdegree <- degree(Scomp_SPONGEdb_lncRpseudo2lncRpseudo_graph, mode="out")
Scomp_SPONGEdb_lncRpseudo2lncRpseudo_hub <- names(sort(Scomp_SPONGEdb_lncRpseudo2lncRpseudo_outdegree[which(Scomp_SPONGEdb_lncRpseudo2lncRpseudo_outdegree!=0)], decreasing=TRUE))[1:ceiling(0.1*length(which(Scomp_SPONGEdb_lncRpseudo2lncRpseudo_outdegree!=0)))]

## Network module identification, the putative ceRNA network is from SPONGEdb
Scomp_SPONGEdb_lncRpseudo2lncRpseudo_module <- netModule(Scomp_SPONGEdb_lncRpseudo2lncRpseudo[, 1:2], modulesize = 10)

## Suivival analysis of network modules, the putative ceRNA network is from SPONGEdb
Scomp_SPONGEdb_lncRpseudo2lncRpseudo_module_survival <- moduleSurvival(Scomp_SPONGEdb_lncRpseudo2lncRpseudo_module, cbind(Pancancer_lncRNA_Exp, Pancancer_pseudogene_Exp), Pancancer_survival)

## Classification analysis of network modules, the putative ceRNA network is from SPONGEdb
Scomp_SPONGEdb_lncRpseudo2lncRpseudo_classify_baseline <- module.classify(cbind(Pancancer_lncRNA_Exp, Pancancer_pseudogene_Exp), tumor_type, Scomp_SPONGEdb_lncRpseudo2lncRpseudo_module, method = "baseline")
Scomp_SPONGEdb_lncRpseudo2lncRpseudo_classify_br <- module.classify(cbind(Pancancer_lncRNA_Exp, Pancancer_pseudogene_Exp), tumor_type, Scomp_SPONGEdb_lncRpseudo2lncRpseudo_module, method = "br")
Scomp_SPONGEdb_lncRpseudo2lncRpseudo_hub_classify_br <- module.classify(cbind(Pancancer_lncRNA_Exp, Pancancer_pseudogene_Exp), tumor_type, list(Scomp_SPONGEdb_lncRpseudo2lncRpseudo_hub), method = "br")

save.image("Scomp.RData")
