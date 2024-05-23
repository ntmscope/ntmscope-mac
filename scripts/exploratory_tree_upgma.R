#!/usr/bin/env Rscript
#Simple upgma tree from MTBseq tree distance matrix
#Input: distance matrix from MTBseq
#Load libraries
library(ggtree)
library(ape)
library(phangorn)
library(gridExtra)
library(tidyverse)

#Import distance matrix from MTBseq
distance_matrix <- readLines("distance.matrix")
distance_matrix <- strsplit(distance_matrix,"\t")
distance_matrix <- lapply(distance_matrix,function(x){
  length(x) <- max(unlist(lapply(distance_matrix,length)));
  return(x)
})
distance_matrix <- do.call(rbind,distance_matrix)
distance_matrix <- as.data.frame(distance_matrix)
distance_matrix <- head(distance_matrix,-1)
distance_matrix$last_column <- NA
row.names(distance_matrix) <- distance_matrix[,1]
colnames(distance_matrix) <- c("X",distance_matrix$V1)
distance_matrix$X <- NULL

#Delete everything afte first underscore
rownames(distance_matrix) <- gsub("_.+","",rownames(distance_matrix))
colnames(distance_matrix) <- gsub("_.+","",colnames(distance_matrix))

#Create tree from distance matrix without metadata------------------------------------------------------------
all_upgma <- upgma(distance_matrix)
#Outputtree as newick file
write.tree(all_upgma,file="all_upgma.nwk")