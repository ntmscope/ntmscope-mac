#!/usr/bin/env Rscript
#ntm_profiler_extractor
#script for extracting species from ntm_profiler output
#input folder should be folder with results / subfolders should bear sequence names

#Load libraries
library(rjson)

#Let use arguments
args = commandArgs(trailingOnly=TRUE)
folder <- args[1]

#Create list with all files within folder
strain_list <- list.dirs(path=folder,full.names=FALSE,recursive=FALSE) #list.dirs(path="YOUR_PATH")
folder_list <- list.dirs(path=folder,recursive=FALSE)

#Create output frame
output_frame <- data.frame(c(""),c(""))
colnames(output_frame) <- c("strain","species")
output_frame <- output_frame[-1,]

#CHANGE to for loop
#Import data from json file    
for (i in 1:length(folder_list)){
  try({
    actual_folder <- folder_list[i]
    strain <- strain_list[i]
    # grep json file
    result <- fromJSON(file=paste(actual_folder,"/",strain,".results.json",sep=""))
    species <- ""
    species <- as.data.frame(result)$species.prediction.species
    #append to output frame
    output_frame[i,]$strain <- strain
    output_frame[i,]$species <- species
  }
  )
}
#Write output table
write.csv(file="ntm_profiler_output.csv",output_frame)