#!/usr/bin/env Rscript
#MAC Ariba Analyzer, R script for macrolide and aminoglycoside resistance prediction from MAC_ariba output and its graphical depiction
#All MAC ariba output should be in one folder (ariba_folder), subfolders should bear sequence names

#load libraries---------------------------------------------------------------------
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)

#DATACOMBINATION--------------------------------------------------------------------
#Combine all ariba outputs into one dataframe
#load paths for tsv-files
ariba_folder <- "ariba_output" #args[1] #Put in folder with ariba output, one folder per strain required

all_paths <- list.files(path=ariba_folder,
                        pattern = "^report.tsv",
                        full.names = TRUE, 
                        recursive = TRUE)

# read file content and leave out comment sign
all_content <-
  all_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = "\t",
         encoding = "UTF-8",
         comment.char = "")

# read dir name of file
all_dirnames <- all_paths %>%
  dirname() %>%
  as.list()

#Replace folder path
all_sequence_names <- as.list(str_replace(all_dirnames,"YOUR_FOLDER",""))

# combine file content list and file name list
all_lists <- mapply(c, all_content, all_sequence_names, SIMPLIFY = FALSE)

# unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)

# change column name
names(all_result)[32] <- "isolate_sequence"

#Create output folder
dir.create("ariba_summary")

#Write csv of joined data
write.csv2(file="ariba_summary/ariba_output_joined.csv",all_result,row.names = FALSE)

#Create first line of output_frame with sequence names
output_frame <- data.frame(as.character(all_sequence_names))
colnames(output_frame) <- "isolate_sequence"


#KNOWN AND UNKNOWN VARIANTS---------------------------------------------------------
#Filter known variants--------------------------------------------------------------
known_variants <- all_result[all_result$has_known_var ==1,]

#Exploratory graphs
pdf("ariba_summary/known_variants.pdf")
ggplot(known_variants)+
  geom_bar(aes(x=known_var_change, fill=ref_name))+
  facet_grid(vars(ref_name), scales="free_y",space="free",drop=TRUE)+
  coord_flip()
dev.off() 

#Filter unknown variants------------------------------------------------------------
unknown_variants <- all_result[all_result$known_var ==0,]

#Exploratory graphs
pdf("ariba_summary/unknown_variants.pdf",height=12)
ggplot(unknown_variants)+
  geom_bar(aes(x=ref_ctg_change, fill=ref_name))+
  facet_grid(vars(ref_name), scales="free_y",space="free",drop=TRUE)+
  coord_flip()
dev.off()

#RRL-MUTATION AND INTERPRETATION----------------------------------------------------
#Filter known variants of rrl
rrl_known_variants <- known_variants[known_variants$ref_name=="rrl",]

#Check if strains have more than one variant in rrl
for (i in duplicated(rrl_known_variants$isolate_sequence)){
  if (i == FALSE){
  }else{
    print("Warning the following isolate has more than one known variant in rrl:")
    print(rrl_nown_variants$isolate_sequence[i])  
  }
}  

#Calculate coverage of rrl gene
rrl_known_variants$rrl_cov <- rrl_known_variants$ref_base_assembled/rrl_known_variants$ref_len
rrl_known_variants$rrl_SNP <- rrl_known_variants$ref_ctg_change
rrl_known_variants$rrl_function <- rrl_known_variants$var_description

#rrl-interpretation
rrl_interpretation <- rrl_known_variants[,c(32:35)]
write.csv(file="ariba_summary/rrl_interpretation.csv",rrl_interpretation,row.names = FALSE)

#Join rrl_interpretation to outputframe
output_frame <- left_join(output_frame,rrl_interpretation, by ="isolate_sequence")

#Code macrolide resistance
output_frame$rrl_code <- 0
output_frame[grep("res",output_frame$rrl_function),]$rrl_code <- 1


#RRS MUTATIONS AND INTERPRETATION---------------------------------------------------
#Filter know variants of rrs
rrs_known_variants <- known_variants[known_variants$ref_name=="rrs",]

#Check if strains have more than one variant in rrs
for (i in duplicated(rrs_known_variants$isolate_sequence)){
  if (i == FALSE){
  }else{
    print("Warning the following isolate has more than one known variant in rrs:")
    print(rrs_known_variants$isolate_sequence[i])  
  }
}    

#Calculate coverage of rrs gene
rrs_known_variants$rrs_cov <- rrs_known_variants$ref_base_assembled/rrs_known_variants$ref_len
rrs_known_variants$rrs_SNP <- rrs_known_variants$ref_ctg_change
rrs_known_variants$rrs_function <- rrs_known_variants$var_description

#rrs-interpretation
rrs_interpretation <- rrs_known_variants[,c(32:35)]
write.csv(file="ariba_summary/rrs_interpretation.csv",rrs_interpretation,row.names = FALSE)

#Join rrs_interpretation to outputframe
output_frame <- left_join(output_frame,rrs_interpretation, by ="isolate_sequence")

#Code aminoglycoside resistance
output_frame$rrs_code <- 0
output_frame[grep("res",output_frame$rrs_function),]$rrs_code <- 1


#CREATE HEATMAP WITH CODED MACROLIDE AND AMINOGLYCOSIDE RESISTANCE----------------
#Extract coded sequences
heatmap_frame <- output_frame[,c("isolate_sequence","rrl_code","rrs_code")]
heatmap_frame_long <- melt(heatmap_frame)

#Create heatmap
pdf("ariba_summary/resistance_heatmap.pdf", height =20, width = 8)
ggplot(heatmap_frame_long)+
  geom_tile(aes(y=isolate_sequence,x=variable,fill=as.character(value)))+
  scale_fill_manual(values=c("white","#F8766D"),labels=c("susceptible","resistant"))+
  labs(fill="Resistance prediction")
dev.off()

#Summary output of heatmap_frame
summary(as.factor(heatmap_frame$rrl_code))
summary(as.factor(heatmap_frame$rrs_code))

#Write output frame to csv
write.csv(file="ariba_summary/output_frame.csv",output_frame,row.names = FALSE)