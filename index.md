# Used bioinformatic tools and scripts
## Quality control
- Fastqc v. 0.11.9 (https://github.com/s-andrews/FastQC)
- fastqc_all.sh (custom script)
```
#!/usr/bin/env bash 

#Define working directory
folder=$(pwd)

#Run fastqc loop
mkdir QC
for file in ${folder}/*.fastq.gz;do 
	echo "Analyzing $file ..."
	fastqc $file -o ${folder}/QC
	echo "Analysis completed."
done
```
- Multiqc v. 1.13 (https://multiqc.info/)

## Removal of adapter sequences
- Fastp v. 0.23.2 (https://github.com/OpenGene/fastp)

## Simulation of reads for assembled sequence data
- dwgsim v. 0.1.12-13 (https://github.com/nh13/DWGSIM)

## Species designation
- NTMprofiler (https://github.com/jodyphelan/NTM-Profiler)
- Type Strain Genome Server, TYGS (https://tygs.dsmz.de/)
- fastANI (https://github.com/ParBLiSS/FastANI)

## Phylogeny
- MTBseq v. XXXX (https://github.com/ngs-fzb/MTBseq_source)
```
MTBseq XXXXX
```
- exploratory_tree.R UPGMA (custom script)
```
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

```
- Raxml v. 8.2.12 (https://github.com/stamatak/standard-RAxML)
```
Raxml XXX
```
- shovill v. 1.1.0 (https://github.com/tseemann/shovill)
- shovil_all.sh (custom script)
```
#!/usr/bin/env bash 

#Set output folder
folder=$1

#Get lists for sequences in folder
ls *1.fastq.gz  > R1_list.txt

#Generate assemblies from all fastq.gz files in folger
while read -r sequence; do
	echo ${sequence} > sequence.txt
	R2_sequence=$(sed -re "s/1.fastq.gz/2.fastq.gz/" sequence.txt)
	output_folder=$(sed -re "s/_.+//" sequence.txt)
	echo "Generating assembly from ${sequence} and ${R2_sequence} into ${folder}/${output_folder}"
	shovill --outdir ${folder}/${output_folder} --R1 ${sequence} --R2 ${R2_sequence} --force	
	rm sequence.txt
done < R1_list.txt

#Remove sequence lists
rm R1_list.txt
```
- Mashdistance (???)
- transmission_cluster_identifier.R
```
#Load libraries:
library(tidyverse)
library(reshape2)

#d12 clusters-----------------------------------------------------------------------
#Read in d12 groups
d12_groups <- read.delim("d12_groups.csv")

#Delete ungrouped isolates
d12_groups <- filter(d12_groups,group!="ungrouped")

#Join metadata to groups
d12_groups <- left_join(d12_groups,metadata_all,by=c("sequence"="nrz_id_1"))

#Identify groups containing patient IDs
d12_groups$group <- as.factor(d12_groups$group)
length(levels(d12_groups$group)) 

#For loop to identify clusters with more than one patient
transmission_frame <- d12_groups[1,]  
transmission_frame <- transmission_frame[-1,]  
for (i in levels(d12_groups$group)){
  print(i)
  group_frame <- filter(d12_groups,group==i)
  if (length(unique(group_frame$patient_id))>1){
    transmission_frame <- rbind(transmission_frame,group_frame)
  }
}

#How many clusters?
length(unique(transmission_frame$group))
length(unique(transmission_frame$patient_id))
length(unique(transmission_frame$site))
length(unique(transmission_frame$country))
write.csv(transmission_frame[,c("sequence","group")],file="transmission_frame.csv",row.names=FALSE)

#For loop to identify clusters with more than one site
transmission_frame_site <- d12_groups[1,]  
transmission_frame_site <- transmission_frame_site[-1,]  
for (i in levels(d12_groups$group)){
  print(i)
  group_frame <- filter(d12_groups,group==i)
  if (length(unique(group_frame$site))>1){
    transmission_frame_site <- rbind(transmission_frame_site,group_frame)
  }
}

#How many clusters?
length(unique(transmission_frame_site$group))
length(unique(transmission_frame_site$patient_id))
length(unique(transmission_frame_site$site))

write.csv(transmission_frame_site, file="transmission_frame_diff_site.csv",row.names=FALSE)

#d12 clusters by country------------------------------------------------------------
#For loop to identify clusters with more than one country
transmission_frame_country <- d12_groups[1,]  
transmission_frame_country <- transmission_frame_country[-1,]  
for (i in levels(d12_groups$group)){
  print(i)
  group_frame <- filter(d12_groups,group==i)
  if (length(unique(group_frame$country))>1){
    transmission_frame_country <- rbind(transmission_frame_country,group_frame)
  }
}

#How many clusters?
length(unique(transmission_frame_country$group))
length(unique(transmission_frame_country$patient_id))
length(unique(transmission_frame_country$country))
length(unique(transmission_frame_country$site))

#d12 clusters by continent----------------------------------------------------------
#For loop to identify clusters with more than one continent
transmission_frame_continent <- d12_groups[1,]  
transmission_frame_continent <- transmission_frame_continent[-1,]  
for (i in levels(d12_groups$group)){
  print(i)
  group_frame <- filter(d12_groups,group==i)
  if (length(unique(group_frame$continent))>1){
    transmission_frame_continent <- rbind(transmission_frame_continent,group_frame)
  }
}

#How many clusters?
length(unique(transmission_frame_continent$group))
length(unique(transmission_frame_continent$patient_id))
length(unique(transmission_frame_continent$continent))
length(unique(transmission_frame_continent$site))

write.csv(transmission_frame_continent, file="transmission_frame_continent.csv",row.names=FALSE)

```
- Microreact (https://microreact.org/)

## Plasmid analysis
- SRST v.0.2.0 (https://katholt.github.io/srst2/)
- Roary v. XXXXX (https://github.com/sanger-pathogens/Roary/blob/master/README.md)
- blastP

## Genome annotation, resistance genes
- ariba v. XXXX (https://github.com/sanger-pathogens/ariba)
- Mab_ariba (https://github.com/samlipworth/Mab_ariba)
- MAC_ariba_analyser.R (custom script)
```
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

```
- AMRfinderPlus v.3.11.2
- Prokka v. 1.14.6. (https://github.com/tseemann/prokka)
- prokka_all.sh
```
#!/usr/bin/env bash 

#Set output folder
folder=$1

#Activate conda environemnt
#conda activate prokka

#Get lists for sequences in folder
ls *.fasta  > R1_list.txt

#Generate assemblies from all fastq.gz files in folger
while read -r sequence; do
	echo ${sequence} > sequence.txt
	output_folder=$(sed -re "s/\\.fasta//" sequence.txt)
	echo "Generating annotation from ${sequence} into ${folder}/${output_folder}"
	prokka --outdir ${folder}/${output_folder} ${sequence} --force 	
	rm sequence.txt
done < R1_list.txt

#Remove sequence lists
rm R1_list.txt
```

## Network analysis
- patient_combinations.R (custom script)
```
#Script for possible combinations of patients

#Load libraries
library(reshape2)
library(lubridate)
library(DescTools)


#Possible combinations without repetition
patients_combs <- t(combn(single_patients,2))
patients_combs <- data.frame(pat1=patients_combs[,1],pat2=patients_combs[,2])
#patients_combs$transmission <- 1 #Cave put coding here for possibility of transmission
patients_combs
patients_combs_inverse <- patients_combs
patients_combs_inverse$pat1 <- patients_combs$pat2
patients_combs_inverse$pat2 <- patients_combs$pat1
patients_combs_double <- rbind(patients_combs,patients_combs_inverse)

```
- transmission_forloop.R (custom script)
```
#Load libraries
library(lubridate)
library(DescTools)

#For loop for possible transmission occasion between all possible patients-------
setwd("/home/nils/Desktop/NTM/10_ntmscope_MAC/ntmscope_MAC_network/output_tma")

#Create file for output
patients_combs$tm_no <- NULL #reset results
sink("ta_results.txt", split=TRUE)
#Start for loop
for (x in 1:nrow(patients_combs)){
  print(paste(patients_combs[x,1],"vs.",patients_combs[x,2]), quote=FALSE)
  print("---------------------------------", quote=FALSE)
  print("", quote=FALSE)
  #Object for the two patients of the comparison
  pat1_time <- sum_frame[sum_frame$newkey1==paste(patients_combs[x,1]),]
  pat2_time <- sum_frame[sum_frame$newkey1==paste(patients_combs[x,2]),]
  count <- 0 # Set transmission opportunity count to zero
  #Nested for loop for the two regarding patients
  for (i in 1:length(pat1_time$date1)){
    for (y in 1:length(pat2_time$date1)){
      result <- as.Date(c(pat1_time$date1[i],pat1_time$date2[i])) %overlaps% as.Date(c(pat2_time$date1[y],pat2_time$date2[y]))
      if(result){
        count <- count+1
        print(count)
        print(paste(pat1_time$date1[i]), quote=FALSE)
        print(paste(pat1_time$date2[i]), quote=FALSE)
        print(paste(pat1_time$ward[i]), quote=FALSE)
        print(pat2_time$date1[y], quote=FALSE)
        print(pat2_time$date2[y], quote=FALSE)
        print(paste(pat2_time$ward[y]), quote=FALSE)
        print(result, quote=FALSE)
        print("", quote=FALSE)
      }
    }
    patients_combs$tm_no[x] <- count
  }
}
sink()
View(patients_combs)
write.csv(file="patients_combs_result.csv",patients_combs)
```
- extraction_script.py