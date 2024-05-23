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