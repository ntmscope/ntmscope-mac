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