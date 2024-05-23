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