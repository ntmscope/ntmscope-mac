#Bash script to conduct NTM-profiler analysis on all sequences in WD
#Attention: reads need to end with R1/R2 and to be situated in working directory
#Use in base environment!!!
#!/usr/bin/env bash 

#Read in input folder with sequence data
folder=$(pwd)

#Generate strainlist from R1-reads
ls *_R1.fastq.gz | sed -re 's/_R1.fastq.gz//' > strains.txt

#Read in strainlist from generated file
strainlist=strains.txt

#Create and go to ariba output folder, move strains.txt to ariba_output_folder
mkdir ntm_profiler_output
cp strains.txt ntm_profiler_output/strains.txt
cd ntm_profiler_output

#Ariba loop over strain list, output to ariba_output
while read strain; do
  echo "Running NTM-Profiler analysis for ${strain}"
  mkdir "${strain}"
  cd "${strain}"
  ntm-profiler profile -1 "${folder}/${strain}_R1.fastq.gz" -2 "${folder}/${strain}_R2.fastq.gz" -p "${strain}"
  cd ..
  echo ""
done < ${strainlist}

#Leave ariba outputfolder
cd ..