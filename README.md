# NTMscope-MAC
This repository comprizes essential scripts and code that have been used in the bioinformatical analyses of the whole genome sequencing project **NTMscope-MAC**. Single bioinformatic tools are summarized below. Sequencing data can be accessed in the European Nucleotide Archive under project number **PRJEB70863**. 
\
\
An additional interactive visualisation of our dataset can be accessed via https://microreact.org/project/ntmscope-mac.
\
\
Further information on NTMscope-related projects can be found on https://www.ntmscope.org. 

# Used bioinformatic tools and scripts
Dependencies of the used bioinformatical tools are not shown.

## NTMseq pipeline
- Most of the scripts below have been implemented in the NTMseq pipeline (Most recent version on https://github.com/ngs-fzb/NTMtools/tree/main/scripts/NTMseq)

## Quality control
- Fastqc v. 0.11.9 (https://github.com/s-andrews/FastQC)
- fastqc_all.sh (custom script)
- Multiqc v. 1.13 (https://multiqc.info/)

## Removal of adapter sequences
- Fastp v. 0.23.2 (https://github.com/OpenGene/fastp)
- starter_fastp_DefaultSettings.sh

## Simulation of reads for assembled sequence data
- dwgsim v. 0.1.12-13 (https://github.com/nh13/DWGSIM)
- starter_dwgsim_artificial_fastq.v01_MD3.sh

## Species designation
- NTM-Profiler v. 0.2.0 (https://github.com/jodyphelan/NTM-Profiler)
- ntm_profiler_all.sh (custom script)
- ntm_profiler_extractor.R (custom script)
- Type Strain Genome Server, TYGS (https://tygs.dsmz.de/)
- fastANI (https://github.com/ParBLiSS/FastANI)

## Phylogeny
- MTBseq v. 1.1.0 (https://github.com/ngs-fzb/MTBseq_source)
```
MTBseq --step TBfull --ref M._avium_DSM44156 --threads YOUR_THREADS --continue --sample YOUR_SAMPLE_LIST
```
- exploratory_tree_upgma.R (custom script)
- Raxml v. 8.2.12 (https://github.com/stamatak/standard-RAxML)
```
raxmlHPC -f a -m GTRGAMMA -s XX_amended_u95_phylo.fasta -n YOUR_NAME -x 12345 -p 12345 -# 500

```
- shovill v. 1.1.0 (https://github.com/tseemann/shovill)
- shovil_all.sh (custom script)
- Mashtree v.1.2.0
- starter_Mashtree.sh
- transmission_cluster_identifier.R (custom script)
- Microreact (https://microreact.org/)
- Microreact project (https://microreact.org/project/wUatE2dwuoHnCH7G4petvb-ntmscope-mac)
  
## Genome annotation
- Prokka v. 1.14.6. (https://github.com/tseemann/prokka)
- prokka_all.sh (custom script)
  
## Plasmid analysis
- SRST2 v.0.2.0 (https://katholt.github.io/srst2/) using a custom database of 152 NTM plasmids downloaded from PLSDB v.2021-06-23-v2
- starter_customDB_SRST2.sh and customDB_to_SRST2.db
- Roary v.3.13.0  (https://github.com/sanger-pathogens/Roary/blob/master/README.md)

## Resistance prediction
- AMRfinderPlus v.3.11.2
- ariba v. 2.14.4 (https://github.com/sanger-pathogens/ariba)
- Mab_ariba (https://github.com/samlipworth/Mab_ariba)
- MAC_ariba_analyser.R (custom script)

## Network analysis
- patient_combinations.R (custom script)
- transmission_forloop.R (custom script)
