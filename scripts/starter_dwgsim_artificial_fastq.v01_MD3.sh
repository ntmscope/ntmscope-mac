#!/bin/bash
#
###	Tool to generate artificial fastq files from a reference
#
#	An alternative would be randomreads.sh from the BBmap suite
###
#
#Program: dwgsim (short read simulator)
#Version: 0.1.11
#Contact: Nils Homer <dnaa-help@lists.sourceforge.net>
#
#Usage:   dwgsim [options] <in.ref.fa> <out.prefix>
# 1)	Change parameters in script if necessary e.g. in notepad++
# You definitely need to change input (where the initial fasta files are stored) and output (where you want to move fastQ files to), rest is up to you (i.e error rate, coverage, read length, quality,…)


# In command line:
# 2)	Cd to folder where FASTA files are stored
# 3)	bash $PathToScript/ScriptName.sh

#
#Options:
#         -e FLOAT      per base/color/flow error rate of the first read [from 0.020 to 0.020 by 0.000]
#         -E FLOAT      per base/color/flow error rate of the second read [from 0.020 to 0.020 by 0.000]
#         -i            use the inner distance instead of the outer distance for pairs [False]
#         -d INT        outer distance between the two ends for pairs [500]
#         -s INT        standard deviation of the distance for pairs [50.000]
#         -N INT        number of read pairs (-1 to disable) [-1]
#         -C FLOAT      mean coverage across available positions (-1 to disable) [100.00]
#         -1 INT        length of the first read [70]
#         -2 INT        length of the second read [70]
#         -r FLOAT      rate of mutations [0.0010]
#         -F FLOAT      frequency of given mutation to simulate low fequency somatic mutations [0.5000]
#                           NB: freqeuncy F refers to the first strand of mutation, therefore mutations
#                           on the second strand occour with a frequency of 1-F
#         -R FLOAT      fraction of mutations that are indels [0.10]
#         -X FLOAT      probability an indel is extended [0.30]
#         -I INT        the minimum length indel [1]
#         -y FLOAT      probability of a random DNA read [0.05]
#        -n INT        maximum number of Ns allowed in a given read [0]
#        -c INT        generate reads for [0]:
#                           0: Illumina
#                           1: SOLiD
#                           2: Ion Torrent
#         -S INT        generate reads [0]:
#                           0: default (opposite strand for Illumina, same strand for SOLiD/Ion Torrent)
#                           1: same strand (mate pair)
#                           2: opposite strand (paired end)
#         -f STRING     the flow order for Ion Torrent data [(null)]
#         -B            use a per-base error rate for Ion Torrent data [False]
#         -H            haploid mode [False]
#         -z INT        random seed (-1 uses the current time) [-1]
#         -M            generate a mutations file only [False]
#         -m FILE       the mutations txt file to re-create [not using]
#         -b FILE       the bed-like file set of candidate mutations [(null)]
#         -v FILE       the vcf file set of candidate mutations (use pl tag for strand) [(null)]
#         -x FILE       the bed of regions to cover [not using]
#         -P STRING     a read prefix to prepend to each read name [not using]
#         -q STRING     a fixed base quality to apply (single character) [not using]
#         -Q FLOAT      standard deviation of the base quality scores [2.00]
#         -s INT        standard deviation of the distance for pairs [50.000]
#         -h            print this message

#Regarding Error rates
#If you let -e and -E default: (0.02) --> you get a Phred (quality) score of 17 ( log10(0.02) * -10 ) for every position
#-e and -E at 0.001 --> Phred quality of 30

Input="$HOME/auto/Thecus_SeqData/NTMscope_Nwetzstein/Drafts/MAC/Assemblies/ToCheck"
Output="$HOME/auto/Thecus_SeqData/NTMscope_Nwetzstein/Drafts/MAC/Assemblies/ToCheck"
coverage="100"
readlenght="150"
OuterDistance="500"
sdOuterDistance="50"
quality=0.001
date=$(date +%Y-%m-%d)


for file in $Input/*.fasta 
do 
echo $file
id=`echo $file | sed -n 's/.fasta//p'`;
id2=$(basename $file .fasta)
echo $id;
echo $id2

if [ -f $id.bfast.fastq ]
	then
		echo "Sample was already analysed"
	else
		dwgsim -1 $readlenght -2 $readlenght -d $OuterDistance -s $sdOuterDistance -C $coverage -e $quality -E $quality -r 0 -F 0 -R 0 -X 0 -y 0 -H $file $id
		echo "Sample was not yet analysed"
fi
done

rm *bfast.fastq
rm *mutations.txt
rm *mutations.vcf


gzip $Input/*.fastq
mv *.fastq.gz $Output/
cd $Output
mmv "*.bwa.read*.fastq.gz" "#1_artificial_dwgsim_"$date"_"$readlength"bp_R#2.fastq.gz"
#HA-GCF-001679305.1.bwa.read1.fastq.gz

#mmv "*_AP*_*read*.fastq.gz" "#1-AP#2_artificial_dwgsim_"$date"_"$readlenght"bp_R#4.fastq.gz"
#mmv "*_NZ-*_*read*.fastq" "#1-NZ-#2_artificial_dwgsim_"$date"_"$readlenght"bp_R#4.fastq"
#mmv "*_NC-*_*read*.fastq" "#1-NC-#2_artificial_dwgsim_"$date"_"$readlenght"bp_R#4.fastq"
#mmv "*_artificial_*read*.fastq" "#1_artificial_dwgsim_"$date"_"$readlenght"bp_R#3.fastq"

#gzip *.fastq

exit 
 
