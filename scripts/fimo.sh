#/!bin/bash
#This script conducts FIMO analysis
motif=$1

##
##Make fasta of regions
samtools faidx ../data/raw/ref/hg38.fa --region-file ../data/outputs/common/eval/fimo/regions.txt -o ../data/outputs/common/eval/fimo/hg38RegRegions.fa 

##
##Mask repeats with Dust
dust ../data/outputs/common/eval/fimo/hg38RegRegions.fa > ../data/outputs/common/eval/fimo/hg38RegRegionsDust.fa

##
##Run FIMO
fimo --oc ../data/outputs/common/eval/fimo/ $motif ../data/outputs/common/eval/fimo/hg38RegRegionsDust.fa
