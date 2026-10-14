#/!bin/bash
#This script conducts FIMO analysis
motif=$1

##
##Make fasta of regions
samtools faidx ../data/raw/ref/hg38.fa --region-file ../data/outputs/eval/fimo/regions.txt -o ../data/outputs/eval/fimo/hg38RegRegions.fa 

##
##Mask repeats with Dust
dust ../data/outputs/eval/fimo/hg38RegRegions.fa > ../data/outputs/eval/fimo/hg38RegRegionsDust.fa

##
##Run FIMO
fimo --oc ../data/outputs/eval/fimo/ $motif ../data/outputs/eval/fimo/hg38RegRegionsDust.fa
