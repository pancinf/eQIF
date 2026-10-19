#!/bin/bash
#This script conducts TFBS enrichment analysis

##
##Make fasta of regions
samtools faidx ../data/raw/ref/hg38.fa --region-file ../data/outputs/eval/ame/regions.txt -o ../data/outputs/eval/ame/hg38RegRegions.fa 

##
##Mask repeats with Dust
dust ../data/outputs/eval/ame/hg38RegRegions.fa > ../data/outputs/eval/ame/hg38RegRegionsDust.fa

##
##Apply AME
ame --seed 1234 --evalue-report-threshold 0.05 -oc ../data/outputs/eval/ame/ameOutJasparShuf --control --shuffle-- ../data/outputs/eval/ame/hg38RegRegionsDust.fa ../data/raw/ame/jaspar.meme
