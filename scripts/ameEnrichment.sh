#/!bin/bash
#This script conducts TFBS enrichment analysis
mkdir -p ../data/raw/ref/
mkdir -p ../data/raw/ame
mkdir -p ../data/outputs/ame
wget -O ../data/raw/ame/jaspar.meme https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt

##
##Make fasta of regions
samtools faidx ../data/raw/ref/hg38.fa --region-file ../data/outputs/eval/ame/regions.txt -o ../data/outputs/eval/ame/hg38RegRegions.fa 

##
##Mask repeats with Dust
dust ../data/outputs/eval/ame/hg38RegRegions.fa > ../data/outputs/eval/ame/hg38RegRegionsDust.fa

##
##Apply AME
ame --seed 1234 --evalue-report-threshold 0.05 -oc ../data/outputs/eval/ame/ameOutJasparShuf --control --shuffle-- ../data/outputs/eval/ame/hg38RegRegionsDust.fa ../data/raw/ame/jaspar.meme
