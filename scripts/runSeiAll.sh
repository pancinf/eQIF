#!/bin/bash
#This script computes and extracts Sei scores
cd sei-framework/
date
for chr in {1..22}; do
	bcftools view ../../data/outputs/gnomAD/gn_v3_red_inf_id2.vcf.gz -o ../../data/outputs/Sei/inVcf/gn_v3_red_inf_id2_seiAll_chr${chr}.vcf -Ov -r chr${chr}

	##
	##Run Sei
	CUDA_VISIBLE_DEVICES=1 sh 1_variant_effect_prediction.sh ../../data/outputs/Sei/inVcf/gn_v3_red_inf_id2_seiAll_chr${chr}.vcf hg38 ../../data/outputs/Sei/ --cuda
	sh 2_varianteffect_sc_score.sh ../../data/outputs/Sei/chromatin-profiles-hdf5/gn_v3_red_inf_id2_seiAll_chr${chr}.ref_predictions.h5 \
	../../data/outputs/Sei/chromatin-profiles-hdf5/gn_v3_red_inf_id2_seiAll_chr${chr}.alt_predictions.h5 \
	../../data/outputs/Sei/ --no-tsv
done
date

