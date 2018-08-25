#!/usr/bin/bash
# -S /bin/bash
#$ -cwd
#/opt/conda/bin/bedtools intersect -wao -v -a gorGor3_VS_nomLeu3_R_variants.csv.vcf -b exomeD.bed > filtered/gorGor3_VS_nomLeu3_R_variants.filt.vcf

for file in ../data/vcf_variants/*.vcf

do

	OUT=../data/vcf_variants/`basename -s .vcf $file`.filt.vcf
	echo $OUT
	/opt/conda/bin/bedtools intersect -wao -v -a $file  -b  exome_cmp.bed >  $OUT
done
