#!/usr/bin/bash
echo Launched terminal part
for file in ../data/variants/vcf_variants/*.vcf
do
    OUT=../data/variants/filt_variants/`basename -s .vcf $file`.filt.vcf
    echo file=$file
	echo out=$OUT
	bedtools intersect -wao -v -a $file  -b  ../data/exome_beds/exomeD.bed >  $OUT

done
