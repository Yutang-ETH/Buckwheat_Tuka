#!/bin/bash

# 11.03.2025

# map hifi reads to Fabian h1 h2 and dip
h1=/home/yutachen/public/Yutangchen/Fabian/juice_buckwheat/h1/Fabian_h1_final.fasta
h2=/home/yutachen/public/Yutangchen/Fabian/juice_buckwheat/h2/Fabian_h2_final.fasta
hifi=/home/yutachen/public/Yutangchen/Fabian/Fabian_hifi_all.fastq.gz


# make a diploid assembly
cat ${h1} ${h2} > Fabian_dip_final.fasta

# map all hifi reads to the diploid assembly
minimap2 -ax map-hifi -t 30 -2 -K 1g --secondary=no Fabian_dip_final.fasta ${hifi} | samtools sort -O bam -@ 30 -o all_HiFi_to_Fabian_dip.bam
samtools index all_HiFi_to_Fabian_dip.bam

