#!/bin/bash

# 07.29.2024

asm_h1=/home/yutachen/public/Yutangchen/Fabian/juice_buckwheat/h1/Fabian_h1_final.fasta
asm_h2=/home/yutachen/public/Yutangchen/Fabian/juice_buckwheat/h2/Fabian_h2_final.fasta

hifi_read_1=/home/yutachen/public/Yutangchen/Fabian/PacBio_Buckwheat/m64141e_240116_144010.hifi_reads.fastq.gz
hifi_read_2=/home/yutachen/public/Yutangchen/Fabian/Revio_Buckwheat/m84151_240405_194459_s3.hifi_reads.bam.fastq.gz

kat comp -t 30 -m 23 -h -H 50000000000 -o ./h1 "${hifi_read_1} ${hifi_read_2}" ${asm_h1}
kat comp -t 30 -m 23 -h -H 50000000000 -o ./h2 "${hifi_read_1} ${hifi_read_2}" ${asm_h2}

kat plot spectra-cn -x 200 -o ./Fabian_h1.png ./h1-main.mx
kat plot spectra-cn -x 200 -o ./Fabian_h2.png ./h2-main.mx
