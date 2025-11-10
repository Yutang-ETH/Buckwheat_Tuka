#!/bin/bash

# 25.04.2023
# assemble the buckwheat gneome using Hifiasm with PacBio HiFi data, Hi-C data, parental data
# done by Yutang


# run hifiasm Hi-C mode
myhic=/home/yutachen/public/Yutangchen/Fabian/HiC_Buckwheat
myhifi=/home/yutachen/public/Yutangchen/Fabian/PacBio_Buckwheat
myrevio=/home/yutachen/public/Yutangchen/Fabian/Revio_Buckwheat

# try hifiasm (Hi-C)
hifiasm -o ./buckwheat.asm -t 48 --h1 $myhic/342771_1-HiC117_S1_R1_001.fastq.gz --h2 $myhic/342771_1-HiC117_S1_R2_001.fastq.gz $myhifi/m64141e_240116_144010.hifi_reads.fastq.gz $myrevio/m84151_240405_194459_s3.hifi_reads.bam.fastq.gz 

# extract fasta 
awk '/^S/{print ">"$2;print $3}' buckwheat.asm.hic.hap1.p_ctg.gfa > buckwheat.asm.hic.hap1.p_ctg.gfa.fasta
awk '/^S/{print ">"$2;print $3}' buckwheat.asm.hic.hap2.p_ctg.gfa > buckwheat.asm.hic.hap2.p_ctg.gfa.fasta
awk '/^S/{print ">"$2;print $3}' buckwheat.asm.hic.p_ctg.gfa > buckwheat.asm.hic.p_ctg.gfa.fasta
awk '/^S/{print ">"$2;print $3}' buckwheat.asm.hic.p_utg.gfa > buckwheat.asm.hic.p_utg.gfa.fasta

# get the assembly statistics 
seqkit stats -N 50,90 -a -b buckwheat.asm.hic.p_ctg.gfa.fasta buckwheat.asm.hic.p_utg.gfa.fasta buckwheat.asm.hic.hap1.p_ctg.gfa.fasta buckwheat.asm.hic.hap2.p_ctg.gfa.fasta > assembly_hifiasm_hic.stats

# try trio binning
hifiasm -o ./buckwheat.asm -t 48 -1 pat.yak -2 mat.yak $myhifi/m64141e_240116_144010.hifi_reads.fastq.gz $myrevio/m84151_240405_194459_s3.hifi_reads.bam.fastq.gz

# extract fasta
awk '/^S/{print ">"$2;print $3}' buckwheat.asm.dip.hap1.p_ctg.gfa > buckwheat.asm.dip.hap1.p_ctg.gfa.fasta
awk '/^S/{print ">"$2;print $3}' buckwheat.asm.dip.hap2.p_ctg.gfa > buckwheat.asm.dip.hap2.p_ctg.gfa.fasta

# get the assembly statistics
seqkit stats -N 50,90 -a -b buckwheat.asm.dip.hap1.p_ctg.gfa.fasta buckwheat.asm.dip.hap2.p_ctg.gfa.fasta > assembly_hifiasm_trio.stats
