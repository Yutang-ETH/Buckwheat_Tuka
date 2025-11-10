```
# By Yutang
# scaffold the assembly from Hifiasm (Trio) to chromosome-level

# As the heterozygosity of the buckwheat genome is not very high, each haplome needs to be scaffolded separately instead of being scaffolded as a diploid assembly
# For scaffolding, the idea would be first anchor the haplome to a chromosome-level reference to construct pseudo-chromosomes
# Then map Hi-C data to the resulting pseudo-chromosomes to do manual curation

# I took this genome (PL4) from this publication, https://doi.org/10.1038/s41477-023-01474-1, as the reference since in this publication they used a genetic linkage map to construct the pseudo-chromosomes
# The PL4 reference genome was downloaded from here, https://drive.google.com/drive/u/1/folders/1faT2M7-KCBA-hbHebpb29BCHAYZUKMmu
# The path to the assembly in our P drive, P:\Yutangchen\Fabian\FesPL4_r1.1_pseudomolecule_contigs

###############################################################################

# Anchor to PL4
# running the command in the ragtag conda env on The Dude
# Extract only chromosomes from the PL4 assembly
# I also renamed the header in the fasta
sed 's/FesPL4_r1\.1_//' FesPL4_r1.1_pseudomolecule_contigs.fa | sed 's/Chr/chr/' | seqkit grep -p 'chr' -r > FesPL4_chr.fa
# FesPL4_chr.fa contains 8 pseudo-chromosomes with headrs: >chr1, >chr2, ..., >chr8

# Make a new folder for scaffolding, P:\Yutangchen\Fabian\scaffold_buckwheat
# In this folder, create two folders named ref and contig, with ref to store the FesPL4_chr.fa and contig to store the buckwheat haplotype-resolved assembly
# In the contig folder, create folder h1 and h2
# Copy buckwheat.asm.dip.hap1.p_ctg.gfa.fasta to h1 folder and rename the file as Fabian_h1.fasta
# Copy buckwheat.asm.dip.hap2.p_ctg.gfa.fasta to h2 folder and rename the file as Fabian_h2.fasta
# I named the buckwheat genome as Fabian

# Edit the scaffold_buckwheat snakemake pipeline to run ragtag
# Contigs shorter than 50 kb will not be anchored

###############################################################################

# build Hi-C contact map with juicer for manual curation
# running the analysis in scratch folder of The Dude
# Path to the folder, /scratch/yutang/Fabian/juice_buckwheat
# in juice_buckwheat, create folder h1 and h2 to generate a Hi-C contact map for each haplome
# take h1 as an example, in folder h1, create required folders, mkdir fastq references restriction_sites
# copy hic data to fastq folder
# cp ~/public/Yutangchen/Fabian/HiC_Buckwheat/342771_1-HiC117_S1_R1_001.fastq.gz fastq/Hic_reads_R1.fastq.gz
# cp ~/public/Yutangchen/Fabian/HiC_Buckwheat/342771_1-HiC117_S1_R2_001.fastq.gz fastq/Hic_reads_R2.fastq.gz  
# copy the scaffolds to references folder
# cp ~/public/Yutangchen/Fabian/scaffold_buckwheat/ragtag_h1/ragtag.scaffold.renamed.fasta references/Fabian_h1_ragtag.fasta
# Hi-C contact map was built with the snakemake pipeline, juice_buckwheat, for each haplome
# After manually checking the Hi-C contact map, no manual curation was made. Therefore the ragtag scaffolds taken as the final assemblies, with the headers changed
# The BUSCO analyses was done on the ragtag scaffolds, but since these scaffolds are considered the final assemblies, the BUSCO results should also indicate the completeness of the final assemblies
```
