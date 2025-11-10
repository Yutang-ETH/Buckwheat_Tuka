#This is a script to run evidence-based gene prediction in one step.

#the following dependencies were installed before hand 
#mamba install bioconda::agat bioconda::mmseqs2 bioconda::transdecoder bioconda::hisat2 bioconda::minimap2=2.28 bioconda::stringtie bioconda::td2 bioconda::gffread bioconda::samtools conda-forge::singularity bioconda::snakemake bioconda::miniprot bioconda::evidencemodeler=2.1.0

# Version of the used tools
# EvidenceModeler 2.1.0 
# Hisat2 2.2.1 
# Minimap2 2.28
# Mmseqs2 17.b804f
# TD2 5.7.1
# Stringtie 3.0.0
# Samtools 1.22
# Miniprot 0.16

# inputs needed
# genome (fasta file of soft-masked genome)
# lr (ISO-seq long reads)
# sr1 (RNA-seq short reads)
# sr2 (RNA-seq short reads)
# protein (protein reference database)


###################
##start of script##
###################


#!/bin/bash

genome=/scratch/axelle/gene_prediction/h2_soft_masked.fasta
lr=/scratch/axelle/gene_prediction/ISO_seq/iso_seq.fq
sr1=/scratch/axelle/gene_prediction/RNA_seq/rna_seq_1.fq.gz
sr2=/scratch/axelle/gene_prediction/RNA_seq/rna_seq_2.fq.gz
protein=/scratch/axelle/gene_prediction/Protein/uniref_90_embryophyta_renamed.fasta
transcripts=/scratch/axelle/gene_prediction/h2_transcripts.fasta

mkdir long_read_aln
mkdir short_read_aln

#minimap2 

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting minimap2"
minimap2 -a -x splice:hq --cs=long -t 16 -L --eqx -2 -G 20000 --secondary=no $genome $lr | samtools sort -@ 16 -O bam -o long_read_aln/h2_lr_minimap2_sorted.bam &

#hisat2
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting hisat2"
hisat2 -x /scratch/axelle/gene_prediction/short_read_aln/h2_soft_masked -1 $sr1 -2 $sr2 --max-intronlen 20000 --dta -p 16 --rna-strandness RF --fr --no-mixed --no-discordant --omit-sec-seq -S short_read_aln/h2_sr_hisat2.sam --summary-file short_read_aln/h2_rna_short_hisat2.report
samtools sort -@ 16 -O bam -o short_read_aln/h2_sr_hisat2_sorted.bam short_read_aln/h2_sr_hisat2.sam

wait

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Minimap2 done"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Hisat2 done"

#indexing
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting indexing"
samtools index -@ 16 long_read_aln/h2_lr_minimap2_sorted.bam &
samtools index -@ 16 short_read_aln/h2_sr_hisat2_sorted.bam

wait
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Indexing done"

#stringtie
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting stringtie mix"
stringtie --mix -f 0.1 -p 8 -c 3 -j 3 -t -m 100 --rf -o h2_stringtie_mix_assembly.gtf short_read_aln/h2_sr_hisat2_sorted.bam long_read_aln/h2_lr_minimap2_sorted.bam &
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting stringtie long read"
stringtie -L -f 0.1 -p 8 -c 3 -j 3 -t -m 100 -s 5 -g 0 -o h2_stringtie_lr_assembly.gtf long_read_aln/h2_lr_minimap2_sorted.bam

wait
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Stringtie done"

gtf_genome_to_cdna_fasta.pl /scratch/axelle/gene_prediction/h2_stringtie_mix_assembly.gtf $genome > h2_transcripts.fasta
#TD2
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting TD2"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting TD2.LongOrfs"
TD2.LongOrfs -t $transcripts -@ 16

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting mmseqs"
mmseqs easy-search /scratch/axelle/gene_prediction/h2_transcripts/longest_orfs.pep $protein /scratch/axelle/gene_prediction/h2_mmseqresult.m8 temp -s 7.0 --threads 16

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting TD2.predict"
TD2.Predict -t $transcripts --retain-mmseqs-hits /scratch/axelle/gene_prediction/h2_mmseqresult.m8

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting gtf to alignment"
gtf_to_alignment_gff3.pl /scratch/axelle/gene_prediction/h2_stringtie_mix_assembly.gtf > h2_transcripts_alignment_format_TD2.gff3

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting cdna alignment"
cdna_alignment_orf_to_genome_orf.pl h2_transcripts.fasta.TD2.gff3 h2_transcripts_alignment_format_TD2.gff3 /scratch/axelle/gene_prediction/h2_transcripts.fasta > h2_transcripts.fasta.TD2.genome.gff3

echo "[`date '+%Y-%m-%d %H:%M:%S'`] TD2 done"

sed 's/*//' /scratch/axelle/gene_prediction/h2_transcripts.fasta.TD2.pep > h2_transcripts.fasta.TD2.nostar.pep
#galba
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting galba"
umask 002 && singularity exec -B ${PWD}:${PWD} /scratch/axelle/gene_prediction/galba.sif galba.pl --genome=$genome --prot_seq=/scratch/axelle/gene_prediction/h2_transcripts.fasta.TD2.nostar.pep --threads=16 --AUGUSTUS_CONFIG_PATH=/scratch/axelle/gene_prediction/Augustus/config --workingdir=/scratch/axelle/gene_prediction/galba_h2 --species=buckwheat_galba_h2 &


#braker3
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting braker3"
umask 002 && singularity exec -B ${PWD}:${PWD} /scratch/axelle/gene_prediction/braker3.sif braker.pl --genome=$genome --prot_seq=$protein --bam=/scratch/axelle/gene_prediction/short_read_aln/h2_sr_hisat2_sorted.bam --workingdir=/scratch/axelle/gene_prediction/braker3_h2 --AUGUSTUS_CONFIG_PATH=/scratch/axelle/gene_prediction/Augustus/config --species=buckwheat_braker_h2 --AUGUSTUS_ab_initio --threads 16 

wait
echo "[`date '+%Y-%m-%d %H:%M:%S'`] galba done"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] braker3 done"

#miniprot
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting miniprot"
miniprot -t 16 -G 20000 --gff-only --outc 0.8 $genome $protein > h2_miniprot_08.gff3
echo "[`date '+%Y-%m-%d %H:%M:%S'`] miniprot done"


#filter miniprot by identity
grep -v -w -f <(grep 'mRNA' h2_miniprot_08.gff3 | cut -f9 | awk -F '[;=]' -v OFS=" " '{{if($6 < 0.8) print $2}}' | sort | uniq) h2_miniprot_08.gff3 > /scratch/axelle/gene_prediction/h2_miniprot_08_filtered.gff3



###prepare evidence for EVM

#miniprot
awk -F '\t' -v OFS='\t' '
$3 == "CDS" {
    gsub("Parent", "ID", $9);
    split($9, a, ";");
    for (i in a) {
        split(a[i], b, "=");
        if (b[1] == "ID" || b[1] == "Parent") {
            score = b[2] * 100;
            break;
        }
    }
    print $1, "miniprot_protAln", "nucleotide_to_protein_match", $4, $5, score, $7, $8, $9
}' h2_miniprot_08_filtered.gff3 | sort -k1,1 -k4,4n --parallel=16 -S 10G > h2_protein_aln_evidence.gff3

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Successfully created h2_protein_aln_evidence.gff3"

#braker
/scratch/axelle/miniforge3/envs/gene_prediction/bin/EvmUtils/misc/braker_GTF_to_EVM_GFF3.pl braker3_h2/braker.gtf | awk -F '\t' -v OFS='\t' '{ if (NF) print $1, "braker3", $3, $4, $5, $6, $7, $8, $9 }' > h2_braker.gff3
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Successfully created h2_braker.gff3"

#augustus
/scratch/axelle/miniforge3/envs/gene_prediction/bin/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl braker3_h2/Augustus/augustus.ab_initio.gtf | awk -F '\t' -v OFS='\t' '{ if (NF) print $1, "augustus", $3, $4, $5, $6, $7, $8, $9 }' > h2_braker_augustus.gff3
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Successfully created h2_braker_augustus.gff3"

#genemark
/scratch/axelle/miniforge3/envs/gene_prediction/bin/EvmUtils/misc/GeneMarkHMM_GTF_to_EVM_GFF3.pl braker3_h2/GeneMark-ETP/genemark_supported.gtf | awk -F '\t' -v OFS='\t' '{ if (NF) print $1, "gmst", $3, $4, $5, $6, $7, $8, $9 }' > h2_braker_gmst.gff3
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Successfully created h2_braker_gmst.gff3"

#galba
/scratch/axelle/miniforge3/envs/gene_prediction/bin/EvmUtils/misc/braker_GTF_to_EVM_GFF3.pl galba_h2/galba.gtf | awk -F '\t' -v OFS='\t' '{ if (NF) print $1, "galba", $3, $4, $5, $6, $7, $8, $9 }' > h2_galba.gff3
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Successfully created h2_galba.gff3"

#transcript alignment evidence
/scratch/axelle/miniforge3/envs/gene_prediction/bin/EvmUtils/misc/align_GTF_to_align_GFF3.pl h2_stringtie_mix_assembly.gtf stringtie_mix > h2_stringtie_mix_aln.gff3

/scratch/axelle/miniforge3/envs/gene_prediction/bin/EvmUtils/misc/align_GTF_to_align_GFF3.pl h2_stringtie_lr_assembly.gtf stringtie_lr > h2_stringtie_lr_aln.gff3

cat h2_stringtie_mix_aln.gff3 h2_stringtie_lr_aln.gff3 > h2_transcript_aln_evidence.gff3
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Successfully created h2_transcript_aln_evidence.gff3"


#combine gene prediction alignment evidence
cat h2_braker.gff3 h2_braker_augustus.gff3 h2_braker_gmst.gff3 h2_galba.gff3 h2_transcripts.fasta.TD2.genome.gff3 > h2_gene_prediction_evidence.gff3
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Successfully created h2_gene_prediction_evidence.gff3"

#run EVM

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting EVM"
umask 002 && /scratch/axelle/miniforge3/envs/gene_prediction/bin/EVidenceModeler --sample_id h2 --genome /scratch/axelle/gene_prediction/h2_soft_masked.fasta --gene_predictions h2_gene_prediction_evidence.gff3 --protein_alignments h2_protein_aln_evidence.gff3 --transcript_alignments h2_transcript_aln_evidence.gff3 --weights evm_weight.txt --segmentSize 900000 --overlapSize 50000 --CPU 16 -S > h2_evm.log 2> h2_evm.err
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Successfully finished EVM"


