#!/bin/bash

awk '$3 == "gene" {
  split($9,a,";");
  gene_id="NA";
  for (i in a) {
    if (a[i] ~ /^ID=/) {
      gene_id=substr(a[i],4)
    }
  }
  print $1, $4-1, $5, gene_id, ".", $7
}' OFS='\t' data/annotated_h1.gff > genes_h1.bed

awk '$3 == "gene" {
  split($9,a,";");
  gene_id="NA";
  for (i in a) {
    if (a[i] ~ /^ID=/) {
      gene_id=substr(a[i],4)
    }
  }
  print $1, $4-1, $5, gene_id, ".", $7
}' OFS='\t' data/annotated_h2.gff > genes_h2.bed

bedtools coverage \
  -a genes_h1.bed \
  -b data/h1_lr_minimap2_sorted.bam \
  -s \
  -mean \
  > h1_gene_mean_coverage.tsv

bedtools coverage \
  -a genes_h2.bed \
  -b data/h2_lr_minimap2_sorted.bam \
  -s \
  -mean \
  > h2_gene_mean_coverage.tsv

# calculate the average coverage per gene
awk '{sum += $7; n++} END {if (n > 0) print sum / n}' h1_gene_mean_coverage.tsv
awk '{sum += $7; n++} END {if (n > 0) print sum / n}' h2_gene_mean_coverage.tsv

# calculate average gene legnth
awk '{sum += ($3 - $2); n++} END {if (n > 0) print sum / n}' h1_gene_mean_coverage.tsv
awk '{sum += ($3 - $2); n++} END {if (n > 0) print sum / n}' h2_gene_mean_coverage.tsv

# calculate average exon number
awk '$3 == "exon" {
  parent="NA"
  split($9,a,";")
  for (i in a) {
    if (a[i] ~ /^Parent=/) {
      parent=substr(a[i],8)
    }
  }
  exon_count[parent]++
}
END {
  for (p in exon_count) {
    sum += exon_count[p]
    n++
  }
  if (n > 0) print "Average_exon_number =", sum / n
}' data/annotated_h1.gff

awk '$3 == "exon" {
  parent="NA"
  split($9,a,";")
  for (i in a) {
    if (a[i] ~ /^Parent=/) {
      parent=substr(a[i],8)
    }
  }
  exon_count[parent]++
}
END {
  for (p in exon_count) {
    sum += exon_count[p]
    n++
  }
  if (n > 0) print "Average_exon_number =", sum / n
}' data/annotated_h2.gff
