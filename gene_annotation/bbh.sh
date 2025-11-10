#This is a script to identify the Best Blast Hit with the Uniref100 (Embryophyta) database using MMseqs2 (v17.b804f).
#target: Uniref100 database
#query: Protein file output from EVidenceModeler
#tag: Name of created file

#!/bin/bash

target=/scratch/axelle/funanno/protein/uniref100_embryophyta.fasta.gz

query=h2.EVM.pep

mytmp=/scratch/axelle/funanno/temp

tag=h2Uniref100BBH

mmseqs easy-search ${query} ${target} ${tag} ${mytmp} --threads 16 --comp-bias-corr 0 --mask 0 -s 0.7 --format-output "query,target,pident,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov,evalue,theader"
