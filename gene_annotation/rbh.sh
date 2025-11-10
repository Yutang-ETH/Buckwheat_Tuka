#This is a script to identify the Reciprocal Best Hit with the Uniref100 (Embryophyta) database using MMseqs2 (v17.b804f).
#target: Uniref100 database
#query: Protein file output from EVidenceModeler
#tag: Name of created file

#!/bin/bash

target=/scratch/axelle/funanno/protein/uniref_100_embryophyta.fasta.gz

query=h2.EVM.pep

mytmp=/scratch/axelle/funanno/tmp

tag=h2Uniref100RBH

mmseqs easy-rbh ${query} ${target} ${tag} ${mytmp} --format-output "query,target,pident,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov,evalue,theader"
