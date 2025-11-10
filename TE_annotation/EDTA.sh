#This is the script that was used to run TE annotation with EDTA

 #!/bin/bash
  perl ../EDTA_repo/EDTA.pl --genome h1.fasta --cds cds.fa --sensitive 1 --threads 16 --anno 1 2>&1 | tee h1_EDTA.log &
  perl ../EDTA_repo/EDTA.pl --genome h2.fasta --cds cds.fa --sensitive 1 --threads 16 --anno 1 2>&1 | tee h2_EDTA.log &
