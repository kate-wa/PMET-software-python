#!/bin/bash

#inputs:

sh /home/username/PMET/PMET-master/PMET_index/PMETindexgenome.sh \
-r /home/username/PMET/PMET-master/PMET_index/scripts \
-o /home/username/PMET/ds1_output_k5n5k_uYes_index \
-i gene_id= -k 5 -n 5000 -p 1000 -j 10 -u Yes -v NoOverlap \
/home/username/PMET/input/Arabidopsis_thaliana.TAIR10.dna_rm.toplevel.fa \
/home/username/PMET/input/Arabidopsis_thaliana.TAIR10.39.gff3 \
/home/username/PMET/input/ArabidopsisPBM_20140210.meme

