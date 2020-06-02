import pandas as pd
import os

### Add headers to motif_found_final
with open('motif_found_final.txt','r') as original: data = original.read()

#with open('motif_found_final.txt','w') as modified:
#    modified.write('Module\tMotif1\tMotif2\tno.intersect\tno.promoters_with_paired_motifs\tno.testGenes\tRaw p-value\tAdjusted p-value (BH)\tAdjusted p-value (Bonf)\tAdjusted p-value (Global Bonf)\tGenes\n' + data)
# new headers for webtool
with open('motif_found_final.txt','w') as modified:
    modified.write('Module\tMotif1\tMotif2\tnointersect\tnopromoterswithpairedmotifs\tnotestGenes\tRawpvalue\tAdjustedpvalueBH\tAdjustedpvalueBonf\tAdjustedpvalueGlobalBonf\tGenes\n' + data)

