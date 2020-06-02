import argparse
import numpy as np
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    # p-value threshold
    parser.add_argument('threshold',type=float)
    args = parser.parse_args()
    return args

def bh(pvalues):
    '''
    Computes the Benjamini-Hochberg FDR correction.

    Input:
        * pvals - vector of p-values to correct
    '''
    pvalues = np.array(pvalues)
    n = pvalues.shape[0]
    new_pvalues = np.empty(n)
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues

args = get_args()
# read in motif_found which contains all raw p-values.
motif_found = pd.read_csv('motif_found.txt', sep='\t', index_col=None, header=None)
motif_found = motif_found.sort_values([0,1,2])
motif_found = motif_found.values

### global Bonferroni over all pairs tested across all cluster
pvalues = np.asarray(motif_found[:,3])
# empty vectors for pvals
pvalues_Gbonf = np.zeros(motif_found.shape[0])
for i in range(0,motif_found.shape[0]):
    pvalues_Gbonf[i] = pvalues[i] * motif_found.shape[0]
    if pvalues_Gbonf[i] > 1:
        pvalues_Gbonf[i] = 1

# performed corrections for each cluster in the motif_found file (designated by values in first column)
# pd.unique retains order where np.unique does not
for clust in pd.unique(motif_found[:,0]):
    print(clust)
    motifs = motif_found[motif_found[:,0]==clust]
    pvalues = np.asarray(motifs[:,3])
    # empty vectors for pvals
    pvalues_bonf = np.zeros(motifs.shape[0])
    pvalues_bh = np.zeros(motifs.shape[0])
    # calculate Bonferroni correction
    for i in range(0,motifs.shape[0]):
        pvalues_bonf[i] = pvalues[i] * motifs.shape[0]
        if pvalues_bonf[i] > 1:
            pvalues_bonf[i] = 1
    #calculate BH correction
    pvalues_bh = bh(pvalues)
    motifs = np.column_stack((motifs[:,[0,1,2,4,5,6,3]],pvalues_bh,pvalues_bonf))
    # the first time through create motif_corrected, otherwise append the array.
    if clust==pd.unique(motif_found[:,0])[0]:
        motif_corrected = motifs
    else:
        motif_corrected = np.append(motif_corrected,motifs,axis=0)

motif_corrected = np.column_stack((motif_corrected,pvalues_Gbonf,motif_found[:,7]))
    # print to file with headers this time

with open('motif_found_final.txt','wb') as f:
    np.savetxt(f, motif_corrected,fmt='%s\t%s\t%s\t%d\t%d\t%d\t%1.3e\t%1.3e\t%1.3e\t%1.3e\t%s')

f.close()

