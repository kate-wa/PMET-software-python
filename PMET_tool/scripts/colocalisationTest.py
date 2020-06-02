## colocalisation Test
## This version tests pairs of motifs for colocalisation using the contents of
## the fimohits file (the for loop that goes through this folder is within this
## script now)
## Charlotte Rich 22.1.18

## Last edit: 5.2.18 for loop through motif pairs is now within this script.
## LAst upload to server 5.2.18 13.29

import argparse
import numpy as np
import pandas as pd
# import scipy.stats as sps
import math
import time
import os
import pickle

# need to add a binomial test to presence absence test.
def get_args():
    # get the arguments (file1,fimo2,tar,n)
    parser = argparse.ArgumentParser()
    #fimo1 path/filename
    parser.add_argument('ICthreshold',type=int)
    parser.add_argument('universe',type=int)
    parser.add_argument('pathToIndex',type=str)
    parser.add_argument('testgenesfile',type=str)
    args = parser.parse_args()
    return args


def overlapcheck(span1, span2):
    #helper function to assess whether span1 and span2 overlap
    if span1[0] >= span2[0] and span1[0] <= span2[1]:
        return True
    if span1[1] >= span2[0] and span1[1] <= span2[1]:
        return True
    return False


def intersect(a, b):
    return list(set(a) & set(b))

def binomial_cdf(x,n,p):
    """
    calculates the binomial-pvalue much more efficiently than sps.binom.sf
    """
    cdf = 0
    b = 0
    for k in range(x+1):
        if k > 0:
            b += + np.log(n-k+1) - np.log(k)
        log_pmf_k = b + k * np.log(p) + (n-k) * np.log(1-p)
        cdf += np.exp(log_pmf_k)
    return cdf

def geom_binom_test(coin,promsize,mot_len):
    binom_p=[]
    flips = 2*(promsize - mot_len + 1)
    for k in range(0,len(coin)):
        vals=np.array(coin[:k+1],dtype=float)
        geom = geo_mean(vals)
        #binom_p.append(sps.binom.sf(k,flips,geom))
        binom_p.append(1-binomial_cdf(k,flips,geom))
    binom_p = np.asarray(binom_p)
    #print(binom_p)
    return(binom_p)

def geo_mean(iterable):
    a = np.log(iterable)
    return np.exp(a.sum()/len(a))


def log_hypergeometric(x, G1, G2, U):
    """
    The computation of the log-scale p-value of the pairwise hypergeometric test.
    Mirrored after Matlab function proposed by Meng et al. (2009).

    Input:
        * x - the size of the overlap of the two sets you're comparing
        * G1 - the size of the first set
        * G2 - the size of the second set
        * U - the size of the test universe (encompassing both sets and all other elements)
    """
    #if any of the input arguments are zeros, then return log p value of 0
    if (x*G1*G2*U == 0):
        return 0
    #just gonna compile their logftable thing within the function and use that as needed
    #it would appear that the largest index of it they ever access is U+1
    logf = np.log(np.arange(1,U+2))
    logf = np.insert(logf,0,0)
    logf = np.cumsum(logf)
    #there we go. one logftable ready. can commence computation of actual teststuff
    minG = min(G1,G2)
    a = []
    for ind in range(1,(minG-x+2)):
        hold = x+ind-1
        a.append(-logf[hold]-logf[G1-hold]-logf[G2-hold]-logf[U+hold-G1-G2])
    a = np.asarray(a)
    amax = max(a)
    return logf[G1]+logf[U-G1]+logf[G2]+logf[U-G2]-logf[U]+amax+np.log(sum(np.exp(a-amax)))

### we could multiprocess this, across the genes.
def subset_fimo(fimo1_1,fimo2_1,i):
    mot1=fimo1_1[fimo1_1[:,1]==i[1]]
    #mot1=fimo1_1[np.nonzero(np.in1d(fimo1_1[:,1],i[1]))[0],:]
    mot1 = mot1[np.argsort(mot1[:,6])]
    mot2=fimo2_1[fimo2_1[:,1]==i[1]]
    #mot2=fimo2_1[np.nonzero(np.in1d(fimo2_1[:,1],i[1]))[0],:]
    mot2 = mot2[np.argsort(mot2[:,6])]
    return(mot1,mot2)

def defineOverlapsBetweenTwoMotifs(mot1,mot2,IC1,IC2,ICthreshold):
    del_inds1 = []
    del_inds2 = []
    #for i in range(0,min(5,mot1.shape[0])):
    for i in range(0,mot1.shape[0]):
        #for j in range(0,min(5,mot2.shape[0])):
        for j in range(0,mot2.shape[0]):
            if overlapcheck(mot1[i,[2,3]],mot2[j,[2,3]]):
                # if we're in here then there is an overlap
                # how long is the overlap?
                # turns out some motifs can entirely eat other motifs, chuck these out first
                if mot1[i,2] - mot2[j,2] > 0 and mot2[j,3] - mot1[i,3] > 0:
                    #mot1 is inside mot2 - chuck it out
                    del_inds1.append(i)
                    del_inds2.append(j)
                elif mot2[j,2] - mot1[i,2] > 0 and mot1[i,3] - mot2[j,3] > 0:
                    #mot2 is inside mot1 - chuck it out
                    del_inds1.append(i)
                    del_inds2.append(j)
                elif mot2[j,2] - mot1[i,2] > 0:
                    # overlap is in the last part of motif1 and beginning of motif2
                    overlaplen = mot1[i,3] - mot2[j,2] + 1
                    # remember to make this an argument later
                    # fwdIC is defined as the sum of the IC of the overlap in the
                    # beginning of the second motif
                    # revIC is defined as the sum of IC in the overlapping end of motif 1
                    if (min(sum(IC1[:overlaplen]),sum(IC2[mot_len1-overlaplen:])) > ICthreshold):
                        # if we're in here it means the overlap is too large and we need to ditch both of these hits
                        del_inds1.append(i)
                        del_inds2.append(j)
                elif mot2[j,2] - mot1[i,2] < 0:
                    # overlap is in the first part of motif1 and end of motif2
                    overlaplen = mot2[j,3] - mot1[i,2] + 1
                    if (min(sum(IC1[:overlaplen]),sum(IC2[mot_len2-overlaplen:])) > ICthreshold):
                        # if we're in here it means the overlap is too large and we need to ditch both of these hits
                        del_inds1.append(i)
                        del_inds2.append(j)
                elif mot2[j,2] - mot1[i,2] ==0:
                    del_inds1.append(i)
                    del_inds2.append(j)
    # del_inds = np.row_stack((del_inds1,del_inds2))
    return(del_inds1,del_inds2)

def find_genes_with_both_motifs(fimo1,fimo2,IC1,IC2,promsize,topNthreshold,mot1_len,mot2_len,intersecting_genes,ICthreshold):
    intersecting_genes = intersect(fimo1[:,1],fimo2[:,1])
    fimo1_1= fimo1[np.nonzero(np.in1d(fimo1[:,1],intersecting_genes))[0],:]
    # do the same for fimo2
    fimo2_1= fimo2[np.nonzero(np.in1d(fimo2[:,1],intersecting_genes))[0],:]
    bm_genes=[]
    iter=enumerate(intersecting_genes)
    for i in iter:
#        print(i[1])
        (mot1,mot2)=subset_fimo(fimo1_1,fimo2_1,i)
        print('mot1 ' + mot1[0,0])
        print('mot2 ' + mot2[0,0])
        (del_inds1,del_inds2) = defineOverlapsBetweenTwoMotifs(mot1,mot2,IC1,IC2,ICthreshold)
        if ((len(del_inds1) & len(del_inds2))==0):
            # we didn't need to remove any overlaps so significance threshold won't change
            # we don't have to compute the binomial score again
            bm_genes.append(mot1[0,1])
        else:
            #delete overlaps that exceed the IC threshold
            mot1 = np.delete(mot1,del_inds1,0)
            mot2 = np.delete(mot2,del_inds2,0)
            # does it still contain both motifs? or have we just deleted all hits?
            if mot1.shape[0] != 0:
                # promoter still contains M1 hits
                if mot2.shape[0] != 0:
                    # promoter still contains M1 and M2 hits
                    # recalculate the binom tests
                    binom_p1 = geom_binom_test(mot1[:,6],promsize[fimo1_1[0,1]],mot1_len)
                    binom_p2 = geom_binom_test(mot2[:,6],promsize[fimo1_1[0,1]],mot2_len)
                    if any(binom_p1 <= topNthreshold[mot1[0,0]]) and any(binom_p2 <= topNthreshold[mot2[0,0]]):
                        bm_genes.append(mot1[0,1])
    return(bm_genes)




def upload_meta_files(indexpath,genefile):
    promsize = {}
    with open((indexpath+'/promoter_lengths.txt'),'r') as fid:
        for line in fid:
            line=line.split()
            promsize[line[0]] = int(line[1])

    topNthreshold = {}
    with open((indexpath+'/binomial_thresholds.txt'),'r') as fid:
        for line in fid:
            line=line.split()
            topNthreshold[line[0]] = float(line[1])

    # load in pickle obkect that contains dictionary of IC vectors for each motif
    with open((indexpath+'/ICdict.pickle'), 'rb') as handle:
        IC_dict = pickle.load(handle)

    input = pd.read_csv(genefile,sep='\t', index_col=None, header=None)
    input = input.values

    # list files in fimohits file
    fimofiles = os.listdir((indexpath+'/fimohits'))
    # remove (invisible) dot files from list
    fimofiles = [f for f in fimofiles if not f[0] == '.']
    return(promsize,topNthreshold,IC_dict,input,fimofiles)

def read_fimo_file(file,indexpath):
    fimo = pd.read_csv((indexpath+'/fimohits/'+file), sep='\t', index_col=None, header=None)
    fimo = fimo.values
    mot_len = fimo[0,3] - fimo[0,2] + 1
    return(fimo,mot_len)



def coloc_test(fimo1,fimo2,motifICstring1,motifICstring2,promsize,topNthreshold,ICthreshold,mot_len1,mot_len2,input,universe):
    ## which promoters contain both motifs
    intersecting_genes = intersect(fimo1[:,1],fimo2[:,1])
    if len(intersecting_genes)>0:
        bm_genes=find_genes_with_both_motifs(fimo1,fimo2,motifICstring1,motifICstring2,promsize,topNthreshold,mot_len1,mot_len2,intersecting_genes,ICthreshold)
    else:
        bm_genes=[]
    ## for each cluster of promoters that you want to test
    clusters=list(set(input[:,0]))
    allresults=np.empty((len(clusters),8),dtype="object")
    a=0
    for cc in clusters:
        promoters_to_test=input[input[:,0]==cc,1]
        # find the overlap between your input and the the list of promoters with both motifs.
        if (len(intersect(promoters_to_test,bm_genes)) > 0):
            pval = np.exp(log_hypergeometric(len(intersect(promoters_to_test,bm_genes)),len(promoters_to_test),len(bm_genes),universe))
            allresults[a,:]=[str(cc),fimo1[0,0],fimo2[0,0],str(pval),len(intersect(promoters_to_test,bm_genes)),len(bm_genes),len(promoters_to_test),';'.join(intersect(promoters_to_test,bm_genes))]
        else:
            allresults[a,:]=[str(cc),fimo1[0,0],fimo2[0,0],str(1),len(intersect(promoters_to_test,bm_genes)),len(bm_genes),len(promoters_to_test),'na']
        a=a+1
    return(allresults)

##################
##################
##################


if __name__ == "__main__":
    args = get_args()
    (promsize,topNthreshold,IC_dict,input,fimofiles)=upload_meta_files(args.pathToIndex,args.testgenesfile)
    total_clusters=int(len(list(set(input[:,0]))))
    for a in range(0,len(fimofiles)):
        # initialise the empty results array
        total_rows=int(len(fimofiles)*total_clusters)
        results=np.empty((total_rows,8),dtype="object")
        i=0
        # load files
        m1=fimofiles[a]
        print(m1)
        (fimo1,mot_len1)=read_fimo_file(m1,args.pathToIndex)
        IC1=IC_dict[fimo1[0,0]]
        for b in range(a+1,len(fimofiles)):
            #load files
            m2=fimofiles[b]
            (fimo2,mot_len2)=read_fimo_file(m2,args.pathToIndex)
            IC2=IC_dict[fimo2[0,0]]
            # do the test
            results[i:i+total_clusters,:]=coloc_test(fimo1,fimo2,IC1,IC2,promsize,topNthreshold,args.ICthreshold,mot_len1,mot_len2,input,args.universe)
            i=i+total_clusters
            del fimo2
        df = pd.DataFrame(results)
        df =df.dropna()
        df.to_csv(('motif_found.txt'),mode='a',sep='\t',header=False,index=False)
