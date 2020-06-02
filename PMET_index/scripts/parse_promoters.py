import numpy as np
import pandas as pd
import argparse
from os.path import basename

def get_args():
    #get the arguments
    parser = argparse.ArgumentParser()
    #input file
    parser.add_argument('fasta_path',type=str)
    parser.add_argument('bed_path',type=str)
    args = parser.parse_args()
    return args

args = get_args()

reader = open(args.fasta_path,'r')
bedfile = pd.read_csv(args.bed_path, sep='\t', index_col=None, header=None)

bedfile = bedfile.values
bedfile = bedfile[:,3]
bedfile = np.asarray(['>'+bed+'\n' for bed in bedfile])
fafile = np.asarray(reader.readlines())
reader.close()
fafile[np.arange(0,len(fafile),2)] = bedfile
writer = open('promoters.fa','w')
writer.writelines(fafile)
writer.close()
