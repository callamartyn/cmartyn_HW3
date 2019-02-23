import os
from Bio import SeqIO
import pandas as pd

def get_sequence(directory, filename):
    filepath = directory + filename
    sequence = str(SeqIO.read(filepath, "fasta").seq.upper())
    return sequence

def get_test_pairs(directory):
    with open(directory + 'Pospairs.txt') as pos:
       positives =  pos.readlines()
    with open(directory + 'Negpairs.txt') as neg:
       negatives =  neg.readlines()
    positives = [x.rstrip() for x in positives]
    positives = [x.split(' ') for x in positives]
    negatives = [x.rstrip() for x in negatives]
    negatives = [x.split(' ') for x in negatives]
    for pair in positives:
        pair[0] = get_sequence(directory, pair[0])
        pair[1] = get_sequence(directory, pair[1])
    for pair in negatives:
        pair[0] = get_sequence(directory, pair[0])
        pair[1] = get_sequence(directory, pair[1])
    return positives, negatives

def get_sub_mtx(path):
    df = pd.read_fwf(path)
    AAs = df.columns.values
    df.set_index(AAs, inplace = True)
    sub_dict = {}
    for i in range(len(df.index)):
        for j in range(len(df.columns)):
            aa1 = df.index[i]
            aa2 = df.columns[j]
            score = df.iloc[i,j]
            sub_dict[(aa1, aa2)] = score
            sub_dict[(aa2, aa1)] = score
    return sub_dict
