import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib import pyplot as plt
from Bio import pairwise2

# this function will create a scoring matrix which can be used to determine the alignment score
# it also returns two other matrices that in conjunction can be used to get the alignment
def mk_weighted_mtx(seq1, seq2, sub_mtx, g_open, g_extend):
    gap_score = {'m': g_open, 'g':g_extend, '':0} # if the previous value was a match, there is a new gap with penalty
    # set dimension for matrices
    m = len(seq1)
    n = len(seq2)
    c = np.zeros((m+1,n+1), dtype=(int)) # c will be the score matrix
    b = np.zeros((m,n), dtype=(int, 2)) # b will keep track of the direction (for traceback)
    d = np.empty((m,n), dtype=str) # d will keep track of whether we assigned a gap or match/mismatch
    for i in range(m):
        for j in range(n):
            match = sub_mtx.loc[seq1[i]].loc[seq2[j]] # get the match/mismatch score from the matrix
            g1 = gap_score[d[i-1, j]]
            g2 = gap_score[d[i, j-1]]
            # find the max value of the possible gaps or mismatches
            possible_scores = (c[i-1, j-1] + match, c[i-1,j] - g1, c[i,j-1] - g2)
            # add this value to score matrix
            c[i,j], max_index = max(possible_scores), np.argmax(possible_scores)
            # now figure out which one was the max and fill in the direction and gap matrices
            if max_index == 0:
                b[i,j] = [i-1,j-1]
                d[i,j] = 'm'
            elif max_index == 1:
                b[i,j] = [i-1,j]
                d[i,j] = 'g'
            else:
                b[i,j] = [i,j-1]
                d[i,j] = 'g'
            # if we've gone negative, reset to 0 (local alignment)
            if c[i,j] < 0:
                c[i,j] = 0
    return b, c, d

# this function will use the matrices created to find the optimal alignment
def find_align(b, c, d, seq1, seq2):
    align1 = []
    align2 = []
    # need to get the coordinates of the last occurence of the highest value in c
    # flip the matrix (so we can look for the first occurence of the highest value)
    c_flip = np.flip(np.flip(c, 0), 1)
    # get coordinates of highest value
    flip_coords = np.unravel_index(c_flip.argmax(), c_flip.shape)
    # convert into coordinates of original matrix
    last = np.subtract(c_flip.shape, np.add(flip_coords, (1, 1)))
    # set starting position as last, will begin walk through matrix here
    pos = last
    next_p = b[tuple(pos)]
    done = False
    while not done:
        next_p = b[tuple(pos)]    # from current position get the direction from matrix b
        if np.array_equal(np.subtract(pos, 1), next_p): # if you moved at an angle to the left (e.g. a match)
            align1.append(seq1[pos[0]]) # get the base and add it to each alignment
            align2.append(seq2[pos[1]])
        elif pos[0] == next_p[0]:
            align1.append('-')
            align2.append(seq2[pos[1]])
        # else add base to seq1 and - to seq2
        else:
            align1.append(seq1[pos[0]])
            align2.append('-')
        pos = next_p    # now move in the direction we were pointing
        done = (c[tuple(pos)] == 0)
    align1.reverse()
    align2.reverse()
    return ''.join(align1), ''.join(align2)

# this function returns only the alignment score for a pair of sequences
# for use in part 1 when we don't need the alignment

def get_score(x, y, sub_mtx, g_open, g_extend, score_type):
    # get only the score matrix
    c = mk_weighted_mtx(x, y, sub_mtx, g_open, g_extend)[1]
    # choose raw or normalized score
    if score_type == 'normalized':
        align_score = np.amax(c) / (min(len(x), len(y)))
    else:
        align_score = np.amax(c)
    return align_score

# this wrapper function will return the alignment given two sequences

def get_align(x, y, sub_mtx, g_open, g_extend):
    # get all the SW matrices
    b, c, d = mk_weighted_mtx(x, y, sub_mtx, g_open, g_extend)
    # find the alignment
    align1, align2 = find_align(b, c, d, x, y)
    return align1, align2

# none of the above functions will be used in main because they run too slowly
# instead I am using biopython local alignment from here on out

# this function will get the alignment score from the biopython alignment
def get_score_bp(a, b, sub_mtx, g_open, g_extend, score_type):
    score = pairwise2.align.localds(a, b, sub_mtx, g_open, g_extend, one_alignment_only = True, score_only = True)
    if score_type == 'normalized':
        return score / min(len(a), len(b))
    else:
        return score

# this function will return the alignment score threshold associated with a given
# true positive rate and gap penalties
def set_cutoff(positives, tp_rate, sub_mtx, g_open, g_extend, score_type):
    # first run through all the positive pairs and record their score
    pos_scores = []
    for i in range(len(positives)):
        score = get_score_bp(positives[i][0], positives[i][1], sub_mtx, g_open, g_extend, score_type)
        pos_scores.append(score)
    # convert the tp rate to a percentile
    # e.g. if we want a 0.7 tp rate, what score represents the 30th percentile
    percent = (1 - tp_rate) * 100
    cutoff = np.percentile(pos_scores, percent)
    return cutoff, pos_scores

# given an alignment score threshold, this function returns the associated
# false positive rate
def fp_rate(negatives, cutoff, sub_mtx, g_open, g_extend, score_type):
    # start with true negatives and false positives at 0
    tn = 0
    fp = 0
    for i in range(len(negatives)):
        # get the alignment score of each pair
        align_score = get_score_bp(negatives[i][0], negatives[i][1], sub_mtx, g_open, g_extend, score_type)
        # compare it to the cutoff and add one to the appropriate bucket
        if align_score >= cutoff:
            fp += 1
        else:
            tn += 1
    rate = fp / (fp + tn)
    return rate

# this function will iterate through a range of penalties for a constant true
# positive rate and return the combination of penalities that gives the lowest
# false positive rate
def find_penalties(positives, negatives, sub_mtx, max_open, max_extend, score_type):
    # this is where I keep track of the best rate and associated penalties
    best_rate = 1
    open_penalty = 0
    extend_penalty = 0
    for i in range(max_open, 0):
        for j in range(max_extend, 0):
            print(i, j)
            # the biopython implementation won't run if the extension penalty is greater
            # than the open penalty
            if i <= j:
                print('calculating true positive threshold...')
                # for this gap combination, find the score cutoff that represents a true positive rate of 70%
                cutoff = set_cutoff(positives, 0.7, sub_mtx, i, j, score_type)[0]
                print ('calculating false positive rate...')
                # for this gap combination and tp rate, find the false positive rate
                rate = fp_rate(negatives, cutoff, sub_mtx, i, j, score_type)
                # if this is better than the previous best rate, replace it along with the corresponding
                # gap penalty combination
                if rate < best_rate:
                    best_rate = rate
                    open_penalty = i
                    extend_penalty = j
    return best_rate, open_penalty, extend_penalty

# this function will return the false positive rates associated with a list of
# true positive rates, these lists can be use to plot an roc curve
def roc_values(positives, negatives, tp_rates, sub_mtx, score_type):
    fp_rates = []
    # get all the positive scores
    pos_scores = set_cutoff(positives, 0.7, sub_mtx, best_open, best_extend, score_type)[1]
    # get all the cutoffs for the range of tp rates
    for n in tp_rates:
        percent = (1 - n) * 100
        cutoff = np.percentile(pos_scores, percent)
        rate = fp_rate(negatives, cutoff, sub_mtx, best_open, best_extend, score_type)
        fp_rates.append(rate)
    return fp_rates
