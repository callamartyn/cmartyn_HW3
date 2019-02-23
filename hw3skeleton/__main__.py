import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib import pyplot as plt
import sys
from .io import get_test_pairs, get_sub_mtx
from .SW import fp_rate, find_penalties
from .optimize import mutate


print('running main script')
data_dir = sys.argv[2]

# import the pairs of positive and negative sequences
pospairs, negpairs = get_test_pairs(sys.argv[2])

# import the matrices, I didn't know how to skip the comments so I just removed them
print('importing substitution matrices')
blosum50 = get_sub_mtx(sys.argv[2] + 'BLOSUM50_nc')
blosum62 = get_sub_mtx(sys.argv[2] + 'BLOSUM62_nc')
MATIO = get_sub_mtx(sys.argv[2] + 'MATIO_nc')


AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
       'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']

# here I will vary the gap penalities and save the combination that performs the best
best_rate, best_open, best_extend = find_penalties(pospairs, negpairs, blosum50, -20, -5, 'raw')

tp_rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
fp_rates_bl50 = roc_values(pospairs, negpairs, tp_rates, blosum50, 'raw')
fp_rates_bl62 = roc_values(pospairs, negpairs, tp_rates, blosum62, 'raw')

# plotting to visualize the roc curves
fig1 = plt.figure(dpi = 300)
plt.title('Receiver Operating Characteristic')
plt.plot(fp_rates_bl50, tp_rates, 'b', label = 'BLOSUM50')
plt.plot(fp_rates_bl62, tp_rates, 'g', label = 'BLOSUM62')
plt.plot([0, 1], [0, 1],'r--', label = 'Random')
plt.legend(loc = 'lower right')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
fig1.savefig('matrix_comparison.png')

# because blosum50 performed better, I will use that one to check the normalized scores
fp_rates_bl50_norm = roc_values(pospairs, negpairs, tp_rates, blosum50, 'normalized')

# and now I will plot that compared to the regular scores
fig2 = plt.figure(dpi = 300)
plt.title('Receiver Operating Characteristic Normalized')
plt.plot(fp_rates_bl50, tp_rates, 'b', label = 'Raw Score')
plt.plot(fp_rates_bl50_norm, tp_rates, 'g', label = 'Normalized Score')
plt.plot([0, 1], [0, 1],'r--', label = 'Random')
plt.legend(loc = 'lower right')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
fig2.savefig('normalize_comparison.png')

# now I will try to optimize the blosum50 matrix using the raw scores
optimize(blosum50, 1.9)
