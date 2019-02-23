import numpy as np
from .SW import  get_score_bp

# this function will take an imput matrix and determine the sum of the true positives
# rates for false positive rates of 0, 0.1, 0.2, and 0.3
def evaluate(sub_mtx):
    neg_scores = []
    tp = 0
    fn = 0
    rate_sum = 0
    for i in range(len(negpairs)):
        # get the alignment score of each pair
        align_score = get_score_bp(negpairs[i][0], negpairs[i][1], sub_mtx, best_open, best_extend, 'raw')
        neg_scores.append(align_score)
    for n in [0, 0.1, 0.2, 0.3]:
        percent = (1 - n) * 100
        cutoff = np.percentile(neg_scores, percent)
        for i in range(len(pospairs)):
            align_score = get_score_bp(pospairs[i][0], pospairs[i][1], sub_mtx, best_open, best_extend, 'raw')
        # compare it to the cutoff and add one to the appropriate bucket
            if align_score >= cutoff:
                tp += 1
            else:
                fn += 1
        tp_rate = tp / (tp + fn)
        rate_sum += tp_rate
    return rate_sum

# I'm using a genetic algorithm to optimize my matrix
# This function will randomly change some values and see if the score has improved
def mutate(matrix, goal):
    # use our input matrix to benchmark, if we don't do better than this, just return
    # the starting matrix
    best_score = evaluate(matrix)
    print('best_score', best_score)
    best_matrix = matrix
    # will randomly change some values by this amount
    values = [-2, -1, 1, 2]
    # make a copy of the matrix
    temp_matrix = matrix.copy()
    for j in range (50):
        # randomly sample 2 amino acids and a value
        aa1, aa2 = np.random.choice(AAs, 2)
        change = np.random.choice(values)
        # change both occurences of that combination of amino acids
        temp_matrix[aa1, aa2] = matrix[aa1, aa2] + change
        temp_matrix[aa2, aa1] = matrix[aa2, aa1] + change
    score = evaluate(temp_matrix)
    print('score', score)
    if score > best_score:
        best_score = score
        best_matrix = temp_matrix
    done = False
    if best_score >= goal:
        done = True
    print(done)
    while not done:
        new_mtx, new_score = mutate(best_matrix, goal)
        print('new_score', new_score)
    return best_matrix, best_score
    
