from hw2skeleton import SW
from hw2skeleton import io
import os
import numpy as np

seq1 = 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'
seq2 = 'ANKTRELCMKSLEHAKVDTSNEARQDGIDLYKHMFENYPPLRKYFKSREEYTAEDVQNDPFFAKQGQKILLACHVLCATYDDRETFNAYTRELLDRHARDHVHMPPEVWTDFWKLFEEYLGKKTTLDEPTKQAWHEIGREFAKEINK'
blosum50 = get_sub_mtx(sys.argv[2] + 'BLOSUM50_nc')

def test_score():
    score = get_score_bp(seq1, seq2, blosum50, -5, -1, 'raw')
    assert score == 123

def test_cutoff():
