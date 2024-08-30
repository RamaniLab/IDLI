from __future__ import print_function, division
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

ref_csv = pd.read_csv('/home/myang/subnucleosome/Subnucleosome_samples.csv')

def eat_pickle(hmmFile):
    with open(hmmFile, 'rb') as fin:
        hmmdat = pickle.load(fin)
    return hmmdat

monodi_ratio = []

for idx, line in tqdm(ref_csv.iterrows()):
    foot_1 = line['FOOTPRINT_1']
    frames['1'] = pd.read_table(foot_1, sep=',', header=0)
    ref_tot = frames['1']['length']

    hist, bin_edges = np.histogram(ref_tot, bins=range(min(ref_tot), max(ref_tot) + 1), density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    peaks, _ = find_peaks(hist)

    peaks_in_range_120_180 = [p for p in peaks if 120 <= bin_centers[p] <= 180]
    peaks_in_range_300_360 = [p for p in peaks if 300 <= bin_centers[p] <= 360]

    peak_height_120_180 = max(hist[peaks_in_range_120_180]) if peaks_in_range_120_180 else 0
    peak_height_300_360 = max(hist[peaks_in_range_300_360]) if peaks_in_range_300_360 else 0

    ratio = peak_height_120_180 / peak_height_300_360 if peak_height_300_360 > 0 else np.inf

    monodi_ratio.append(ratio)
    print(ratio)
