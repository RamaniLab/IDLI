import os
import sys
import pandas as pd
import numpy as np
import scipy as sp
from tqdm import tqdm
import boost_histogram as bh
import pickle

target = 0.9 * len(clusters)
cusum = 0

fho = open('/home/myang/subnucleosome/unsorted_hepatocyte/footprint_centric/plotme_vplots_clusters_res_60_foxa2ko_footprintcentric.txt', 'w')
fho2 = open('/home/myang/subnucleosome/unsorted_hepatocyte/footprint_centric/plotme_accplots_clusters_res_60_foxa2ko_footprintcentric.txt', 'w')
fho3 = open('/home/myang/subnucleosome/unsorted_hepatocyte/footprint_centric/plotme_lengths_clusters_res_60_foxa2ko_footprintcentric.txt', 'w')

for i in range(143):
    cusum += len(clusters[clusters == str(i)])
    if cusum >= target:
        break
    print(i)
    print(len(clusters[clusters == str(i)]))
    sub = maximum[clusters == str(i)]
    summed = np.sum(sub, axis=0)
    reshaped = summed.reshape(40, 40)
    plotter = sp.stats.zscore(reshaped[::-1], axis=1)
    acc = try_rplot[clusters == str(i)][:, :140]
    acc2 = try_rplot[clusters == str(i)][:, 140:280]
    acc3 = try_rplot[clusters == str(i)][:, 280:]
    acc_mean = np.nanmean(acc, axis=0)
    acc2_mean = np.nanmean(acc2, axis=0)
    acc3_mean = np.nanmean(acc3, axis=0)
    sublens = tot_lens[clusters == str(i)]

    for j in range(21, 200):
        for k in range(-100, 100):
            ci = int((200 - j) * 40 / 180)
            cj = int((k + 100) / 5)
            print("%s\t%s\t%s\t%s" % (j, k, plotter[ci][cj], i+1), file=fho)
    
    for j in range(140):
        x_coord = j - 70
        y_coord1 = acc_mean[j]
        y_coord2 = acc2_mean[j]
        y_coord3 = acc3_mean[j]
        print("%s\t%s\t%s\t31" % (x_coord, y_coord1, i+1), file=fho2)
        print("%s\t%s\t%s\t51" % (x_coord, y_coord2, i+1), file=fho2)
        print("%s\t%s\t%s\t71" % (x_coord, y_coord3, i+1), file=fho2)
    
    for j in range(len(sublens)):
        print("%s\t%s" % (sublens[j], i+1), file=fho3)

fho.close()
fho2.close()
fho3.close()
