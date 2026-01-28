import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
from numba import njit
from numba.typed import Dict
from numba.core import types
import glob
import pysam
from matplotlib import pyplot as plt
import os, sys, re

def HMMv2(hmmInput, t=1/1000):
    d = np.zeros((2, 35000))
    d[0, :] = (1 - t) ** np.arange(35000)
    d[1, :] = 1 - d[0, :]
    d = np.log(d)
    
    binarized_molecules = {}

    for k, v in tqdm(hmmInput.items()):
        seqlen = v['inDat'].shape[0]
        binarized_molecules[k] = np.full(
            shape=(v['cclen']),
            fill_value=np.nan,
            dtype=float)

        tb = np.zeros((seqlen, 2), dtype=float)
        tb_ptr = np.zeros((seqlen, 2), dtype=int)

        _viterbi_v2(
            v['inDat'].pos.values,
            v['inDat'].methPred.values,
            v['inDat'].posProb.values,
            v['inDat'].negProb.values,
            binarized_molecules[k],
            tb,
            tb_ptr,
            d
        )
        
    return binarized_molecules

@njit
def _viterbi_v2(pos, sequence, P_M_Ac, P_M_InAc, path, tb, tb_ptr, d):
    seqlen = sequence.shape[0]
    InAc_vals = np.zeros(2)
    Ac_vals = np.zeros(2)
    
    c = sequence[0]
    log_05 = np.log(0.5)
    tb[0, 0] = log_05 + np.log(1 - c + (2 * c - 1) * P_M_InAc[0])
    tb[0, 1] = log_05 + np.log(1 - c + (2 * c - 1) * P_M_Ac[0])

    for i in range(1, seqlen):
        dist = pos[i] - pos[i - 1]
        dist_stay = d[0, dist]
        dist_leave = d[1, dist]

        c = sequence[i]
        l1 = 1 - c
        l2 = 2 * c - 1
        
        em_ac = np.log(l1 + l2 * P_M_Ac[i])
        em_inac = np.log(l1 + l2 * P_M_InAc[i])

        InAc_vals[0] = tb[i - 1, 0] + dist_stay + em_inac
        InAc_vals[1] = tb[i - 1, 1] + dist_leave + em_inac
        Ac_vals[0] = tb[i - 1, 0] + dist_leave + em_ac
        Ac_vals[1] = tb[i - 1, 1] + dist_stay + em_ac

        if InAc_vals[0] > InAc_vals[1]:
            tb_ptr[i, 0] = 0
            tb[i, 0] = InAc_vals[0]
        else:
            tb_ptr[i, 0] = 1
            tb[i, 0] = InAc_vals[1]

        if Ac_vals[0] > Ac_vals[1]:
            tb_ptr[i, 1] = 0
            tb[i, 1] = Ac_vals[0]
        else:
            tb_ptr[i, 1] = 1
            tb[i, 1] = Ac_vals[1]

    k = 0 if tb[-1, 0] > tb[-1, 1] else 1
    
    for o in range(seqlen - 1, -1, -1):
        curpos = pos[o]
        nextpos = pos[o - 1]
        path[curpos] = k
        k = tb_ptr[o, k]
        
        dist = curpos - nextpos
        k_dist = path[curpos] - k
        slope = (k_dist) / (dist)
        
        for z in range(curpos - 1, nextpos, -1): 
            path[z] = slope * (z - nextpos) + k
            
def call_footprints(hmmFile, outfile):
    with open(hmmFile, 'rb') as fin:
        hmmdat = pickle.load(fin)

    regions = {'zmw': [], 'length': [], 'start': [], 'end': []}

    for zmw in hmmdat:
        hmm = hmmdat[zmw]
        inacregion = hmm
        inacregion[np.isfinite(inacregion)] = inacregion[np.isfinite(inacregion)] > 0.5

        inacswitch = np.diff(inacregion)
        switchp = np.where(np.logical_or(inacswitch == 1, inacswitch == -1))[0]

        if len(switchp) < 1:
            if len(hmm) > 50 and hmm[50] == 0:
                regions['zmw'].append(zmw)
                regions['length'].append(len(hmm) - 0)
                regions['start'].append(np.nan)
                regions['end'].append(np.nan)
            continue
        if inacswitch[switchp[0]] == -1:
            inInacReg = False
            regStart = -1
            regEnd = -1
        if inacswitch[switchp[0]] == 1:
            inInacReg = True
            regStart = np.nan
            regEnd = -1
        for point in switchp:
            if inacswitch[point] == -1 and not inInacReg:
                inInacReg = True
                regStart = point + 1
            if inacswitch[point] == 1 and inInacReg:
                inInacReg = False
                regEnd = point
                regions['zmw'].append(zmw)
                if np.isnan(regStart):
                    regions['length'].append(regEnd - 0)
                else:
                    regions['length'].append(regEnd - regStart)
                regions['start'].append(regStart)
                regions['end'].append(regEnd)
        if inInacReg:
            regions['zmw'].append(zmw)
            regions['length'].append(len(hmm) - regStart)
            regions['start'].append(regStart)
            regions['end'].append(np.nan)
    regionD = pd.DataFrame(regions)
    regionD['start'] = regionD['start'].fillna(0)
    regionD['end'] = regionD['end'].fillna(regionD['start'] + regionD['length'])
    regionD['midpoint'] = (regionD['end'] + regionD['start']) / 2
    regionD.to_csv(outfile)

def main():
    # Usage: script.py your_metadata.csv
    csv_path = sys.argv[1] 
    ref_df = pd.read_csv(csv_path)
    
    t_range = [10, 310, 510, 710, 1010]

    # We use .unique() so we don't process the same sample twice
    for hmm_input_path in ref_df['FOR_HMM'].unique():
        if not os.path.exists(hmm_input_path):
            print(f"File not found: {hmm_input_path}")
            continue

        # Create a base name for output files based on the input filename
        base_output = hmm_input_path.replace('_forHMM.pickle', '')

        print(f"Processing: {hmm_input_path}")

        for i in t_range:
            hmm_fho = f"{base_output}_HMMres_{i}.pickle"
            foot_out = f"{base_output}_{i}_footprint_out.csv"

            with open(hmm_input_path, 'rb') as fin:
                hmmInput = pickle.load(fin)

            hmm_out = HMMv2(hmmInput, t=i / 10000)

            with open(hmm_fho, 'wb') as fout:
                pickle.dump(hmm_out, fout)

            call_footprints(hmm_fho, foot_out)
            
            print(f"   - Saved {i} sensitivity results.")

if __name__ == "__main__":
    main()