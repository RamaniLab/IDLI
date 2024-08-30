import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import boost_histogram as bh
import pickle

def convert_keys_for_tn5(r0_acc):
    converted_acc = {}
    for key in r0_acc.keys():
        if '/' in key:
            middle_key = key.split('/')[1]
            converted_acc[int(middle_key)] = r0_acc[key]
        else:
            converted_acc[int(key)] = r0_acc[key]
    return converted_acc

def process_accessibility(arr):
    arr = np.where(arr > 0.5, 1, arr)
    arr = np.where(arr <= 0.5, 0, arr)
    arr = np.where(np.isnan(arr), 2, arr)
    return arr

def get_central_idx(zmw, idx_list):
    length = len(r0_acc.get(zmw, []))
    central_pos = length / 2
    distances = np.abs(idx_list - central_pos)
    return idx_list[np.argmin(distances)]

def filter_central_most_idx(group):
    central_idx = get_central_idx(group['zmw'].iloc[0], group['idx'].values)
    return group[group['idx'] == central_idx]

def idx_register_motifs(motif_file):
    sitefile = pd.read_table(motif_file, sep='\t')
    if len(sitefile.columns) == 8:
        names = ['zmw', 'chrid', 'rstart', 'rend', 'site', 'score', 'rstrand', 'mstrand']
    elif len(sitefile.columns) == 7:
        names = ['zmw', 'chrid', 'rstart', 'rend', 'site', 'rstrand', 'mstrand']
    else:
        names = ['zmw', 'chrid', 'rstart', 'rend', 'rstrand', 'mstrand']
    sitefile.columns = names

    site_plus = sitefile[sitefile.rstrand == '+']
    site_minus = sitefile[sitefile.rstrand == '-']

    site_plus['idx'] = site_plus['site'] - site_plus['rstart']
    site_minus['idx'] = -1 * (site_minus['site'] - site_minus['rend'])

    sites_tot = pd.concat([site_plus, site_minus])

    def filter_groups(group):
        if group.duplicated(subset=['chrid', 'rstart', 'rend']).any():
            return group
        elif group.shape[0] == 1:
            return group
        else:
            return pd.DataFrame()

    filtered_sites = sites_tot.groupby('zmw').apply(filter_groups).reset_index(drop=True)
    filtered_sites = filtered_sites[(filtered_sites['idx'] >= 100)]
    return filtered_sites

def eat_pickle(hmmFile):
    with open(hmmFile, 'rb') as fin:
        return pickle.load(fin)

ref_csv = pd.read_csv('/home/myang/subnucleosome/Subnucleosome_samples.csv', sep=',')
indices = np.array([62, 63, 64, 65, 124, 125, 126, 208, 209, 210, 211])
hep_ref = ref_csv[ref_csv.INDEX.isin(indices)]

os.chdir('/home/myang/subnucleosome/unsorted_hepatocyte/footprint_centric/')
wd = '/home/myang/subnucleosome/unsorted_hepatocyte/footprint_centric/'
domains = ['BOUND_FOXA2', 'RANDOM_FOXA2']

nuc1_midpoints, nuc2_midpoints, nuc3_midpoints, nuc4_midpoints, nuc5_midpoints = [], [], [], [], []
nuc1_lengths, nuc2_lengths, nuc3_lengths, nuc4_lengths, nuc5_lengths = [], [], [], [], []
domains_list, samples_zmws, feature_idxs, accessibility_patterns = [], [], [], []
mstrands, rstrands, genotypes, reps = [], [], [], []

for idx, line in tqdm(hep_ref.iterrows(), total=hep_ref.shape[0], desc="Processing Samples"):
    label = line['SAMPLE']
    feature_boundfoxa2 = line['CHIP_FOXA2']
    feature_randomfoxa2 = line['RANDOM_FOXA2']
    geno = line['HEP_SAMPLE']
    rep = line['HEP_REP']
    hmmFile = line['HMMRES_1']
    r0_acc = eat_pickle(hmmFile)
    if 'Tn5' in label:
        r0_acc = convert_keys_for_tn5(r0_acc)

    files = [feature_boundfoxa2, feature_randomfoxa2]
    skip_sample = False
    for domain in domains:
        file_path = f'{wd}{label}_{domain}_individual_labels.txt'
        if not os.path.exists(file_path):
            skip_sample = True
            break
    if skip_sample:
        continue

    lbl_cnt = 0
    for file in files:
        try:
            domain = domains[lbl_cnt]
            lengths = pd.read_table(f'{wd}{label}_{domain}_individual_labels.txt', sep='\t', header=None)
            lengths.columns = ['zmw', 'midpoint', 'length', 'info']
            lengths[['site_type', 'mstrand', 'rstrand']] = lengths['info'].str.extract(r"\('(\w+)', '(\+|-)', '(\+|-)'\)")
            lengths_unique = lengths.drop_duplicates()
            lengths_unique['unique_id'] = lengths_unique['zmw'].apply(lambda x: f"{label}_{x}")
            lengths_unique['unique_id'] = lengths_unique['unique_id'].astype(str)
            lbl_cnt += 1
            reference_in = idx_register_motifs(file)
            reference_in['zmw'] = reference_in['zmw'].astype(int)
            reference_in['unique_id'] = reference_in['zmw'].apply(lambda x: f"{label}_{x}")

            for idx, row in tqdm(reference_in.iterrows(), total=reference_in.shape[0]):
                unique_id = str(row['unique_id'])
                idx_value = int(row['idx'])
                zmw_of_interest = row['zmw']
                mstrand = str(row['mstrand'])
                rstrand = str(row['rstrand'])

                matched_lengths = lengths_unique[lengths_unique['unique_id'] == unique_id].copy()
                if len(matched_lengths) < 3:
                    continue
                matched_lengths['idx_value'] = idx_value
                matched_lengths['abs_distance'] = np.abs(matched_lengths['midpoint'] - idx_value)

                if (matched_lengths['abs_distance'] < 500).sum() < 3:
                    continue
                
                matched_lengths = matched_lengths.sort_values(by='abs_distance')
                nuc3_midpoint = matched_lengths.iloc[0]['midpoint']
                nuc3_length = matched_lengths.iloc[0]['length']

                nuc2_candidates = matched_lengths[(matched_lengths['midpoint'] < idx_value) & (matched_lengths['midpoint'] != nuc3_midpoint)]
                if len(nuc2_candidates) > 0:
                    nuc2_midpoint = nuc2_candidates.iloc[0]['midpoint']
                    nuc2_length = nuc2_candidates.iloc[0]['length']
                else:
                    nuc2_midpoint, nuc2_length = None, None

                nuc4_candidates = matched_lengths[(matched_lengths['midpoint'] > idx_value) & (matched_lengths['midpoint'] != nuc3_midpoint)]
                if len(nuc4_candidates) > 0:
                    nuc4_midpoint = nuc4_candidates.iloc[0]['midpoint']
                    nuc4_length = nuc4_candidates.iloc[0]['length']
                else:
                    nuc4_midpoint, nuc4_length = None, None

                if nuc2_midpoint is not None:
                    nuc1_candidates = matched_lengths[(matched_lengths['midpoint'] < nuc2_midpoint)]
                    if len(nuc1_candidates) > 0:
                        nuc1_midpoint = nuc1_candidates.iloc[0]['midpoint']
                        nuc1_length = nuc1_candidates.iloc[0]['length']
                    else:
                        nuc1_midpoint, nuc1_length = None, None
                else:
                    nuc1_midpoint, nuc1_length = None, None

                if nuc4_midpoint is not None:
                    nuc5_candidates = matched_lengths[(matched_lengths['midpoint'] > nuc4_midpoint)]
                    if len(nuc5_candidates) > 0:
                        nuc5_midpoint = nuc5_candidates.iloc[0]['midpoint']
                        nuc5_length = nuc5_candidates.iloc[0]['length']
                    else:
                        nuc5_midpoint, nuc5_length = None, None
                else:
                    nuc5_midpoint, nuc5_length = None, None

                nuc1_midpoints.append(nuc1_midpoint)
                nuc2_midpoints.append(nuc2_midpoint)
                nuc3_midpoints.append(nuc3_midpoint)
                nuc4_midpoints.append(nuc4_midpoint)
                nuc5_midpoints.append(nuc5_midpoint)

                nuc1_lengths.append(nuc1_length)
                nuc2_lengths.append(nuc2_length)
                nuc3_lengths.append(nuc3_length)
                nuc4_lengths.append(nuc4_length)
                nuc5_lengths.append(nuc5_length)
                
                start = int(idx_value) - 1000
                end = int(idx_value) + 1000
        
                acc = r0_acc[int(zmw_of_interest)][max(0, start):min(len(r0_acc[int(zmw_of_interest)]), end)]
        
                if start < 0:
                    acc = np.pad(acc, (abs(start), 0), constant_values=np.nan)
                if end > len(r0_acc[int(zmw_of_interest)]):
                    acc = np.pad(acc, (0, end - len(r0_acc[int(zmw_of_interest)])), constant_values=np.nan)
        
                assert len(acc) == 2000, f"Length of acc is {len(acc)}, expected 2000"
        
                acc_processed = process_accessibility(acc)
                
                domains_list.append(domain)
                samples_zmws.append(unique_id)
                feature_idxs.append(idx_value)
                accessibility_patterns.append(acc_processed)
                mstrands.append(mstrand)
                rstrands.append(rstrand)
                genotypes.append(geno)
                reps.append(rep)

        except FileNotFoundError:
            continue

        except Exception as e:
            print(f"Error processing {label} for domain {domain}: {e}")

np.save('nuc1_midpoints.npy', nuc1_midpoints)
np.save('nuc2_midpoints.npy', nuc2_midpoints)
np.save('nuc3_midpoints.npy', nuc3_midpoints)
np.save('nuc4_midpoints.npy', nuc4_midpoints)
np.save('nuc5_midpoints.npy', nuc5_midpoints)

np.save('nuc1_lengths.npy', nuc1_lengths)
np.save('nuc2_lengths.npy', nuc2_lengths)
np.save('nuc3_lengths.npy', nuc3_lengths)
np.save('nuc4_lengths.npy', nuc4_lengths)
np.save('nuc5_lengths.npy', nuc5_lengths)

np.save('domains_list.npy', domains_list)
np.save('samples_zmws.npy', samples_zmws)
np.save('feature_idxs.npy', feature_idxs)
np.save('accessibility_patterns.npy', accessibility_patterns)
np.save('mstrands.npy', mstrands)
np.save('rstrands.npy', rstrands)
np.save('genotypes.npy', genotypes)
np.save('reps.npy', reps)
