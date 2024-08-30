import os
import pandas as pd
import numpy as np
import boost_histogram as bh
import pickle
from tqdm import tqdm

ref_csv = pd.read_csv('/home/myang/subnucleosome/Subnucleosome_samples.csv', sep=',')
indices = np.array([62, 63, 64, 65, 124, 125, 126, 208, 209, 210, 211])
hep_ref = ref_csv[ref_csv.INDEX.isin(indices)]

os.chdir('/home/myang/subnucleosome/unsorted_hepatocyte/footprint_centric/')

def eat_pickle(hmmFile):
    with open(hmmFile, 'rb') as fin:
        return pickle.load(fin)

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

    sites_tot = pd.concat([site_plus[site_plus.idx >= 100], site_minus[site_minus.idx >= 100]])
    return sites_tot

for idx, line in hep_ref.iterrows():
    label_out = line['SAMPLE']
    frames = {
        '1': pd.read_table(line['FOOTPRINT_1'], sep=',', header=0),
        '31': pd.read_table(line['FOOTPRINT_31'], sep=',', header=0),
        '51': pd.read_table(line['FOOTPRINT_51'], sep=',', header=0),
        '71': pd.read_table(line['FOOTPRINT_71'], sep=',', header=0),
        '101': pd.read_table(line['FOOTPRINT_101'], sep=',', header=0)
    }

    r0_acc = eat_pickle(line['HMMRES_1'])
    r1_acc = eat_pickle(line['HMMRES_31'])
    r2_acc = eat_pickle(line['HMMRES_51'])
    r3_acc = eat_pickle(line['HMMRES_71'])
    r4_acc = eat_pickle(line['HMMRES_101'])

    ref_tot = frames['1']
    r1_tot = frames['31']
    r2_tot = frames['51']
    r3_tot = frames['71']
    r4_tot = frames['101']
    
    ref_frame = ref_tot[(ref_tot.length < 200) & (ref_tot.midpoint > 200)]
    r0_frame = ref_frame[(ref_frame.length < 80) & (ref_frame.length > 20)]
    r1_frame = r1_tot[(r1_tot.length < 200) & (r1_tot.length > 20)]
    r2_frame = r2_tot[(r2_tot.length < 200) & (r2_tot.length > 20)]
    r3_frame = r3_tot[(r3_tot.length < 200) & (r3_tot.length > 20)]
    r4_frame = r4_tot[(r4_tot.length < 200) & (r4_tot.length > 20)]

    files = [line['CHIP_FOXA2'], line['RANDOM_FOXA2']]
    label_set = ['BOUND_FOXA2', 'RANDOM_FOXA2']

    for lbl_cnt, file in enumerate(files):
        label = label_set[lbl_cnt]
        reference_in = idx_register_motifs(file)
        try:
            overlap_df = ref_tot[ref_tot.length < 200].merge(reference_in, on='zmw')
        except ValueError:
            continue
        overlap_df = overlap_df[np.abs(overlap_df.midpoint - overlap_df.idx) < 500]
        overlap_df['tf'] = label

        zmws = {line['zmw']: (line['tf'], line['mstrand'], line['rstrand']) for idx, line in overlap_df.iterrows()}

        zmw_ids, nuc_ids, label_ids = [], [], []
        nuc_profiles_r0, nuc_profiles_r1, nuc_profiles_r2, nuc_profiles_r3, nuc_profiles_r4 = [], [], [], [], []
        r0_mat, r1_mat, r2_mat, r3_mat, r4_mat = [], [], [], [], []
        nuc_sizes = []
        fiber_counter = 0
        unique_refs = set()

        for i in tqdm(zmws):
            mdir, rdir = zmws[i][1], zmws[i][2]
            ref = int(i)
            nucs = overlap_df[overlap_df.zmw == ref]['midpoint'].values
            if len(nucs) == 0:
                continue
            lengths = overlap_df[overlap_df.zmw == ref]['length'].values

            r0_sub = r0_frame[r0_frame.zmw == ref]
            r1_sub = r1_frame[r1_frame.zmw == ref]
            r2_sub = r2_frame[r2_frame.zmw == ref]
            r3_sub = r3_frame[r3_frame.zmw == ref]
            r4_sub = r4_frame[r4_frame.zmw == ref]

            dist_r0 = np.subtract.outer(r0_sub['midpoint'].values, nucs)
            dist_r0[~((dist_r0 < 100) & (dist_r0 > -100))] = np.nan
            mlen_r0 = np.subtract.outer(r0_sub['length'].values, np.zeros(len(nucs)))
            mlen_r0[~((dist_r0 < 100) & (dist_r0 > -100))] = np.nan    

            dist_r1 = np.subtract.outer(r1_sub['midpoint'].values, nucs)
            dist_r1[~((dist_r1 < 100) & (dist_r1 > -100))] = np.nan
            mlen_r1 = np.subtract.outer(r1_sub['length'].values, np.zeros(len(nucs)))
            mlen_r1[~((dist_r1 < 100) & (dist_r1 > -100))] = np.nan

            dist_r2 = np.subtract.outer(r2_sub['midpoint'].values, nucs)
            dist_r2[~((dist_r2 < 100) & (dist_r2 > -100))] = np.nan
            mlen_r2 = np.subtract.outer(r2_sub['length'].values, np.zeros(len(nucs)))
            mlen_r2[~((dist_r2 < 100) & (dist_r2 > -100))] = np.nan

            dist_r3 = np.subtract.outer(r3_sub['midpoint'].values, nucs)
            dist_r3[~((dist_r3 < 100) & (dist_r3 > -100))] = np.nan
            mlen_r3 = np.subtract.outer(r3_sub['length'].values, np.zeros(len(nucs)))
            mlen_r3[~((dist_r3 < 100) & (dist_r3 > -100))] = np.nan

            dist_r4 = np.subtract.outer(r4_sub['midpoint'].values, nucs)
            dist_r4[~((dist_r4 < 100) & (dist_r4 > -100))] = np.nan
            mlen_r4 = np.subtract.outer(r4_sub['length'].values, np.zeros(len(nucs)))
            mlen_r4[~((dist_r4 < 100) & (dist_r4 > -100))] = np.nan

            for j in range(len(nucs)):
                midpoint = nucs[j]
                size = lengths[j]
                coord1 = int(midpoint - 100)
                try:
                    new_key = int(i)
                    acc = r0_acc[new_key][coord1:coord1 + 200]
                except KeyError:
                    prefix = list(r0_acc.keys())[0].split('/')[0]
                    suffix = 'ccs'
                    new_key = f"{prefix}/{int(i)}/{suffix}"
                    acc = r0_acc[new_key][coord1:coord1 + 200]
                if len(acc) < 200:
                    continue

                zmw_r0 = acc
                zmw_r1 = r1_acc[new_key][coord1:coord1 + 200]
                zmw_r2 = r2_acc[new_key][coord1:coord1 + 200]
                zmw_r3 = r3_acc[new_key][coord1:coord1 + 200]
                zmw_r4 = r4_acc[new_key][coord1:coord1 + 200]

                if mdir == '+':
                    if rdir == "+":
                        r0_mat.append(zmw_r0)
                        r1_mat.append(zmw_r1)
                        r2_mat.append(zmw_r2)
                        r3_mat.append(zmw_r3)
                        r4_mat.append(zmw_r4)

                        nuc_hist_r0 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r1 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r2 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r3 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r4 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))

                        nuc_hist_r0.fill(mlen_r0[:, j], dist_r0[:, j])
                        nuc_hist_r1.fill(mlen_r1[:, j], dist_r1[:, j])
                        nuc_hist_r2.fill(mlen_r2[:, j], dist_r2[:, j])
                        nuc_hist_r3.fill(mlen_r3[:, j], dist_r3[:, j])
                        nuc_hist_r4.fill(mlen_r4[:, j], dist_r4[:, j])

                        vals_r0 = nuc_hist_r0.values().flatten()
                        vals_r1 = nuc_hist_r1.values().flatten()
                        vals_r2 = nuc_hist_r2.values().flatten()
                        vals_r3 = nuc_hist_r3.values().flatten()
                        vals_r4 = nuc_hist_r4.values().flatten()

                        nuc_profiles_r0.append(vals_r0)
                        nuc_profiles_r1.append(vals_r1)
                        nuc_profiles_r2.append(vals_r2)
                        nuc_profiles_r3.append(vals_r3)
                        nuc_profiles_r4.append(vals_r4)

                        zmw_ids.append(i)
                        label_ids.append(zmws[i])
                        nuc_ids.append(nucs[j])
                        nuc_sizes.append(lengths[j])
                    else:
                        r0_mat.append(np.flip(zmw_r0))
                        r1_mat.append(np.flip(zmw_r1))
                        r2_mat.append(np.flip(zmw_r2))
                        r3_mat.append(np.flip(zmw_r3))
                        r4_mat.append(np.flip(zmw_r4))

                        nuc_hist_r0 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r1 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r2 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r3 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r4 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))

                        nuc_hist_r0.fill(mlen_r0[:, j], -1 * dist_r0[:, j])
                        nuc_hist_r1.fill(mlen_r1[:, j], -1 * dist_r1[:, j])
                        nuc_hist_r2.fill(mlen_r2[:, j], -1 * dist_r2[:, j])
                        nuc_hist_r3.fill(mlen_r3[:, j], -1 * dist_r3[:, j])
                        nuc_hist_r4.fill(mlen_r4[:, j], -1 * dist_r4[:, j])

                        vals_r0 = nuc_hist_r0.values().flatten()
                        vals_r1 = nuc_hist_r1.values().flatten()
                        vals_r2 = nuc_hist_r2.values().flatten()
                        vals_r3 = nuc_hist_r3.values().flatten()
                        vals_r4 = nuc_hist_r4.values().flatten()

                        nuc_profiles_r0.append(vals_r0)
                        nuc_profiles_r1.append(vals_r1)
                        nuc_profiles_r2.append(vals_r2)
                        nuc_profiles_r3.append(vals_r3)
                        nuc_profiles_r4.append(vals_r4)

                        zmw_ids.append(i)
                        label_ids.append(zmws[i])
                        nuc_ids.append(nucs[j])
                        nuc_sizes.append(lengths[j])

                elif mdir == '-': 
                    if rdir == "+":
                        r0_mat.append(np.flip(zmw_r0))
                        r1_mat.append(np.flip(zmw_r1))
                        r2_mat.append(np.flip(zmw_r2))
                        r3_mat.append(np.flip(zmw_r3))
                        r4_mat.append(np.flip(zmw_r4))

                        nuc_hist_r0 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r1 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r2 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r3 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r4 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))

                        nuc_hist_r0.fill(mlen_r0[:, j], -1 * dist_r0[:, j])
                        nuc_hist_r1.fill(mlen_r1[:, j], -1 * dist_r1[:, j])
                        nuc_hist_r2.fill(mlen_r2[:, j], -1 * dist_r2[:, j])
                        nuc_hist_r3.fill(mlen_r3[:, j], -1 * dist_r3[:, j])
                        nuc_hist_r4.fill(mlen_r4[:, j], -1 * dist_r4[:, j])

                        vals_r0 = nuc_hist_r0.values().flatten()
                        vals_r1 = nuc_hist_r1.values().flatten()
                        vals_r2 = nuc_hist_r2.values().flatten()
                        vals_r3 = nuc_hist_r3.values().flatten()
                        vals_r4 = nuc_hist_r4.values().flatten()

                        nuc_profiles_r0.append(vals_r0)
                        nuc_profiles_r1.append(vals_r1)
                        nuc_profiles_r2.append(vals_r2)
                        nuc_profiles_r3.append(vals_r3)
                        nuc_profiles_r4.append(vals_r4)

                        zmw_ids.append(i)
                        label_ids.append(zmws[i])
                        nuc_ids.append(nucs[j])
                        nuc_sizes.append(lengths[j])
                    else:
                        r0_mat.append(zmw_r0)
                        r1_mat.append(zmw_r1)
                        r2_mat.append(zmw_r2)
                        r3_mat.append(zmw_r3)
                        r4_mat.append(zmw_r4)

                        nuc_hist_r0 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r1 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r2 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r3 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))
                        nuc_hist_r4 = bh.Histogram(bh.axis.Regular(40, 20, 200), bh.axis.Regular(40, -100, 100))

                        nuc_hist_r0.fill(mlen_r0[:, j], dist_r0[:, j])
                        nuc_hist_r1.fill(mlen_r1[:, j], dist_r1[:, j])
                        nuc_hist_r2.fill(mlen_r2[:, j], dist_r2[:, j])
                        nuc_hist_r3.fill(mlen_r3[:, j], dist_r3[:, j])
                        nuc_hist_r4.fill(mlen_r4[:, j], dist_r4[:, j])

                        vals_r0 = nuc_hist_r0.values().flatten()
                        vals_r1 = nuc_hist_r1.values().flatten()
                        vals_r2 = nuc_hist_r2.values().flatten()
                        vals_r3 = nuc_hist_r3.values().flatten()
                        vals_r4 = nuc_hist_r4.values().flatten()

                        nuc_profiles_r0.append(vals_r0)
                        nuc_profiles_r1.append(vals_r1)
                        nuc_profiles_r2.append(vals_r2)
                        nuc_profiles_r3.append(vals_r3)
                        nuc_profiles_r4.append(vals_r4)

                        zmw_ids.append(i)
                        label_ids.append(zmws[i])
                        nuc_ids.append(nucs[j])
                        nuc_sizes.append(lengths[j])
            if ref not in unique_refs:
                unique_refs.add(ref)
                fiber_counter += 1

        np.save(f'{label_out}_{label}_individual_nucs_r0.npy', np.vstack(nuc_profiles_r0))
        np.save(f'{label_out}_{label}_individual_nucs_r1.npy', np.vstack(nuc_profiles_r1))
        np.save(f'{label_out}_{label}_individual_nucs_r2.npy', np.vstack(nuc_profiles_r2))
        np.save(f'{label_out}_{label}_individual_nucs_r3.npy', np.vstack(nuc_profiles_r3))
        np.save(f'{label_out}_{label}_individual_nucs_r4.npy', np.vstack(nuc_profiles_r4))

        np.save(f'{label_out}_{label}_acc_r0.npy', np.vstack(r0_mat))
        np.save(f'{label_out}_{label}_acc_r1.npy', np.vstack(r1_mat))
        np.save(f'{label_out}_{label}_acc_r2.npy', np.vstack(r2_mat))
        np.save(f'{label_out}_{label}_acc_r3.npy', np.vstack(r3_mat))
        np.save(f'{label_out}_{label}_acc_r4.npy', np.vstack(r4_mat))

        nuc_profiles_r0, nuc_profiles_r1, nuc_profiles_r2, nuc_profiles_r3, nuc_profiles_r4 = [], [], [], [], []

        with open(f'{label_out}_{label}_individual_labels.txt', 'w') as fho:
            for i in range(len(zmw_ids)):
                print(f"{zmw_ids[i]}\t{nuc_ids[i]}\t{nuc_sizes[i]}\t{label_ids[i]}", file=fho)
