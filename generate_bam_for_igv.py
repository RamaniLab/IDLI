import sys
import os
import numpy as np 
import pandas as pd
import pickle
from tqdm import tqdm
from numba import njit 
from numba.typed import Dict
from numba.core import types
import pysam

cellResults = "/wynton/group/goodarzilab/ramanilab/results/pacbio/240402_MY_BrdUIP_FOXA2_CsLow/"
samples = [2, 3]
merged_bam = cellResults + "aligned/240402_MY_BrdUIP_FOXA2_CsLow.split.FOXA2_F1_Het.ccs.aligned.sorted.bam"
referenceFile = "//wynton/home/ramani/myang/pipeline/sampleReferences/240402_MY_BrdUIP_FOXA2_CsLow.sampleReference.wynton.csv"
output_directory = cellResults

args = (samples, merged_bam, referenceFile, output_directory)
print(args)

int_dict = Dict.empty(
    key_type=types.unicode_type,
    value_type=types.int64
)
for key in range(50000):
    int_dict[str(key)] = key

def compress_state_vector(state_vectors):
    state_strings = dict()
    for k, v in tqdm(state_vectors.items()):
        state_strings[k] = _compress_state_vector_to_state_string(v)
    return state_strings

@njit
def _compress_state_vector_to_state_string(state_vector):
    compressed_string = ""
    state = state_vector[0]
    if state == 0:
        prev_char = "I"
    elif state == 1:
        prev_char = "A"
    elif np.isnan(state):
        prev_char = "N"
    else:
        prev_char = "C"

    count = 1
    current_char = ""
    for state in state_vector[1:]:
        if state == 0:
            current_char = "I"
        elif state == 1:
            current_char = "A"
        elif np.isnan(state):
            current_char = "N"
        else:
            current_char = "C"

        if current_char != prev_char:
            compressed_string += str(count) + prev_char
            count = 1
            prev_char = current_char
        else:
            count += 1

    if current_char == prev_char:
        compressed_string += str(count) + prev_char

    return compressed_string

@njit
def convert_state_string_to_MM_ML(state_string, int_dict):
    cur_count = ""
    ops = ["N", "C", "I", "A"]
    append_dict = {"I": 0, "A": 1}
    MM = []
    ML = []
    previous_skip = 0
    prev_op = ""

    for i in range(len(state_string)):
        char = state_string[i]
        if char in ops:
            cur_count_int = int_dict[cur_count]
            if char == "N":
                previous_skip += cur_count_int
            else:
                MM.append(previous_skip)
                if char == "C":                
                    append_prob = append_dict[prev_op]
                    numerator = 1 - 2 * append_prob
                    slope = numerator / (cur_count_int + 1)
                    ML.append(int(255 * (slope + append_prob)))
                    for j in range(1, cur_count_int):
                        MM.append(0)
                        ML.append(int(255 * (slope * (j + 1) + append_prob)))
                else:
                    append_prob = int(255 * (append_dict[char]))
                    ML.append(append_prob)
                    for j in range(1, cur_count_int):
                        MM.append(0)
                        ML.append(append_prob)
                previous_skip = 0
            prev_op = char
            cur_count = ""
        else:
            cur_count += char

    return MM, ML

def main():
    sampleRef = pd.read_csv(referenceFile, sep=',', index_col='index')
    print(sampleRef)

    hmmDict = {}
    for samp in samples:
        with open('{}processed/binarized/HMMv2out/{}_{}_HMMres.pickle'.format(output_directory, sampleRef['cell'][samp], sampleRef['sampleName'][samp]), 'rb') as f:
            hmmDict.update(pickle.load(f))
    
    state_strings = compress_state_vector(hmmDict)
    
    infile = pysam.AlignmentFile(merged_bam, check_sq=False)
    outfile = pysam.AlignmentFile('{}processed/annot/{}'.format(output_directory, os.path.basename(merged_bam).replace('.bam', '.samosa.bam')), 'wb', check_sq=False, template=infile)

    for read in tqdm(infile):
        zmw = read.get_tag('zmw')
        if zmw in state_strings:
            mm, ml = convert_state_string_to_MM_ML(state_strings[zmw], int_dict)
            read.set_tags(
                read.get_tags() +
                [
                    ("Mm", "N+a,{};".format(",".join(map(str, mm))), "Z"),
                    ("Ml", ml, "C")
                ]
            )
            outfile.write(read)

    infile.close()
    outfile.close()
    print("complete!")

if __name__ == '__main__':
    main()
