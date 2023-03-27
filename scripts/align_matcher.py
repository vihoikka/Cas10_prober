import pandas as pd
import os
import argparse
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import ast
import Bio

parser = argparse.ArgumentParser()
parser.add_argument("--cas10", required = True)
parser.add_argument("--cas7", required = True)
parser.add_argument("--cas5", required = True)
parser.add_argument("--cora", required = True)
parser.add_argument("--rna16s", required = True)
parser.add_argument("--outfolder", required = True)

args = parser.parse_args()
cas10 = args.cas10
cas7 = args.cas7
cas5 = args.cas5
cora = args.cora
rna16s= args.rna16s
outfolder = args.outfolder


import itertools

# Define a function to read a FASTA-format alignment file
def read_alignment(filename):
    sequences = {}
    with open(filename) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip()
                sequences[name] = ""
            else:
                sequences[name] += line.strip()
    return sequences

# Define the list of alignment filenames
filenames = [cas10, cas7, cas5, cora, rna16s]

# Read each alignment file and store the sequences in a dictionary
alignments = {}
for filename in filenames:
    sequences = read_alignment(filename)
    alignments[filename] = sequences

# Find the pairwise combinations of alignments
pairs = list(itertools.combinations(filenames, 2))

# Loop over each pair of alignments
for pair in pairs:
    # Get the sequences from each alignment
    seq1 = alignments[pair[0]]
    seq2 = alignments[pair[1]]

    new_dict_seq1 = {}
    new_dict_seq2 = {}

    for key in seq1.keys():
        new_key = key.split()[0]
        new_dict_seq1[new_key] = seq1[key]

    for key in seq2.keys():
        new_key = key.split()[0]
        new_dict_seq2[new_key] = seq2[key]

    seq1 = new_dict_seq1
    seq2 = new_dict_seq2

    # Get the set intersection of the entry names
    names1 = set(seq1.keys())
    names2 = set(seq2.keys())
    shared_names = names1.intersection(names2)

    # Extract only the filename (without the directory path and the "alignment" part)
    pair_filenames = [os.path.basename(f)[:-4].split("_")[0] for f in pair]

    base_folder_name = f"{outfolder}/{pair_filenames[0]}_{pair_filenames[1]}"
    os.makedirs(base_folder_name)

    # Write new alignment files with only the shared entries
    for name in shared_names:
        with open(f"{base_folder_name}/{pair_filenames[0]}_{pair_filenames[1]}_{pair_filenames[0]}.afa", "a") as f1, \
                open(f"{base_folder_name}/{pair_filenames[0]}_{pair_filenames[1]}_{pair_filenames[1]}.afa", "a") as f2:
            f1.write(">" + name + "\n" + seq1[name] + "\n")
            f2.write(">" + name + "\n" + seq2[name] + "\n")

    with open(f"{base_folder_name}/{pair_filenames[0]}_{pair_filenames[1]}_common.csv", "a") as commonNamesFile:
        commonNamesFile.write("Genome\t" + pair_filenames[0] + "_" + pair_filenames[1] + "\n")
        for name in shared_names:
            commonNamesFile.write(name + ",True\n")

    with open(f"{outfolder}/{pair_filenames[0]}_{pair_filenames[1]}.done", "w") as donefile:
        donefile.write("done")