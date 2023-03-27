import pandas as pd
import os
import argparse
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor


parser = argparse.ArgumentParser()
parser.add_argument("-ms", "--min_spacers", required = True)
parser.add_argument("-ta", "--targets", required = True)
parser.add_argument("-pcs", "--nodes", required = True)
parser.add_argument("-oi", "--out_interesting", required = True)
parser.add_argument("-bp", "--base_path", required=False)
parser.add_argument("-s", "--species", required=False)
parser.add_argument("-th", "--threshold", required=True)
parser.add_argument("-gg", "--gene_to_genome", required=True)
parser.add_argument("-pgb", "--phages_genbank", required=True)

args = parser.parse_args()
targets = args.targets
threshold = float(args.threshold)
min_spacers = int(args.min_spacers)
nodes = args.nodes
out_interesting = args.out_interesting
phages_genbank = args.phages_genbank
#base_path = args.base_path
#species = args.species
gene_to_genome = args.gene_to_genome

#Load phage info into dataframe
phages = pd.read_csv(nodes, sep = '\t', header=0)

final_output = pd.DataFrame({'temp' : []})

temp_all_nodes = phages

targets = targets.split(",")
results = {}

#Subset into interesting phages
counter = 0
for t in targets:
    interesting_phages = phages[(phages["fraction_" + t] >= threshold) & (phages["interactions"] >= min_spacers)]
    #print(interesting_phages)
    if counter == 0:
        final_output = pd.merge(temp_all_nodes, interesting_phages, on="Id", how="right", suffixes=('', '_y'))
        final_output.drop(final_output.filter(regex='_y$').columns, axis=1, inplace=True)
    else:
        final_output = pd.merge(final_output, interesting_phages, on="Id", how="outer", suffixes=('', '_y'))
        final_output.drop(final_output.filter(regex='_y$').columns, axis=1, inplace=True)
    counter += 1

final_output.dropna(axis=1, how='all', inplace=True)

print(final_output)
#make list out of interesting phage accessions

final_output.to_csv(out_interesting, index = False, sep = "\t", header = False)


#
#
# with open(out_interesting, mode='wt', encoding='utf-8') as file:
#     file.write('\n'.join(interesting_proteins))
#
# with open(out_all, mode='wt', encoding='utf-8') as file:
#     file.write('\n'.join(all_proteins))
#
#  #Load gene to genome map file
# gg = pd.read_csv(gene_to_genome, sep=",", header=0)
#
# #Link genes to genomes
# interesting_proteins = pd.merge(gg, interesting_phages, left_on="contig_id", right_on = "phage", how="right")
# all_proteins = pd.merge(gg, phages, left_on="contig_id", right_on = "phage", how="right")
# print(str(interesting_proteins.shape))
# print(str(all_proteins.shape))
#
# interesting_proteins.to_csv(out_interesting, index=False, sep = '\t', header = True)
# all_proteins.to_csv(out_all, index=False, sep = '\t', header = True)
