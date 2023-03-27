import pandas as pd
import os
import argparse
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor


parser = argparse.ArgumentParser()
parser.add_argument("-ms", "--min_spacers", required = True)
parser.add_argument("-ta", "--target", required = True)
parser.add_argument("-pcs", "--nodes", required = True)
parser.add_argument("-oi", "--out_interesting", required = True)
parser.add_argument("-oa", "--out_all", required = True)
parser.add_argument("-bp", "--base_path", required=False)
parser.add_argument("-s", "--species", required=False)
parser.add_argument("-th", "--threshold", required=True)
parser.add_argument("-gg", "--gene_to_genome", required=True)
parser.add_argument("-pgb", "--phages_genbank", required=True)

args = parser.parse_args()
target = args.target
threshold = float(args.threshold)
min_spacers = int(args.min_spacers)
nodes = args.nodes
out_interesting = args.out_interesting
out_all = args.out_all
phages_genbank = args.phages_genbank
#base_path = args.base_path
#species = args.species
gene_to_genome = args.gene_to_genome

#Load phage info into dataframe
phages = pd.read_csv(nodes, sep = '\t', header=0)
print(phages)
#Subset into interesting phages
interesting_phages = phages[(phages["fraction_" + target] >= threshold) & (phages["interactions"] >= min_spacers)]
#make list out of interesting phage accessions
interesting_phages_list = interesting_phages.phage.values.tolist()
print(interesting_phages_list)


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

all_new_records = []
interesting_proteins = []
all_proteins = []

#NOTE: this version only saves proteins from interesting phages. All other
#phage proteins are discarded since scanning through all of them takes too
#much time. TODO: Make a separate rule to gather all phage proteins.

all_phages = SeqIO.parse(open(phages_genbank,"r"), "genbank")
print(all_proteins)

def checkPhage(phage):
    phagename = phage.name
    if phagename in interesting_phages_list:
        for feature in tqdm(phage.features):
            seq = feature.qualifiers["translation"][0]
            id = str(feature.qualifiers["protein_id"][0])
            new_record = SeqRecord(
                Seq(str(seq)),
                id=id,
                description=str(feature.type + "|" + feature.qualifiers.keys()) #this is used later to check if its and actual protein
                )
            return new_record

with ProcessPoolExecutor(max_workers=32) as executor:
    for r in tqdm(executor.map(checkPhage, all_phages)):
        interesting_proteins.append(r)


# for phage in tqdm(all_phages):
#     if phage.name in interesting_phages_list:
#         for feature in phage.features:
#             seq = feature.qualifiers["translation"][0]
#             id = str(feature.qualifiers["protein_id"][0])
#             new_record = SeqRecord(
#                 Seq(str(seq)),
#                 id=id,
#                 description=str(feature.type + "|" + feature.qualifiers.keys()) #this is used later to check if its and actual protein
#                 )
#             interesting_proteins.append(new_record)
        #else:
            #all_proteins.append(new_record)

#we need to filter out pseudogenes or other weird entries with no protein translation
filtered_interesting_proteins = []
for protein in interesting_proteins:
    type = protein.description.split("|")[0]
    qualifiers = protein.description.split("|")[1]
    if (type == "CDS") and ("pseudo" not in qualifiers.keys()) and ("translation" in feature.qualifiers.keys()):
        filtered_interesting_proteins.append(protein)




# for record in tqdm(SeqIO.parse(open(phages_genbank,"r"), "genbank")):
#     for feature in record.features:
#         #we need to filter out pseudogenes or other weird entries with no protein translation
#         if (feature.type == "CDS") and ("pseudo" not in feature.qualifiers.keys()) and ("translation" in feature.qualifiers.keys()):
#                 seq = feature.qualifiers["translation"][0]
#                 id = str(feature.qualifiers["protein_id"][0])
#                 new_record = SeqRecord(
#                     Seq(str(seq)),
#                     id=id,
#                     description=""
#                     )
#                 if record.name in interesting_phages_list:
#                     interesting_proteins.append(new_record)
#                 else:
#                     all_proteins.append(new_record)


print("Interesting proteins: " + str(len(interesting_proteins)))
print("Interesting proteins: " + str(len(filtered_interesting_proteins)))

#print("All proteins: " + str(len(all_proteins)))

#with open(out_all, "w") as output_handle:
#    SeqIO.write(all_proteins, output_handle, "fasta")

with open(out_interesting, "w") as output_handle:
    SeqIO.write(filtered_interesting_proteins, output_handle, "fasta")
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
