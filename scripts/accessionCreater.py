import pandas as pd
import os
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("-ms", "--min_spacers", required = True)
parser.add_argument("-pcs", "--phage_csm6", required = True)
parser.add_argument("-oi", "--out_interesting", required = True)
parser.add_argument("-oa", "--out_all", required = True)
parser.add_argument("-bp", "--base_path", required=False)
parser.add_argument("-s", "--species", required=False)
parser.add_argument("-ct", "--csm6_threshold", required=True) 
parser.add_argument("-gg", "--gene_to_genome", required=True)
parser.add_argument("-pgb", "--phages_genbank", required=True)

args = parser.parse_args()
csm6_threshold = float(args.csm6_threshold)
min_spacers = int(args.min_spacers)
phage_csm6 = args.phage_csm6
out_interesting = args.out_interesting
out_all = args.out_all
phages_genbank = args.phages_genbank
#base_path = args.base_path
#species = args.species
gene_to_genome = args.gene_to_genome

#Load phage csm6 info into dataframe
phages = pd.read_csv(phage_csm6, sep = '\t', header=0)
#Subset into interesting phages
interesting_phages = phages[(phages['fraction_csm6'] >= csm6_threshold) & (phages["interactions"] >= min_spacers)]
#make list out of interesting phage accessions
interesting_phages_list = interesting_phages.phage.values.tolist()
print(interesting_phages_list)


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

all_new_records = []
interesting_proteins = []
all_proteins = []
for record in tqdm(SeqIO.parse(open(phages_genbank,"r"), "genbank")):
    for feature in record.features:
        #we need to filter out pseudogenes or other weird entries with no protein translation
        if (feature.type == "CDS") and ("pseudo" not in feature.qualifiers.keys()) and ("translation" in feature.qualifiers.keys()):
                seq = feature.qualifiers["translation"][0]
                id = str(feature.qualifiers["protein_id"][0])
                new_record = SeqRecord(
                    Seq(str(seq)),
                    id=id,
                    description=""
                    )
                if record.name in interesting_phages_list:
                    interesting_proteins.append(new_record)
                else:
                    all_proteins.append(new_record)

                            
print("Interesting proteins: " + str(len(interesting_proteins)))
print("All proteins: " + str(len(all_proteins)))

with open(out_all, "w") as output_handle:
    SeqIO.write(all_proteins, output_handle, "fasta")

with open(out_interesting, "w") as output_handle:
    SeqIO.write(interesting_proteins, output_handle, "fasta")
     
"""
with open(out_interesting, mode='wt', encoding='utf-8') as file:
    file.write('\n'.join(interesting_proteins))

with open(out_all, mode='wt', encoding='utf-8') as file:
    file.write('\n'.join(all_proteins))

 #Load gene to genome map file
gg = pd.read_csv(gene_to_genome, sep=",", header=0)

#Link genes to genomes
interesting_proteins = pd.merge(gg, interesting_phages, left_on="contig_id", right_on = "phage", how="right")
all_proteins = pd.merge(gg, phages, left_on="contig_id", right_on = "phage", how="right")
print(str(interesting_proteins.shape))
print(str(all_proteins.shape))

interesting_proteins.to_csv(out_interesting, index=False, sep = '\t', header = True)
all_proteins.to_csv(out_all, index=False, sep = '\t', header = True)
 """