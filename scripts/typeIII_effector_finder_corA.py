import argparse
import os
import pandas as pd
import gffutils
import ast
from Bio import SeqIO
import regex

parser = argparse.ArgumentParser(description='Process some inputs.')
parser.add_argument('--sample', type=str, help='Input sample', required=True)
#parser.add_argument('--cctyper', type=str, help='Input CCTyper', required=True)
parser.add_argument('--CorAFinder', type=str, help='Input CorAFinder', required=True)
parser.add_argument('--gff', type=str, help='Input GFF', required=True)
parser.add_argument('--out', type=str, help='Output file', required=True)
parser.add_argument('--cctyper_folder', type=str, help='Cctyper folder', required=True)
parser.add_argument('--this_folder', type=str, help='This folder', required=True)
parser.add_argument('--proteins', type=str, help='This folder', required=True)

finalTable = pd.DataFrame(columns=["protein", "cc_typer_evalue", "cctyper_annotation", "genome_annotation", "seq", "operon", "operon_pos", "host","cas10_cyclase","Cas10_cyclase_match_seq"])

cctyper_evalue_threshold = 1e-10

Cas10_evalue_threshold = 1e-8
cas10_search_range = 300 #when matching cctyper cas10 location with the gff file, expand search range up- and downstream by this amount (bp)
general_gene_coordinate_range = 200 #same as above, but for other genes than cas10

args = parser.parse_args()

sample = args.sample
#cctyper = args.cctyper
CorAFinder = args.CorAFinder
gff = args.gff
out = args.out
cctyper_folder = args.cctyper_folder
this_folder = args.this_folder
proteins = args.proteins

cas_operons = cctyper_folder + "/cas_operons.tab"
crispr_cas = cctyper_folder + "/CRISPR_Cas.tab"
cctyper_gene_positions = cctyper_folder + "/genes.tab"

def getGeneLocationCRISPRCas(cctyper_gene_positions, contig, position_order):
    #print(cctyper_gene_positions)
    #print(contig)
    #print(position_order)
    positions = cctyper_gene_positions #df of where gene position translates to coordinates
    genepos = positions.loc[(positions["Contig"] == contig) & (positions["Pos"] == position_order)] #find the row that corresponds to current operon
    #print(genepos)
    return [genepos["Start"].values[0],genepos["End"].values[0]]

def getCas10CyclaseDomain(cas10):
    seq = str(cas10.seq)
    #print(seq)
    #regex = "r'G.DD'"
    #motif = re.search(r"\G.DD", seq)
    motif_fuzzy = regex.findall("(GGDD){e<=1}", seq, overlapped=True)
    return(motif_fuzzy)


#print(cas_operons)
#test if cctyper results exist
if not os.path.isfile(crispr_cas):
    print("No CRISPR-Cas loci associated with this sample")
    exit()
else: #store type III loci into a dictionary
    cas_operons = pd.read_csv(cas_operons, sep="\t", header = 0, converters={"Genes":pd.eval, "E-values":pd.eval, "Positions":pd.eval})
    crispr_cas = pd.read_csv(crispr_cas, sep="\t", header = 0)
    cctyper_gene_positions = pd.read_csv(cctyper_gene_positions, sep="\t", header = 0)

    type_III_loci = {}
    for index, row in crispr_cas.iterrows():
        if "III-" in row["Prediction"]:
            operon_information = {}
            operon_information["contig"] = row["Contig"]
            operon_information["pos"] = row["Operon_Pos"]
            operon_information["start"] = row["Operon_Pos"].strip("[]").split(",")[0]
            operon_information["end"] = row["Operon_Pos"].strip("[]").split(",")[1]
            operon_information["type"] = row["Prediction"]
            operon_information["coa_genes"] = []
            operon_information["coa_type"] = []
            operon_information["coa_acc"] = []
            operon_information["has_effector"] = False
            operon_information["Cas10_cyclase"] = False
            operon_information["Cas10_cyclase_match_seq"] = []
            type_III_loci[row["Operon"]] = operon_information
            #print("Found CRISPR-Cas type " + str(type_III_loci[row["Operon"]]))

#don't continue if type III CRISPR-Cas loci not found
if len(type_III_loci) == 0:
    print("No type III CRISPR-Cas loci found. Exiting type III effector finder")
    exit()

#Load CorAFinder data
CorAFinder = pd.read_csv(CorAFinder, delim_whitespace=True, header = 0, index_col=False)

#read gff using gffutils
gff_db = gffutils.create_db(gff, dbfn=this_folder + '_gff.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
gff = gffutils.FeatureDB(this_folder + '_gff.db', keep_order=True)

#load proteins with biopython
proteins = SeqIO.to_dict(SeqIO.parse(proteins, "fasta"))


#lookup CorAFinder results. Are any effectors located in a type III locus?
#continue_script = False #used as a check later
for key, crispr_locus in type_III_loci.items(): #loop through CRISPR-Cas III loci
    for index, row in CorAFinder.iterrows():
        protein = row["query_name"]
        gff_protein = gff["cds-" + protein]
        #print(row["query_name"])
        contig = gff_protein.seqid

        #if cOA-gene within CRISPR-Cas locus and on same contig
        if (int(gff_protein.start) >= int(crispr_locus["start"])) & (int(gff_protein.start) <= int(crispr_locus["end"])) & (contig == crispr_locus["contig"]):
            type_III_loci[key]["coa_genes"].append(row["target_name"])
            type_III_loci[key]["coa_type"].append(row["target_name"].split("_")[0])
            type_III_loci[key]["coa_acc"].append(protein)
            type_III_loci[key]["has_CorA"] = True
            #continue_script = True
        #is_on_crispr(row[""])

#if continue_script == False:#
#    exit()

print(type_III_loci)


#Go through all genes in the cctyper output
for key, value in type_III_loci.items():
    print("Current CRISPR-Cas locus is " + str(value["type"]))
    crispr_locus_table = finalTable
    genes = {}
    row = cas_operons.loc[cas_operons["Operon"] == key] #find the row that corresponds to current operon
    genelist = row["Genes"].tolist()[0]
    gene_pos_list = row["Positions"].tolist()[0]
    evaluelist = row["E-values"].tolist()[0]

    for id, gene in enumerate(genelist):
        #print(gene)
        if ("Cas10" in gene) & (float(evaluelist[id]) <= Cas10_evalue_threshold): #get information on Cas10 cyclase motif, if we trust that this is Cas10
            cas10_cctyper_pos = gene_pos_list[id]
            cas10_global_start_end_coordinates = getGeneLocationCRISPRCas(cctyper_gene_positions, value["contig"], cas10_cctyper_pos) #where in the genome is Cas10, according to cctyper
            cas10_start = int(cas10_global_start_end_coordinates[0]) - cas10_search_range
            cas10_end = int(cas10_global_start_end_coordinates[1]) + cas10_search_range
            print("The Cas10 start: " + str(cas10_start))
            print("The Cas10 end: " + str(cas10_end))
            gff_cas10 = gff.region(seqid=value["contig"], start = cas10_start, end = cas10_end, completely_within = True, featuretype = "CDS") #find corresponding coordinates from gff file
            print("The Cas10 coordinates: " + str(gff_cas10))
            gff_cas10 = list(gff_cas10)[0] #convert generator to list
            print("gff cas 10 type: " + str(type(gff_cas10)))
            print("Opened generator object: " + str(gff_cas10))
            try:
                cas10_genomic_id = gff_cas10.attributes["Name"][0] #pick the code for the Cas10 protein from the gff object
                cas10_seq = proteins[cas10_genomic_id] #find this protein in the original .faa file
                cas10_cyclase_domains = getCas10CyclaseDomain(cas10_seq)

                if len(cas10_cyclase_domains) > 0: #if GGDD motif (with one mismatch allowed) is found, mark it as so
                    print("GGDD found")
                    type_III_loci[key]["Cas10_cyclase"] = True
                    type_III_loci[key]["Cas10_cyclase_match_seq"] = cas10_cyclase_domains
            except:
                print("Could extract Cas10 (" + str(gff_cas10) + ")")

        gene_info = {}
        gene_info["evalue"] = float(evaluelist[id])
        gene_info["flagged"] = False
        print("Gene E-value: " + str(gene_info["evalue"]) + " (" + str(gene) + ")")

        if gene_info["evalue"] > cctyper_evalue_threshold:
            print("E-value crossed")
            gene_info["flagged"] = True

            global_start_end_coordinates = getGeneLocationCRISPRCas(cctyper_gene_positions, value["contig"], gene_pos_list[id]) #returns cctyper gene coordinates (start and stop as list, e.g. [4,235])
            gff_protein_object = gff.region(seqid=value["contig"], start = int(global_start_end_coordinates[0])-general_gene_coordinate_range, end = int(global_start_end_coordinates[1])+general_gene_coordinate_range, completely_within = True, featuretype = "CDS") #find corresponding coordinates from gff file
            print(gff_protein_object)
            gff_protein = list(gff_protein_object)[0] #convert generator to list
            print("-----")
            print(gff_protein.attributes)
            print("-----")
            protein_id = gff_protein.attributes["Name"][0]
            protein_annotation = gff_protein.attributes["product"][0]
            protein_seq = str(proteins[protein_id].seq)


            df = pd.DataFrame({
                "protein": [protein_id],
                "cc_typer_evalue": [gene_info["evalue"]],
                "cctyper_annotation": [genelist[id]],
                "genome_annotation": [protein_annotation],
                "seq": [protein_seq],
                "locus_pos": [gene_pos_list[id]],
                "host": [sample],
                "cas10_cyclase": ["-"],
                "cas10_cyclase_match_seq": ["-"]
            })


            crispr_locus_table = pd.concat([crispr_locus_table,df])

    crispr_locus_table = crispr_locus_table.assign(cas10_cyclase = type_III_loci[key]["Cas10_cyclase"])
    crispr_locus_table = crispr_locus_table.assign(cas10_cyclase_match_seq = str(type_III_loci[key]["Cas10_cyclase_match_seq"]))

    finalTable = pd.concat([finalTable, crispr_locus_table])

print("Finished. Final table:")
print(finalTable)

finalTable.to_csv(out, index=True, sep = '\t', header = True)