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
parser.add_argument('--out_cora_info', type=str, help='Output file', required=True)
parser.add_argument('--out_cora', type=str, help='Output file', required=True)
parser.add_argument('--out_cas10', type=str, help='Output file', required=True)
parser.add_argument('--out_cas5', type=str, help='Output file', required=True)
parser.add_argument('--out_cas7', type=str, help='Output file', required=True)
parser.add_argument('--out_cora_cctyper', type=str, help='Output file', required=True)
parser.add_argument('--out_info', type=str, help='Output file', required=True)
parser.add_argument('--cctyper_folder', type=str, help='Cctyper folder', required=True)
parser.add_argument('--this_folder', type=str, help='This folder', required=True)
parser.add_argument('--proteins', type=str, help='This folder', required=True)


cctyper_evalue_threshold = 1e-10

Cas10_evalue_threshold = 1e-20
cas10_search_range = 300 #when matching cctyper cas10 location with the gff file, expand search range up- and downstream by this amount (bp)
general_gene_coordinate_range = 1000 #same as above, but for other genes than cas10
corA_max_distance_from_crispr = 2000
cas10_length_cutoff = 500

args = parser.parse_args()

sample = args.sample
#cctyper = args.cctyper
CorAFinder = args.CorAFinder
out_cora_info = args.out_cora_info
out_cora = args.out_cora
out_cora_cctyper = args.out_cora_cctyper
out_cas10 = args.out_cas10
out_cas5 = args.out_cas5
out_cas7 = args.out_cas7
out_info = args.out_info
gff = args.gff
cctyper_folder = args.cctyper_folder
this_folder = args.this_folder
proteins = args.proteins

cas_operons = cctyper_folder + "/cas_operons.tab"
crispr_cas = cctyper_folder + "/CRISPR_Cas.tab"
cctyper_gene_positions = cctyper_folder + "/genes.tab"

finalTable = pd.DataFrame(columns=["type_iii", "CorA_in_iii_locus", "host", "CorA_seq", "Cas10_seq"])

info_table = pd.DataFrame({"Sample" : sample,
                    "CorA" : False,
                    "Cas10" : False,
                    "GGDD" : False,
                    "GGDD_seq" : "",
                    "Subtype" : ""
                    },index=[0])

def getGeneLocationCRISPRCas(cctyper_gene_positions, contig, position_order):
    #print(cctyper_gene_positions)
    #print(contig)
    #print(position_order)
    positions = cctyper_gene_positions #df of where gene position translates to coordinates
    genepos = positions.loc[(positions["Contig"] == contig) & (positions["Pos"] == position_order)] #find the row that corresponds to current operon
    #print(genepos)
    return [genepos["Start"].values[0],genepos["End"].values[0]]

def getCas10CyclaseDomain(cas10):
    seq = str(cas10)
    #regex = "r'G.DD'"
    #motif = re.search(r"\G.DD", seq)
    motif_fuzzy = regex.findall("(GGDD){e<=1}", seq, overlapped=True)
    if "GGDD" in motif_fuzzy:
        return "GGDD"
    else:
        return str(motif_fuzzy[0]).strip("][")


#print(cas_operons)
#test if cctyper results exist
if (not os.path.isfile(crispr_cas) or (not os.path.isfile(cas_operons))):
    print("No CRISPR-Cas loci associated with this sample")
    finalTable.to_csv(out_cora_info, sep = '\t', header = True, index = False)
    info_table.to_csv(out_info, sep = '\t', header = False, index = False)
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
            operon_information["CorA_genes"] = []
            operon_information["CorA_type"] = []
            operon_information["CorA_acc"] = []
            operon_information["has_CorA"] = False
            operon_information["Cas10_cyclase"] = False
            operon_information["Cas10_cyclase_match_seq"] = []
            operon_information["CorA_seq"] = []
            operon_information["CorA_object"] = []
            operon_information["CorA_seq_cctyper"] = []
            operon_information["CorA_object_cctyper"] = []
            operon_information["Cas10_boolean"] = False
            operon_information["Cas10_object"] = []
            operon_information["Cas10_seq"] = []
            operon_information["Cas5_object"] = []
            operon_information["Cas5_seq"] = []
            operon_information["Cas7_object"] = []
            operon_information["Cas7_seq"] = []
            type_III_loci[row["Operon"]] = operon_information
            #print("Found CRISPR-Cas type " + str(type_III_loci[row["Operon"]]))

#don't continue if type III CRISPR-Cas loci not found
if len(type_III_loci) == 0:
    print("No type III CRISPR-Cas loci found. Exiting type III effector finder")
    info_table.to_csv(out_info, sep = '\t', header = False, index = False)
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
        gffkey = "cds-" + protein
        #if str(gffkey) in gff: #usually proteins are marked with a cds-prefix in the gff file
        cdsNotFound = False
        coraNotFound = False
        try:
            gff_protein = gff["cds-" + protein]
            contig = gff_protein.seqid
        except: 
            cdsNotFound = True
            pass
        if cdsNotFound == True:
            try:
                gff_protein = gff[protein]
                contig = gff_protein.seqid
            except:
                coraNotFound = True
                pass
        if coraNotFound == True:
            print("This CorA protein was not found in the gff file: " + str(protein))
            gff_protein = "Not found"
        
        #if CorA within CRISPR-Cas locus and on same contig
        if (gff_protein != "Not found") and (int(gff_protein.start) >= int(crispr_locus["start"])) & (int(gff_protein.start) <= int(crispr_locus["end"])+int(corA_max_distance_from_crispr)) & (contig == crispr_locus["contig"]):
            CorA_genomic_id = gff_protein.attributes["Name"][0]
            CorA_object = proteins[CorA_genomic_id]
            CorA_object.id = sample
            CorA_object.description = ""
            CorA_seq = CorA_object.seq
            type_III_loci[key]["CorA_genes"].append(row["target_name"])
            type_III_loci[key]["CorA_type"].append(row["target_name"].split("_")[0])
            type_III_loci[key]["CorA_acc"].append(protein)
            type_III_loci[key]["has_CorA"] = True
            type_III_loci[key]["CorA_seq"] = CorA_seq
            type_III_loci[key]["CorA_object"] = CorA_object

        #if CorA is not found, the dictionary remains unchanged (the default values are for CorA not being found)

#see naming scheme at https://github.com/Russel88/CRISPRCasTyper/blob/99719a135c36c03ce6d2b84a8e4b28022e6104d6/data/interference.json
cas5_names = ["Cas5", "Cmr3", "Csm4"]
cas7_names = ["Cas7", "Csm3", "Cmr4"]

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

        if any(name in gene for name in cas5_names) & (float(evaluelist[id]) <= Cas10_evalue_threshold): #get information on Cas10 cyclase motif, if we trust that this is Cas10
            cas5_cctyper_pos = gene_pos_list[id]
            cas5_global_start_end_coordinates = getGeneLocationCRISPRCas(cctyper_gene_positions, value["contig"], cas5_cctyper_pos) #where in the genome is cas5, according to cctyper
            cas5_start = int(cas5_global_start_end_coordinates[0]) - cas10_search_range
            cas5_end = int(cas5_global_start_end_coordinates[1]) + cas10_search_range

            try:
                gff_cas5 = gff.region(seqid=value["contig"], start = cas5_start, end = cas5_end, completely_within = True, featuretype = "CDS") #find corresponding coordinates from gff file
                print("The cas5 coordinates: " + str(gff_cas5))
                gff_cas5 = list(gff_cas5)[0] #convert generator to list
                print("gff cas5 type: " + str(type(gff_cas5)))
                print("Opened generator object: " + str(gff_cas5))

                gff_cas5 = gff.region(seqid=value["contig"], start = cas5_start, end = cas5_end, completely_within = True, featuretype = "CDS") #find corresponding coordinates from gff file
                gff_cas5 = list(gff_cas5)[0] #convert generator to list

            except:
                print("Could not pinpoint cas5 genomic position")
                continue
            try:
                cas5_genomic_id = gff_cas5.attributes["Name"][0] #pick the code for the cas5 protein from the gff object
                cas5_object = proteins[cas5_genomic_id]
                cas5_object.id = sample
                cas5_object.description = ""
                type_III_loci[key]["Cas5_object"] = cas5_object
                type_III_loci[key]["Cas5_seq"] = cas5_object.seq #find this protein in the original .faa file

            except:
                print("Could not extract Cas5 (" + str(gff_cas5) + ")")
                continue

        if any(name in gene for name in cas7_names) & (float(evaluelist[id]) <= Cas10_evalue_threshold): #get information on Cas10 cyclase motif, if we trust that this is Cas10
            cas7_cctyper_pos = gene_pos_list[id]
            cas7_global_start_end_coordinates = getGeneLocationCRISPRCas(cctyper_gene_positions, value["contig"], cas7_cctyper_pos) #where in the genome is cas7, according to cctyper
            cas7_start = int(cas7_global_start_end_coordinates[0]) - cas10_search_range
            cas7_end = int(cas7_global_start_end_coordinates[1]) + cas10_search_range

            try:
                gff_cas7 = gff.region(seqid=value["contig"], start = cas7_start, end = cas7_end, completely_within = True, featuretype = "CDS") #find corresponding coordinates from gff file
                print("The cas7 coordinates: " + str(gff_cas7))
                gff_cas7 = list(gff_cas7)[0] #convert generator to list
                print("gff cas7 type: " + str(type(gff_cas7)))
                print("Opened generator object: " + str(gff_cas7))

                gff_cas7 = gff.region(seqid=value["contig"], start = cas7_start, end = cas7_end, completely_within = True, featuretype = "CDS") #find corresponding coordinates from gff file
                gff_cas7 = list(gff_cas7)[0] #convert generator to list

            except:
                print("Could not pinpoint cas7 genomic position")
                continue
            try:
                cas7_genomic_id = gff_cas7.attributes["Name"][0] #pick the code for the cas7 protein from the gff object
                cas7_object = proteins[cas7_genomic_id]
                cas7_object.id = sample
                cas7_object.description = ""
                type_III_loci[key]["Cas7_object"] = cas7_object
                type_III_loci[key]["Cas7_seq"] = cas7_object.seq #find this protein in the original .faa file

            except:
                print("Could not extract Cas7 (" + str(gff_cas7) + ")")
                continue

        if ("CorA".lower() in gene.lower()) & (float(evaluelist[id]) <= Cas10_evalue_threshold): #get information on Cas10 cyclase motif, if we trust that this is Cas10
            cora_cctyper_cctyper_pos = gene_pos_list[id]
            cora_cctyper_global_start_end_coordinates = getGeneLocationCRISPRCas(cctyper_gene_positions, value["contig"], cora_cctyper_cctyper_pos) #where in the genome is cora_cctyper, according to cctyper
            cora_cctyper_start = int(cora_cctyper_global_start_end_coordinates[0]) - cas10_search_range
            cora_cctyper_end = int(cora_cctyper_global_start_end_coordinates[1]) + cas10_search_range

            try:
                gff_cora_cctyper = gff.region(seqid=value["contig"], start = cora_cctyper_start, end = cora_cctyper_end, completely_within = True, featuretype = "CDS") #find corresponding coordinates from gff file
                print("The cora_cctyper coordinates: " + str(gff_cora_cctyper))
                gff_cora_cctyper = list(gff_cora_cctyper)[0] #convert generator to list
                print("gff cora_cctyper type: " + str(type(gff_cora_cctyper)))
                print("Opened generator object: " + str(gff_cora_cctyper))

                gff_cora_cctyper = gff.region(seqid=value["contig"], start = cora_cctyper_start, end = cora_cctyper_end, completely_within = True, featuretype = "CDS") #find corresponding coordinates from gff file
                gff_cora_cctyper = list(gff_cora_cctyper)[0] #convert generator to list

            except:
                print("Could not pinpoint cora_cctyper genomic position")
                continue
            try:
                cora_cctyper_genomic_id = gff_cora_cctyper.attributes["Name"][0] #pick the code for the cora_cctyper protein from the gff object
                cora_cctyper_object = proteins[cora_cctyper_genomic_id]
                cora_cctyper_object.id = sample
                cora_cctyper_object.description = ""
                type_III_loci[key]["CorA_cctyper_object"] = cora_cctyper_object
                type_III_loci[key]["CorA_cctyper_seq"] = cora_cctyper_object.seq #find this protein in the original .faa file

            except:
                print("Could not extract Cora_cctyper (" + str(gff_cora_cctyper) + ")")
                continue

        if ("Cas10" in gene) & (float(evaluelist[id]) <= Cas10_evalue_threshold): #get information on Cas10 cyclase motif, if we trust that this is Cas10
            cas10_cctyper_pos = gene_pos_list[id]
            cas10_global_start_end_coordinates = getGeneLocationCRISPRCas(cctyper_gene_positions, value["contig"], cas10_cctyper_pos) #where in the genome is Cas10, according to cctyper
            cas10_start = int(cas10_global_start_end_coordinates[0]) - cas10_search_range
            cas10_end = int(cas10_global_start_end_coordinates[1]) + cas10_search_range
            print("The Cas10 start: " + str(cas10_start))
            print("The Cas10 end: " + str(cas10_end))
            try:
                gff_cas10 = gff.region(seqid=value["contig"], start = cas10_start, end = cas10_end, completely_within = True, featuretype = "CDS") #find corresponding coordinates from gff file
                print("The Cas10 coordinates: " + str(gff_cas10))
                gff_cas10 = list(gff_cas10)[0] #convert generator to list
                print("gff cas 10 type: " + str(type(gff_cas10)))
                print("Opened generator object: " + str(gff_cas10))

                gff_cas10 = gff.region(seqid=value["contig"], start = cas10_start, end = cas10_end, completely_within = True, featuretype = "CDS") #find corresponding coordinates from gff file
                gff_cas10 = list(gff_cas10)[0] #convert generator to list

            except Exception as errormsg:
                print("Could not pinpoint Cas10 genomic location: " + str(errormsg))
                #continue
            try:
                cas10_genomic_id = gff_cas10.attributes["Name"][0] #pick the code for the Cas10 protein from the gff object
                cas10_object = proteins[cas10_genomic_id]
                cas10_object.id = sample
                cas10_object.description = ""
                if len(cas10_object.seq) >= cas10_length_cutoff:
                    type_III_loci[key]["Cas10_boolean"] = True
                    type_III_loci[key]["Cas10_object"] = cas10_object
                    type_III_loci[key]["Cas10_seq"] = cas10_object.seq #find this protein in the original .faa file
                    cas10_cyclase_domains = getCas10CyclaseDomain(type_III_loci[key]["Cas10_seq"])

                    if len(cas10_cyclase_domains) > 0: #if GGDD motif (with one mismatch allowed) is found, mark it as so
                        print("GGDD found")
                        type_III_loci[key]["Cas10_cyclase"] = True
                        type_III_loci[key]["Cas10_cyclase_match_seq"] = cas10_cyclase_domains
            except Exception as errormsg:
                print("Could not extract Cas10: " + str(errormsg))
                #continue

        gene_info = {}
        gene_info["evalue"] = float(evaluelist[id])
        gene_info["flagged"] = False
        print("Gene E-value: " + str(gene_info["evalue"]) + " (" + str(gene) + ")")

    df = pd.DataFrame({"type_iii" : [str(value["type"])],
                        "CorA_in_iii_locus" : [type_III_loci[key]["has_CorA"]],
                        "CorA_seq" : [type_III_loci[key]["CorA_seq"]],
                        "Cas10_seq" : [type_III_loci[key]["Cas10_seq"]],
                        "Cas5_seq" : [type_III_loci[key]["Cas5_seq"]],
                        "host" : [sample]})
    
    crispr_locus_table = pd.concat([crispr_locus_table,df])


    # info_table = pd.DataFrame({"sample" : [sample],
    #                     "CorA" : type_III_loci[key]["has_CorA"],
    #                     "Cas10" : type_III_loci[key]["Cas10_boolean"],
    #                     "GGDD" : type_III_loci[key]["Cas10_cyclase"],
    #                     "GGDD_seq" : type_III_loci[key]["Cas10_cyclase_match_seq"],
    #                     "Subtype" : type_III_loci[key]["type"]
    #                     }, index=[0])

    print("Has CorA: " + str(type_III_loci[key]["has_CorA"]))
    info_table.at[0, 'CorA'] = type_III_loci[key]["has_CorA"]
    info_table.at[0, 'Cas10'] = type_III_loci[key]["Cas10_boolean"]
    info_table.at[0, 'GGDD'] = type_III_loci[key]["Cas10_cyclase"]
    info_table.at[0, 'GGDD_seq'] = type_III_loci[key]["Cas10_cyclase_match_seq"]
    info_table.at[0, 'Subtype'] = type_III_loci[key]["type"]


    finalTable = pd.concat([finalTable, crispr_locus_table])

    SeqIO.write(type_III_loci[key]["CorA_object"], out_cora, "fasta")
    SeqIO.write(type_III_loci[key]["CorA_object_cctyper"], out_cora_cctyper, "fasta")
    SeqIO.write(type_III_loci[key]["Cas10_object"], out_cas10, "fasta")
    SeqIO.write(type_III_loci[key]["Cas5_object"], out_cas5, "fasta")
    SeqIO.write(type_III_loci[key]["Cas7_object"], out_cas7, "fasta")


print("Finished. Final table:")


finalTable.to_csv(out_cora_info, sep = '\t', header = True, index = False)
info_table.to_csv(out_info, sep = '\t', header = False, index = False)