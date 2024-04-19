'''
Snakemake pipeline for Cas10 and CorA extraction and aligments (Haotian et al. 2023, NAR).
Actual trees with annotations are produced by the R script "tree.R" in the scripts folder.
'''

project = "run1"

base_path = "/media/volume/st_andrews/cas10_corA_2" + "/" + project
genomes_folder = "/media/volume/st_andrews/cas10_corA/genomes/bacteria/ncbi_dataset/data"

cas10_cluster_threshold = 0.90 #initial Cas10 clustering threshold
crispr_locus_interference_cutoff = 20 #cutoff for CRISPR loci interference completeness. Loci with less than this percentage of interference genes present are discarded

protein_clustering = str(config["protein_clustering"])
getGenomesBy = str(config["getGenomesBy"])
cas10_anchor = config["cas10_anchor"]
catyper_hmm_evalue = "1e-10"

hmm_msa_folder = "/media/volume/st_andrews/databases/cATyper/hmm/msa_050523"
hmm_database_folder = "/media/volume/st_andrews/databases/cATyper/hmm/050523" #contains the cctyper CorA hmm database
hmm_database_file = "effectors_050523.hmm" #hmm database for ca-associated proteins

if cas10_anchor == False:
    prefiltering_wildcards = "05_host_genomes"
    prefiltering_host_genomes = "06_host_genomes"
elif cas10_anchor == True:
    prefiltering_wildcards = "02_host_wildcards"
    prefiltering_host_genomes = "03_host_genomes"

def aggregate_crisprcas(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/07_cctyper/{i}/{i}_renaming.done", i=ivals)

def aggregate_host_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/02_genome_wildcards/{i}.txt", i=ivals)

def aggregate_cas10_sequences(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas10.faa", c=cvals)

    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_CorA.faa", c=cvals)

def aggregate_renamed_crisprs(wildcards):
    checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
    ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    return expand(base_path + "/07_cctyper/{j}/{j}_renaming.done", j=ivals)


def aggregate_download_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/06_host_genomes/{i}/{i}_genome.fna", i=ivals)


    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    print(expand(base_path + "/19_16S_seq/{j}_16S.fna", j=ivals))
    return expand(base_path + "/19_16S_seq/{j}_16S.fna", j=ivals)

def aggregate_cas10_booleans(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/062_genomes_cas10/{i}/{i}_cas10_boolean.tsv", i=ivals)

def aggregate_cas10_seq_postboolean(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/063_cas10_seq/{i}/{i}_cas10.faa", i=ivals)


def aggregate_cas10_sequences_prior_to_1st_clustering(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/063_cas10_seq/{i}/{i}_cas10.faa", i=ivals)


def aggregate_cas10_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/03_postfiltering_genome_wildcards/{j}.txt", j=jvals)


    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/30_clustertable/{j}/{j}.tsv", j=jvals)

def aggregate_typeIII_info(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_crispr_iii_info.tsv", c=cvals)

def aggregate_taxInfo(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/06_host_genomes/{j}/{j}_taxon.txt", j=ivals)

def aggregate_renamed(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/07_cctyper/{j}/{j}_renaming.done", j=ivals)



    checkpoint_output = checkpoints.align_matcher.get(**wildcards).output[0]
    kvals = glob_wildcards(os.path.join(checkpoint_output,"{k}.done")).k
    return expand(base_path + "/42_matching_trees_distance_matrices/{k}/{k}_correlation.tsv", k=kvals)

def aggregate_cATyper_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_results.tsv", c=cvals)

def aggregate_CorA_sequences_catyper(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/10_cATyper_hmm/{c}/CorA.faa", c=cvals)

def aggregate_CorA_sequences(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/CorA.faa", c=cvals)


def aggregate_unknowns(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/30_unknown_effectors/{c}/{c}_unknown_proteins.faa", c=cvals)

#Cas10 blast parameters
blast_headers = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sacc"
cas10_e_value = "1e-20"
corA_hmm_value = "1e-20"

print("Starting Cas10/CorA pipeline. Configs:\nProtein clustering: "+ str(protein_clustering)+ "\n")

rule all: 
    input: base_path + "/done"


genome_count = 51920
subsampling_seed = 666

genomes_json = os.path.join(genomes_folder,"assembly_data_report.jsonl")


rule write_down_genome_info:
    '''
    Outputs a text file genome names based on folders in genomes_folders
    '''
    output: base_path + "/01_genomelist/genomelist.txt"
    shell:
        """
        cd {genomes_folder}
        find * -maxdepth 0 -type d > {output}
        """

if getGenomesBy == "local":
    checkpoint expand_host_genomelist_to_wildcards:
        '''
        Takes as input a list of genome accessions separated by newline. These genomes are then turned into wildcards.
        In random mode, uses subsampling to randomly sample a set of genomes.
        '''
        input: rules.write_down_genome_info.output
        output: directory(base_path + "/" + prefiltering_wildcards)
        run:
            import pandas as pd
            if not os.path.exists(str(output)):
                os.makedirs(str(output))
            genomefiles = pd.read_csv(str(input), sep = "\t", header = None)
            #print("HELLOO " + genomefiles)
            genomes = [] #final list of accessions

            #Random mode. Subsamples n genomes from all genomes randomly with seed.
            if config["genome_mode"] == "random":
                subsample_genomes = genomefiles.sample(n = genome_count, random_state=subsampling_seed)
                genomes = subsample_genomes[0]

            elif config["genome_mode"] == "all":
                genomes = genomefiles[0]

            #Taxid id mode. NCBI datasets json is probed to find representatives of a taxid.
            elif config["genome_mode"] == "taxid":
                json_df = pd.read_json(genomes_json, lines=True)
                json_df = json_df.join(pd.json_normalize(json_df['organism'])).drop('organism', axis='columns') #expand pandas dictionary column into new columns in the actual df
                taxidlistlist = []

                with open(config["taxidlistfile"]) as f: #read taxid list file
                    taxids_comma = f.read()
                    idlist = taxids_comma.split(",")

                taxidlist = pd.DataFrame(idlist) #convert list to df
                #taxidlist = pd.read_csv(config["taxidlistfile"], sep=",", header=None)
                taxidlist.columns = ["taxId"] #add column to the df
                json_df['taxId'] = json_df['taxId'].astype("string")
                taxidlist['taxId'] = taxidlist['taxId'].astype("string")
                chosen = pd.merge(json_df, taxidlist, on="taxId", how="inner")
                #chosen = json_df.loc[json_df["taxId"] == int(config["taxid"])]
                print("Shape of subsampled dataframe using TaxIDs from " + str(config["taxidlistfile"]) + ": " + str(chosen.shape))
                genomes = chosen["accession"].values.tolist()


            #add genome to list if .faa, .gff and .fna files are present
            for i in genomes:
                if (os.path.isfile(os.path.join(genomes_folder,i,"protein.faa")) and (os.path.isfile(os.path.join(genomes_folder,i,"genomic.gff")))): #make sure it's annotated
                    for fname in os.listdir(os.path.join(genomes_folder,i)): #checking fna is more complicated because the filename is not consistent
                        if fname.endswith('.fna'):
                            with open(str(output) + "/" + str(i) + '.txt', 'w') as f:
                                f.write(str(output))
        


    rule download_genomes:
        '''
        Based on the wildcards from expand_host_genomelist_to_wildcards,
        makes symlinks for the desired genomes from the genomes_folder
        into the working folder (prefiltering_host_genomes).
        '''
        input:
            ivalue = base_path + "/" + prefiltering_wildcards + "/{i}.txt",
        output: 
            fna = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_genome.fna",
            faa = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa",
            gff = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_features.gff"
        params: #these are mostly deprecated
            folder = base_path + "/" + prefiltering_host_genomes,
        conda: "envs/ncbidownload.yaml"
        #threads: 4
        log:
            out = base_path + "/logs/" + prefiltering_host_genomes + "/{i}.log",
            err = base_path + "/logs/" + prefiltering_host_genomes + "/{i}.err"
        shell:
            '''
            echo "Creating symlinks for chosen genomes"
            ln -s {genomes_folder}/{wildcards.i}/*.fna {output.fna}
            ln -s {genomes_folder}/{wildcards.i}/*.faa {output.faa}
            ln -s {genomes_folder}/{wildcards.i}/*.gff {output.gff}
            '''


elif getGenomesBy == "remote":
    species_space = str(config["species"])
    checkpoint getHosts:
        '''
        Downloads desired host genomes from NCBI using the wonderful python program ncbi-genome-download.
        Files are then gunzipped and named after their parent folder
        '''
        output: directory(base_path + "/" + prefiltering_host_genomes)
        conda: "envs/ncbidownload.yaml"
        params:
            parallels = "40",
            folder = base_path + "/06_host_genomes"
        shell:
            '''
            mkdir -p {params.folder}
            cd {params.folder}
            count=$(ncbi-genome-download --dry-run --genera "{species_space}" bacteria | wc -l)
            echo "Will download $count genomes from {species_space}. Note that the download will likely fail if the dataset is large (>1000 genomes). In such a case, just restart the script."
            printf 'Do you want to continue? (y/n)? '
            read -p "Do you want to continue? (y/n) " -n 1 -r
            if [[ $REPLY =~ ^[Yy]$ ]] ;then
                ncbi-genome-download --genera "{species_space}" bacteria --parallel {params.parallels} --formats fasta,protein-fasta,gff,assembly-report
                mv refseq/bacteria/* {params.folder}
                find {params.folder} -type f -name '*.gz' | tqdm | xargs gunzip
                echo renaming
                rename 's!(.*)/(.*)genomic\.fna!$1/$1_genome.fna!' */*genomic.fna
                rename 's!(.*)/(.*)protein\.faa!$1/$1_proteins.faa!' */*protein.faa
                rename 's!(.*)/(.*)genomic\.gff!$1/$1_features.gff!' */*genomic.gff
                rename 's!(.*)/(.*)assembly_report\.txt!$1/$1_report.txt!' */*assembly_report.txt
            else
                exit 1
            fi

            '''


if cas10_anchor == True: #if we filter analyzable genomes by the presence of Cas10. Uses Cas10 HMM profiles from CCtyper.
    print("Running Cas10 stuff")
    rule Cas10_genomes:
        '''
        First filter all sequences to get only >500 AA long proteins. This prevents the inclusion of
        truncated Cas10s or those that are cut artifically due to contig ending. Also reduces HMM searches.
        
        Then:
        Searches the filtered proteins of host against a user-specified preprepared hmm profile.
        Strips file of "#" -rows and checks if there is anything left. If there is, the gene has been found.
        Boolean output format:
        i,True/False,cas10_accession
        '''
        input:
            proteins = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa"
        output:
            hmm = base_path + "/062_genomes_cas10/{i}/{i}_cas10.tsv",
            boolean = base_path + "/062_genomes_cas10/{i}/{i}_cas10_boolean.tsv"
        params:
            proteins_filt = base_path + "/062_genomes_cas10/{i}/{i}_proteins_lengthfiltered.faa",
            out = base_path + "/062_genomes_cas10/{i}/{i}_temp.out",
            rows1 = base_path + "/062_genomes_cas10/{i}/{i}_temp_rows_1.out",
            cas10_db = "/media/volume/st_andrews/cas10_corA_2/profiles/Profiles/all_cas10s.hmm",
            rows = base_path + "/062_genomes_cas10/{i}/{i}_temp_rows.out",
            headers = base_path + "/062_genomes_cas10/{i}/{i}_headers.out",
            all_data = base_path + "/062_genomes_cas10/{i}/{i}_all_data.out",
        conda: "envs/hmmer.yaml"
        shell:
            '''
            echo "Running rule Cas10_genomes"
            echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm}
            if [ -s "{input.proteins}" ]; then
                cat {input.proteins} | seqkit seq -m 500 > {params.proteins_filt}
                if [ -s "{params.proteins_filt}" ]; then
                    hmmscan --domtblout {params.out} --cpu 1 -E {corA_hmm_value} {params.cas10_db} {params.proteins_filt} &> /dev/null
                    grep -v "#" {params.out} > {params.rows}||:
                    head -1 {params.rows} > {params.rows1}
                    echo "id,cas10_boolean,cas10_acc" > {output.boolean}
                    if [ -s {params.rows1} ]; then
                        cat {params.rows1} >> {output.hmm}
                        ACC=$(awk -F ' ' '{{print $4}}' {params.rows1})
                        echo "{wildcards.i},True","${{ACC}}" > {output.boolean}
                    else
                        echo "{wildcards.i},False,-" > {output.boolean}
                    fi
                    touch {output.hmm}
                else
                    echo "{wildcards.i},False,-" > {output.boolean}
                fi
            else
                echo "{wildcards.i},False,-" > {output.boolean}
            fi
            '''

    rule extract_cas10_sequence:
        input:
            boolean = rules.Cas10_genomes.output.boolean,
            proteins = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa"
        output: base_path + "/063_cas10_seq/{i}/{i}_cas10.faa"
        shell: 
            '''
            CAS10_ID=$(awk -F ',' '{{print $3}}' {input.boolean})
            if [ "$CAS10_ID" != "-" ]; then
                echo ">"{wildcards.i} > {output}
                seqkit grep -r -p ${{CAS10_ID}} {input.proteins} | seqkit seq -w 0 | tail -n +2 >> {output}
            else
                touch {output}
            fi
            '''

    rule concat_cluster_cas10_proteins:
        '''
        Takes all Cas10 sequences, concatenates them into one fasta and
        clusters this fasta at cas10_cluster_threshold (0-1) similarity.
        Greps only representative lines "marked by *"
        A python script extracts the accession from these lines into output.reps
        '''
        input: aggregate_cas10_seq_postboolean
        output:
            all_cas10 = base_path + "/064_cas10_clusters/cas10_all.faa",
            clusters = base_path + "/064_cas10_clusters/cas10_clust.faa.clstr",
            proteins = base_path + "/064_cas10_clusters/cas10_clust.faa",
            reps = base_path + "/064_cas10_clusters/cas10_unique_genomes.txt",
        params:
            clusterlines = base_path + "/064_cas10_clusters/clusterlines.txt"
        threads: 40
        conda: "envs/trees.yaml
        shell:
            '''
            echo "concatenating cas10 sequences for clustering"
            find '{base_path}/063_cas10_seq' -maxdepth 2 -type f -wholename '*/*_cas10.faa' -print0 | xargs -0 cat > {output.all_cas10}
            echo clustering
            cd-hit -i {output.all_cas10} -o {output.proteins} -c {cas10_cluster_threshold} -n 5 -d 0 -M 16000 -T {threads}
            grep "*" {output.clusters} > {params.clusterlines}
            python scripts/getClusterRep.py --input {params.clusterlines} --output {output.reps}
            '''

    checkpoint expand_host_genomelist_to_wildcards_postcas10:
        '''
        '''
        input:
            genomes = rules.concat_cluster_cas10_proteins.output.reps
        output: directory(base_path + "/03_postfiltering_genome_wildcards")
        run:
            print("Expanding host list to wildcards...")
            if not os.path.exists(str(output)):
                os.makedirs(str(output))

            with open(str(input.genomes)) as filehandle:
                genomes = [line.rstrip('\n') for line in filehandle]

            for i in genomes:
                sample = str(str(i).split(",")[0])
                with open(str(output) + "/" + sample + '.txt', 'w') as f:
                    f.write(sample)



    rule download_genomes_postCas10:
        '''
        After filtering by presence of Cas10, create new symlinks for such genomes.
        '''
        input: 
            genomes = base_path + "/03_postfiltering_genome_wildcards/{j}.txt",
        output: 
            fna = base_path + "/06_host_genomes/{j}/{j}_genome.fna",
            faa = base_path + "/06_host_genomes/{j}/{j}_proteins.faa",
            gff = base_path + "/06_host_genomes/{j}/{j}_features.gff"
        params: #these are mostly deprecated
            folder = base_path + "/06_host_genomes",
            zipped = base_path + "/06_host_genomes/{j}.zip",
            original_fna = base_path + "/06_host_genomes/{j}/*.fna", #for checking successful download
            original_prot = base_path + "/06_host_genomes/{j}/protein.faa",#for checking successful download
            original_gff = base_path + "/06_host_genomes/{j}/genomic.gff", #for checking successful download
        conda: "envs/ncbidownload.yaml"
        #threads: 4
        log:
            out = base_path + "/logs/06_host_genomes/{j}.log",
            err = base_path + "/logs/06_host_genomes/{j}.err"
        shell:
            '''
            ln -s {genomes_folder}/{wildcards.j}/*.fna {output.fna}
            ln -s {genomes_folder}/{wildcards.j}/*.faa {output.faa}
            ln -s {genomes_folder}/{wildcards.j}/*.gff {output.gff}
            '''

rule getTaxInfo:
    input: base_path + "/03_postfiltering_genome_wildcards/{j}.txt"
    output:
        taxon = base_path + "/06_host_genomes/{j}/{j}_taxon.txt"
    run:
        import pandas as pd
        json_df = pd.read_json(genomes_json, lines=True)
        json_df = json_df.join(pd.json_normalize(json_df['organism'])).drop('organism', axis='columns') #expand pandas dictionary column into new columns in the actual df
        json_df = json_df.join(pd.json_normalize(json_df['averageNucleotideIdentity'])).drop('averageNucleotideIdentity', axis='columns') #expand pandas dictionary column into new columns in the actual df 
        try:
            row = json_df.loc[json_df["accession"] == wildcards.j]
            species = row['submittedSpecies'].iloc[0].strip("][")
            genus = species.split(" ")[0]
            with open(output.taxon, "w") as f:
                f.write("Sample\tGenus\tSpecies\n")
                f.write(wildcards.j + "\t" + genus + "\t" + species + "\n")
        except:
            with open(output.taxon, "w") as f:
                            f.write(wildcards.j + "\tUnknown\tUnknown\n")
            pass


rule concat_taxInfo:
    input: aggregate_taxInfo
    output: base_path + "/06_host_genomes/taxInfo.txt"
    shell:
        """
        echo "Sample\tgenus\tspecies\n" > {output}
        find '{base_path}/06_host_genomes/' -maxdepth 2 -type f -wholename '*_taxon.txt' -print0 | xargs -0 cat >> {output}
        """

rule CRISPRCasTyper:
    '''
    Types CRISPR-Cas loci in the genomes based on fasta.
    Takes output of the checkpoint. Automatically scales to the number of files
    produced by the checkpoint, so no need to call for the aggregator function here.

    '''
    input: base_path + "/06_host_genomes/{j}/{j}_genome.fna"
    output: base_path + "/07_cctyper/{j}/{j}.done"
    conda: "envs/CRISPRCasTyper.yaml"
    log: base_path + "/logs/07_cctyper/{j}/{j}_cctyper.log"
    params:
        outdir = base_path + "/07_cctyper/{j}"
    shell:
        '''
        rm -rf {params.outdir}
        cctyper '{input}' '{base_path}/07_cctyper/{wildcards.j}' --prodigal single 2>{log}
        touch '{output}'
        '''

rule CRISPRCasTyper_rename:
    '''
    This uses the CRISPR-Cas.tab file to include only fully functional CRISPR-Cas loci
    The sed adds assembly information to the results file, if it exists (which is originally based on just the nucleotide accessions)
    '''
    input: rules.CRISPRCasTyper.output
    output: base_path + "/07_cctyper/{j}/{j}_renaming.done"
    params:
        outdir = base_path + "/07_cctyper/{j}"
    shell:
        '''
        if test -e "{params.outdir}/CRISPR_Cas.tab";then
            mv {params.outdir}/CRISPR_Cas.tab {params.outdir}/CRISPR_Cas_old.tab
            sed '1s/$/\tassembly/; 2,$s/$/\t{wildcards.j}/' {params.outdir}/CRISPR_Cas_old.tab > {params.outdir}/CRISPR_Cas.tab
            rm {params.outdir}/CRISPR_Cas_old.tab
        fi
        touch {output}
        '''

rule concat_renamed_crisprs:
    '''
    This rule aggregates the outputs of the CRISPRCasTyper_rename rule.
    It is necessary to do this because the CRISPRCasTyper_rename rule is called by the checkpoint,
    and the checkpoint does not automatically aggregate the outputs of the rule.
    '''
    input: aggregate_renamed
    output: base_path + "/07_cctyper/renaming.done"
    shell:
        '''
        touch {output}
        '''

checkpoint type_iii_wildcarder:
    '''
    Extracts only type III loci from the CCTyper output
    This checkpoint examines outputs from CCTyper and creates new wildcards based on complete loci.
    This marks a point where the pipeline no longer looks at individual strains, but rather at the CRISPR-Cas loci individually.
    As outputs, this checkpoint extracts locus-specific information from CCTyper outputs (CRISPR-Cas.tab and cas_operons.tab)
    '''
    input:
        cctyper_done = aggregate_renamed,
        tax_info = rules.concat_taxInfo.output
    output: directory(base_path + "/071_cctyper_loci")
    params:
        interference_cutoff = crispr_locus_interference_cutoff, #0-100. If the interference score is lower than this, the locus is discarded
        cctyper_folder = base_path + "/07_cctyper",
    conda: "envs/pandas.yaml"
    shell:
        '''
        python3 scripts/loci_wildcarder.py --input_folder {params.cctyper_folder} --output_folder {output} --interference_cutoff {params.interference_cutoff} --tax_info {input.tax_info}
        '''


rule cATyper_hmm:
    '''
    Characterizes the CRISPR-Cas locus regarding CorA and/or other effectors

    Steps:
        1. Get contig of locus
        2. Use the gff file to extract names of proteins in that contig. Also extract their coordinates to a csv file.
        3. Create a new protein multifasta that contains said proteins
        4. Run hmmscan on that multifasta
        5. Go through the hmm hits and extract the coordinates from the csv file
        6. Mark down boolean values on whether the effector is near the CRISPR-Cas locus or not

    I reuse the same python script before and after the hmm run. The mode is set to "pre_hmm" before the hmm run, and "post_hmm" after the hmm run.
    '''
    input:
        #this input is just for obtaining the wildcard. The sample name is derived from the wildcard
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
    output:
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        catyper = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_results.tsv",
        cora = base_path + "/10_cATyper_hmm/{c}/CorA.faa",
        #unknown_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_unknown_proteins.faa", #new in version 9
        hmm_targets = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_targets.tsv",    
        hmm_rows = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv",
        temp_rows = base_path + "/10_cATyper_hmm/{c}/{c}_temp_rows.out",
    params:
        outdir = base_path + "/10_cATyper_hmm/{c}",
        hmm_msa_folder = hmm_msa_folder,
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_profile = hmm_database_folder + "/" + hmm_database_file,
        temp_hmm = base_path + "/10_cATyper_hmm/{c}/{c}_temp.out",
        temp_rows = base_path + "/10_cATyper_hmm/{c}/{c}_temp_rows.out",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/10_cATyper_hmm/logs/{c}.out",
        err = base_path + "/10_cATyper_hmm/logs/{c}.err",
        out2 = base_path + "/10_cATyper_hmm/logs/{c}.out2",
        err2 = base_path + "/10_cATyper_hmm/logs/{c}.err2",
    shell:
        '''
        python scripts/catyper_prepper_9.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode pre_hmm 2> {log.err} 1> {log.out}
        echo "Running hmmscan" >> {log.out}
        hmmscan --domtblout {params.temp_hmm} --cpu 8 -E {catyper_hmm_evalue} {params.hmm_profile} {output.contig_proteins}  &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.c}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.c}" >> {log.out}
            touch {output.hmm_rows}
        fi
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*.fa > {output.hmm_targets}
        echo "Running cATyper" >> {log.out}
        python scripts/catyper_prepper_9.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {output.hmm_rows} --catyper_out {output.catyper}  2> {log.err2} 1> {log.out2}
        touch {output.cora}
        '''

rule concatenate_cATyper_hmm:
    input: aggregate_cATyper_hmm
    output: base_path + "/10_cATyper_hmm/cATyper_all.tsv"
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {base_path}/10_cATyper_hmm/*/*_cATyper_results.tsv > {output}
        touch {output}
        """


rule typeIII_characterizer:
    '''
    This rule takes the output of the typeIII_wildcarder rule and searches for:
        - Cas10
        - Cas5
        - Cas7
        - CorA
        - Unknown genes
    The outputs are fasta sequences for each protein as well as an info table (.tsv)
    that concatenates information on said genes
    '''
    input:
        cctyper = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv", #the cas_operons file is trimmed down to only include this specific {c} locus
        #gff = base_path + "/06_host_genomes/{j}/{j}_features.gff",
        #proteins = base_path + "/06_host_genomes/{j}/{j}_proteins.faa",
    output:
        Cas10_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas10.faa",
        Cas5_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas5.faa",
        Cas7_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas7.faa",
        CorA_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_CorA.faa",
        info = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_crispr_iii_info.tsv",
    params:
        this_folder = base_path + "/09_crispr_iii_CorA/loci/{c}",
        outputfolder = base_path + "/09_crispr_iii_CorA/loci/{c}",
        sample_folder = base_path + "/06_host_genomes/", #this is used to find the gff and proteins files that are strain-, not locus-specific
        cctyper_folder = base_path + "/07_cctyper"
    conda: "envs/gff_utils.yaml"
    log:
        out = base_path + "/09_crispr_iii_CorA/logs/{c}.out",
        err = base_path + "/09_crispr_iii_CorA/logs/{c}.err"
    shell:
        '''
        python3 scripts/type_iii_effector_finder_2.0.py --locus_id {wildcards.c} --sample_folder {params.sample_folder} --this_folder {params.this_folder} --outputfolder {params.outputfolder} --cas_operons {input.cctyper} --info_out {output.info} --cctyper_path {params.cctyper_folder} 2> {log.err} 1> {log.out}
        touch {output.Cas10_fasta}
        touch {output.Cas5_fasta}
        touch {output.Cas7_fasta}
        touch {output.CorA_fasta}
        '''

rule concatenate_type_iii_info:
    input: aggregate_typeIII_info
    output: base_path + "/09_crispr_iii_CorA/loci/type_iii_info.tsv"
    shell:
        '''
        echo "Cas10\tCas5\tCas7\tCorA\tLocus\tSample\tCas10_GGDD\tCas10_GGDD_seq\tUnknown_genes\tSubtype" > {output}
        find '{base_path}/09_crispr_iii_CorA/loci' -maxdepth 2 -type f -wholename '*/*_crispr_iii_info.tsv' -print0 | xargs -0 cat >> {output}
        '''

#a rule that uses the output file cora_type_iii_info from the rule above filter samples with CorA and further divide them into CRISPR-Cas subtypes
rule cora_plot_extractor:
    '''
    Gets cctyper generated plots for each sample with a CorA and a type III CRISPR-Cas system.
    '''
    input: 
        info = rules.concatenate_type_iii_info.output[0]
    output:
        done = base_path + "/xtra1_cora_iii_loci_plots/done.done",
    run:
        import pandas as pd
        import shutil
        #using pandas, get the info from the cora_type_iii_info file and filter out samples that do not have a CorA
        cora_type_iii_info = pd.read_csv(input.info, sep = "\t")
        cora_type_iii_info = cora_type_iii_info[cora_type_iii_info["CorA"] == True]

        #create output folder if it does not exist
        if not os.path.exists(base_path + "/xtra1_cora_iii_loci_plots"):
            os.makedirs(base_path + "/xtra1_cora_iii_loci_plots")

        #For each sample in the pandas dataframe, extract the corresponding plot from the cctyper folder (path is base_path + "/07_cctyper/{j}/plot.png) and copy it to the output folder in a CRISPR-Cas subtype subfolder (e.g. III-A or III-B) depending on the Subtype column in the cora_type_iii_info file
        for index, row in cora_type_iii_info.iterrows():
            sample = row["Sample"]
            subtype = row["Subtype"]
            print(sample + "," + subtype)
            #check if subtype does not contain the substring "Hybrid" (this is because some samples have a hybrid subtype, e.g. III-A/B, and we don't want to create a III-A/B folder)
            if "Hybrid" not in subtype:
                #if the subtype folder does not exist, create it
                if not os.path.exists(base_path + "/xtra1_cora_iii_loci_plots/" + subtype):
                    os.makedirs(base_path + "/xtra1_cora_iii_loci_plots/" + subtype)
                shutil.copyfile(base_path + "/07_cctyper/" + sample + "/plot.png", base_path + "/xtra1_cora_iii_loci_plots/" + subtype + "/" + sample + "_plot.png")

        #create the done.done file to indicate that the rule has finished running
        open(output.done, "w").close()


rule CorA_concatenate:
    '''
    This version works on the cATyper outputs instead of CCTyper outputs
    '''
    input: aggregate_CorA_sequences_catyper
    output: base_path + "/09_crispr_iii_CorA/CorAs.faa"
    shell:
        '''
        cat {input} > {output} 
        '''

rule CorA_cluster:
    '''

    '''
    input: rules.CorA_concatenate.output
    output:
        proteins = base_path + "/15_CorA_cluster/CorA_cluster.faa",
        clusterinfo = base_path + "/15_CorA_cluster/CorA_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/15_CorA_cluster/logs/CorA_align.out",
        err = base_path + "/15_CorA_cluster/logs/CorA_align.err"
    threads: 40
    shell:
        '''
        cd-hit -i {input} -o {output.proteins} -c 0.90 -n 5 -d 0 -M 16000 -T {threads}
        '''


rule CorA_align:
    input: rules.CorA_cluster.output.proteins
    output: base_path + "/16_CorA_align/CorA_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/16_CorA_align/logs/CorA_align.out",
        err = base_path + "/16_CorA_align/logs/CorA_align.err"
    threads: 40
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule CorA_align_unclustered:
    input: rules.CorA_concatenate.output
    output: base_path + "/16_CorA_align_unclustered/CorA_alignment_unclustered.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/16_CorA_align/logs/CorA_align.out",
        err = base_path + "/16_CorA_align/logs/CorA_align.err"
    threads: 40
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule CorA_tree:
    '''

    '''
    input: rules.CorA_align.output
    output: base_path + "/17_CorA_tree/CorA_tree.txt",
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''

rule CorA_tree_unclustered:
    '''

    '''
    input: rules.CorA_align_unclustered.output
    output: base_path + "/17_CorA_tree_unclustered/CorA_tree_unclustered.txt",
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''


rule Cas10_concatenate:
    input: aggregate_cas10_sequences
    output: base_path + "/09_crispr_iii_CorA/cas10s.faa"
    shell:
        '''
        cat {input} > {output} 
        '''

rule Cas10_cluster:
    '''
    '''
    input: 
        cas10 = rules.Cas10_concatenate.output
    output:
        proteins = base_path + "/10_Cas10_cluster/cas10_cluster.faa",
        clusterinfo = base_path + "/10_Cas10_cluster/cas10_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/10_Cas10_align/logs/cas10_align.out",
        err = base_path + "/10_Cas10_align/logs/cas10_align.err"
    threads: 40
    shell:
        '''
        echo {protein_clustering}
        if [ {protein_clustering} = "True" ]; then
            cd-hit -i {input.cas10} -o {output.proteins} -c 0.99 -n 5 -d 0 -M 16000 -T {threads}
        else
            cp {input.cas10} {output.proteins}
            touch {output.clusterinfo}
        fi
        '''


rule Cas10_align:
    '''
    
    '''
    input: rules.Cas10_cluster.output.proteins
    output: base_path + "/11_Cas10_align/cas10_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/11_Cas10_align/logs/cas10_align.out",
        err = base_path + "/11_Cas10_align/logs/cas10_align.err"
    threads: 40
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule Cas10_tree:
    '''

    '''
    input: rules.Cas10_align.output
    output: base_path + "/12_Cas10_tree/cas10_tree.txt",
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''

rule final:
    input:
        tax = base_path + "/06_host_genomes/taxInfo.txt",
        tree_Cas10 = rules.Cas10_tree.output,
        type_iii_info = rules.concatenate_type_iii_info.output,
        concat_taxInfo = rules.concat_taxInfo.output,
        catyper = rules.concatenate_cATyper_hmm.output,
        tree_CorA = rules.CorA_tree.output,
        tree_CorA_unclustered = rules.CorA_tree_unclustered.output,
    output: base_path + "/done"
    shell:
        '''
        touch {output}
        mkdir -p /media/volume/st_andrews/cas10_corA_2/R
        cp /media/volume/st_andrews/cas10_corA_2/test1/12_Cas10_tree/cas10_tree.txt /media/volume/st_andrews/cas10_corA_2/test1/R
        cp /media/volume/st_andrews/cas10_corA_2/test1/06_host_genomes/taxInfo.txt /media/volume/st_andrews/cas10_corA_2/test1/R
        cp /media/volume/st_andrews/cas10_corA_2/test1/10_cATyper_hmm/cATyper_all.tsv /media/volume/st_andrews/cas10_corA_2/test1/R
        cp /media/volume/st_andrews/cas10_corA_2/test1/09_crispr_iii_CorA/loci/type_iii_info.tsv /media/volume/st_andrews/cas10_corA_2/test1/R
        '''
