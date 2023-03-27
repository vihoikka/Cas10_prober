'''
This version uses predownloaded genomes and simply characterizes all
for CRISPR-Cas
'''

protein_clustering = str(config["protein_clustering"])
getGenomesBy = str(config["getGenomesBy"])
cas10_anchor = config["cas10_anchor"]

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
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).i
    return expand(base_path + "/09_crispr_iii_CorA/{j}/{j}_Cas10.faa", j=ivals)

def aggregate_cas5_sequences(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/09_crispr_iii_CorA/{j}/{j}_Cas5.faa", j=ivals)

def aggregate_cas7_sequences(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/09_crispr_iii_CorA/{i}/{i}_Cas7.faa", i=ivals)


def aggregate_CorA_sequences(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/09_crispr_iii_CorA/{j}/{j}_CorA.faa", j=ivals)

def aggregate_download_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/06_host_genomes/{i}/{i}_genome.fna", i=ivals)

def aggregate_16s(wildcards):
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

def aggregate_clustertables(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/30_clustertable/{j}/{j}.tsv", j=jvals)

def aggregate_typeIII_info(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/09_crispr_iii_CorA/{j}/{j}_crispr_iii_info.tsv", j=ivals)

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


# def aggregate_matched_trees(wildcards):
#     checkpoint_output = checkpoints.align_matcher.get(**wildcards).output[0]
#     kvals = glob_wildcards(os.path.join(checkpoint_output,"{k}.done")).j
#     return expand(base_path + "/30_clustertable/{j}/{j}.tsv", j=jvals)

def aggregate_dms(wildcards):
    checkpoint_output = checkpoints.align_matcher.get(**wildcards).output[0]
    kvals = glob_wildcards(os.path.join(checkpoint_output,"{k}.done")).k
    return expand(base_path + "/42_matching_trees_distance_matrices/{k}/{k}_correlation.tsv", k=kvals)

base_path = "/media/volume/st_andrews/cas10_corA_2"

#Cas10 blast parameters
blast_headers = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sacc"
cas10_e_value = "1e-20"
corA_hmm_value = "1e-20"



ecoli_k12_16s = "/media/volume/st_andrews/cas10_corA/genomes/16s_ecoli_K12.fasta"

print("Starting Cas10/CorA pipeline. Configs:\nProtein clustering: "+ str(protein_clustering)+ "\n")

rule all: 
    input: base_path + "/done"
    #base_path + "/063_cas10_boolean/cas10_boolean.tsv"
    
    #base_path + "/19_16S_seq/cluster_16s.txt"
    #base_path + "/21_16s_tree/16s_tree.txt"
    #base_path + "/12_Cas10_tree/cas10_tree.txt"
    #base_path + "/17_CorA_tree/CorA_tree.txt"
    #aggregate_CorA
    #aggregate_crisprcas
    #base_path + "/04_genomelist/genomelist.txt"
    #base_path + "/04_host_genome_files/hostlist.txt"
    #aggregate_actual_genomes
    #base_path + "/01_cas10/cas10_blast.out"


genomes_folder = "/media/volume/st_andrews/cas10_corA/genomes/bacteria/ncbi_dataset/data/"
genome_count = 50926
subsampling_seed = 666

genomes_json = os.path.join(genomes_folder,"assembly_data_report.jsonl")





rule write_down_genome_info:
    '''
    Outputs a text file genome names based on folders in genomes_folders
    '''
    output: base_path + "/01_genomelist/genomelist.txt"
    shell:
        """
        ls {genomes_folder} > {output}
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
            ln -s {genomes_folder}{wildcards.i}/*.fna {output.fna}
            ln -s {genomes_folder}{wildcards.i}/*.faa {output.faa}
            ln -s {genomes_folder}{wildcards.i}/*.gff {output.gff}
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
        clusters this fasta at 99% similarity.
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
        shell:
            '''
            echo "concatenating cas10 sequences for clustering"
            find '{base_path}/063_cas10_seq' -maxdepth 2 -type f -wholename '*/*_cas10.faa' -print0 | xargs -0 cat > {output.all_cas10}
            echo clustering
            cd-hit -i {output.all_cas10} -o {output.proteins} -c 0.99 -n 5 -d 0 -M 16000 -T {threads}
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
            ln -s {genomes_folder}{wildcards.j}/*.fna {output.fna}
            ln -s {genomes_folder}{wildcards.j}/*.faa {output.faa}
            ln -s {genomes_folder}{wildcards.j}/*.gff {output.gff}
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

rule makeBlastDBnt:
    input: base_path + "/06_host_genomes/{j}/{j}_genome.fna"
    output:
        blastdb = base_path + "/061_host_genomes_blastDB/genomes/{j}/{j}.nhr"
    params:
        output = base_path + "/061_host_genomes_blastDB/genomes/{j}/{j}"
    conda: "envs/blast.yaml"
    shell:
        '''
        makeblastdb -in {input} -out {params.output} -dbtype nucl
        '''

rule blast_16s:
    input:
        db = rules.makeBlastDBnt.output.blastdb,
        query = ecoli_k12_16s
    output:
        blast = base_path + "/18_16S_blast/{j}.out"
    conda: "envs/blast.yaml"
    params:
        db_base_name = base_path + "/061_host_genomes_blastDB/genomes/{j}/{j}",
        temp_blast = base_path + "/061_host_genomes_blastDB/genomes/{j}/temp.blast"
    shell:
        '''
        blastn -db {params.db_base_name} -task blastn -query {input.query} -outfmt "6 qseqid qlen sseqid stitle staxids saccver sstart send evalue length pident mismatch gapopen gaps sstrand qcovs qcovhsp bitscore" -evalue 0.01 -word_size 11 -dust no > {output.blast}
        echo -e 'Query_id\tQuery_length\tSubject_id\tSubject_sci_title\tSubject_taxID\tSubject_accession_ID_version\tSubject_start\tSubject_end\tEvalue\tLength\tPercentage_identical\tMismatches\tOpen_gaps\tGaps\tSubject_strand\tCoverageSubject\tCoverageHSP\tBitscore' | cat - {output.blast} > {params.temp_blast} && mv {params.temp_blast} {output.blast}
        '''

rule extract16s:
    '''
    NOTE: Here I tested wrapping the Run python script into Shell. It's a bit finicky and hard to read, but works.
    The indentations here are important.
    '''
    input:
        blast = rules.blast_16s.output.blast,
        fna = base_path + "/06_host_genomes/{j}/{j}_genome.fna"
    output:
        seq_16S = base_path + "/19_16S_seq/{j}_16S.fna"
    conda: "envs/trees.yaml"
    shell:
        """
        python -c "import pandas as pd
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
blast = pd.read_csv('{input.blast}', sep = '\t', header = 0)
if blast.size > 0:
    best_hit = blast.loc[blast['Bitscore'].idxmax()]
    print(best_hit)
    coordinates = [best_hit['Subject_start'],best_hit['Subject_end']]

    genome = SeqIO.to_dict(SeqIO.parse('{input.fna}', 'fasta')) #read fasta into biopython object and convert to dictionary
    contig_16s = genome[best_hit['Subject_id']] #access dictionary without looping the list of accessions. Faster!
    seq_16s = contig_16s[coordinates[0]:coordinates[1]]

    if coordinates[0] > coordinates [1]: #if start is larger than end, then hit is on reverse strand. Flip them around for clarity.
        start = coordinates[0]
        end = coordinates[1]
        coordinates[0] = end
        coordinates[1] = start
        seq_16s = contig_16s[coordinates[0]:coordinates[1]].reverse_complement()
    else: #if not, just take the sequence
        seq_16s = contig_16s[coordinates[0]:coordinates[1]]

    seq_16s.id = '{wildcards.j}'
    seq_16s.description = '16S'

    SeqIO.write(seq_16s, '{output.seq_16S}', 'fasta')
"
        touch {output.seq_16S}
        """
        
rule concatenate_16s:
    input: aggregate_16s
    output: base_path + "/19_16S_seq/16s.fna.aggregate"
    shell:
        '''
        find '{base_path}/19_16S_seq/' -maxdepth 2 -type f -wholename '*.fna' -print0 | xargs -0 cat > {output}
        '''
        #cat {input} > {output}
        #printf '{input}' | xargs -r0 cat > {output}

rule cluster_16s:
    input: 
        concat_16s = rules.concatenate_16s.output,
    output: base_path + "/19_16S_seq/cluster_16s.txt",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/19_16S_seq/logs/cluster_16s.out",
        err = base_path + "/19_16S_seq/logs/cluster_16s.err"
    threads: 40
    shell:
        '''
        cd-hit -i {input.concat_16s} -o {output} -c 0.99 -n 5 -d 0 -M 160000 -T {threads}
        '''


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
    UPDATED VERSION. This uses the CRISPR-Cas.tab file to include only fully functional CRISPR-Cas loci
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

rule CorA_finder:
    '''
    Searches all proteins of host against a user-specified preprepared hmm profile.
    Strips file of "#" -rows and checks if there is anything left. If there is, the gene has been found.
    Result in a tab-separated file with just header and one row.
    Header = host\tcustom_gene_profile
    Row = sample\tyes or no

    NOTE: you learned something important here. If the outcome of grep is an empty file, it returns a non-zero exit.
    This causes the pipeline to fail without any error message. You spent hours in debugging this. The solution is
    to add ||: after grepping to catch the error. As a result of this adventure, this shell script is needlessly complicated.
    Prettify it someday if desired.
    '''
    input:
        proteins = base_path + "/06_host_genomes/{j}/{j}_proteins.faa"
    output:
        #boolean =  base_path + "/08_CorA_finder/{j}/{j}_CorA_boolean.tsv",
        hmm = base_path + "/08_CorA_finder/{j}/{j}_CorA_hmm.tsv"
    params:
        out = base_path + "/08_CorA_finder/{j}/{j}_temp.out",
        corA_hmm_db = "/media/volume/st_andrews/cas10_corA_2/profiles/Profiles/all_CorA.hmm", #from cctyper
        #corA_hmm_db = "/media/volume/st_andrews/databases/cATyper/hmm/msa_dec22/hmm/corA_cora.hmm",
        rows = base_path + "/08_CorA_finder/{j}/{j}_temp_rows.out",
        headers = base_path + "/08_CorA_finder/{j}/{j}_headers.out",
        all_data = base_path + "/08_CorA_finder/{j}/{j}_all_data.out",
    log:
        out = base_path + "/logs/08_CorA_finder/{j}.log",
        err = base_path + "/logs/08_CorA_finder/{j}.err"
    conda: "envs/hmmer.yaml"
    shell:
        '''
        if [ -s "{input.proteins}" ]; then
            echo "Proteins found for {wildcards.j}. Now starting hmmer" > {log.out}
            hmmscan --domtblout {params.out} --cpu 8 -E {corA_hmm_value} {params.corA_hmm_db} {input} &> /dev/null
            echo "Hmmscan finished for {wildcards.j}" >> {log.out}
            echo "Looking into {params.out}" >> {log.out}
            if [ -s {params.out} ]; then
                grep -v "#" {params.out} > {params.rows}||:
            else
                touch {params.rows}
                echo "Touching empty result file in {wildcards.j}" >> {log.out}
            fi
            echo "Reverse grepping for {wildcards.j} finished" >> {log.out}
            touch {output.hmm}
            echo "Adding headers for {wildcards.j}" >> {log.out}
            echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {params.headers}
            echo "Merging results with headers for {wildcards.j}" >> {log.out}
            cat {params.headers} > {output.hmm}
            cat {params.rows} >> {output.hmm}
            echo "Done for {wildcards.j}" >> {log.out}
        else
            echo "No proteins for {wildcards.j}" >> {log.out}
            echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm}
        fi
        '''


rule typeIII_new_effectors:
    '''
    Modified from SpacerCrawler's typeIII_new_effectors rule. This version
    focuses on CorA.

    NOTE: This doesn't output protein sequences for multiple type III loci-associated
            Cas10s or CorAs. Only one sequence of both will be output.
            Having multiple type III's is rare, but in may happen. Fix this in the future.

    From SpacerCrawler rule:
    This rule uses a simple algorithm to search for novel effectors candidates in type III CRISPR-Cas loci:
    1. Run cctyper and cATyper
    2. Get coordinates of found CRISPR-Cas III loci from cctyper (if any)
	3. Using cATyper output, check if these coordinates contain known effectors.
    4. Any locus that does not contain a known effector is now "suspicious"
	4. For each such suspicious locus, extract the cas operons png file from cctyper to a folder. Also copy the gbff files to one place for easier download. Then:
		a. Go through the gff file and flag any hypothetical proteins within the cas operon coordinates. Also get their positions in the cas operon using cctyper output
		b. Use this positional information to see how cctyper had annotated it (cas_operons.tab columns Genes and Positions and E-values)

        Final table:

    	protein     cctyper_suspicious	cc_typer_evalue	    cctyper_annotation  genome_annotation	            seq	        locus	locus_pos	host
		GCF_0001	yes	                1e-09	            CasR	            hypothetical protein	        MPAIJEWOIF			
		GCF_0002	yes	                1e-10	            asdas	            saved-domain containing protein	MADFOUSHDF		

    A concatenizing rule might take all sequences from all samples and cluster them. This way we could see how widespread a given new effector is.	

    '''
    input:
        cctyper = rules.CRISPRCasTyper.output, #note that this is just the dummy .done file. If no CRISPR-Cas systems were found, the actual files will not exist
        CorAFinder = rules.CorA_finder.output.hmm,
        gff = base_path + "/06_host_genomes/{j}/{j}_features.gff",
        proteins = base_path + "/06_host_genomes/{j}/{j}_proteins.faa",
    output:
        Cas10_fasta = base_path + "/09_crispr_iii_CorA/{j}/{j}_Cas10.faa",
        Cas5_fasta = base_path + "/09_crispr_iii_CorA/{j}/{j}_Cas5.faa",
        Cas7_fasta = base_path + "/09_crispr_iii_CorA/{j}/{j}_Cas7.faa",
        CorA_fasta = base_path + "/09_crispr_iii_CorA/{j}/{j}_CorA.faa",
        CorA_fasta_cctyper = base_path + "/09_crispr_iii_CorA/{j}/{j}_CorA_cctyper.faa",
        cora_type_iii_info = base_path + "/09_crispr_iii_CorA/{j}/{j}_crispr_iii_CorA.tsv",
        host_cas10_cora_info = base_path + "/09_crispr_iii_CorA/{j}/{j}_crispr_iii_info.tsv"
    params:
        cctyper_folder = base_path + "/07_cctyper/{j}",
        this_folder = base_path + "/09_crispr_iii_CorA/{j}"
    conda: "envs/gff_utils.yaml"
    log:
        out = base_path + "/09_crispr_iii_CorA/logs/{j}.out",
        err = base_path + "/09_crispr_iii_CorA/logs/{j}.err"
    shell:
        '''
        python3 scripts/typeIII_effector_finder_corA_simplified.py --cctyper_folder {params.cctyper_folder} --sample {wildcards.j} --CorAFinder {input.CorAFinder} --gff {input.gff} --out_cas10 {output.Cas10_fasta} --out_cas5 {output.Cas5_fasta} --out_cas7 {output.Cas7_fasta} --out_cora {output.CorA_fasta} --out_cora_cctyper {output.CorA_fasta_cctyper} --out_cora_info {output.cora_type_iii_info} --out_info {output.host_cas10_cora_info} --this_folder {params.this_folder} --proteins {input.proteins} 2> {log.err} 1> {log.out}
        touch {output}
        '''

rule concatenate_type_iii_info:
    input: aggregate_typeIII_info
    output: base_path + "/09_crispr_iii_CorA/type_iii_info.tsv"
    shell:
        '''
        echo "Sample\tCorA\tCas10\tGGDD\tGGDD_seq\tSubtype" > {output}
        find '{base_path}/09_crispr_iii_CorA' -maxdepth 2 -type f -wholename '*/*_crispr_iii_info.tsv' -print0 | xargs -0 cat >> {output}
        '''


rule getCas10_without_CorA:
    '''
    Gathers all type III CRISPR-Cas genomes without CorA,
    and copies their CCTyper plots and info tables to CorA_negatives folder
    '''
    input: rules.concatenate_type_iii_info.output
    output: base_path + "/09_crispr_iii_CorA/CorA_negatives/III_no_CorA.txt"
    params:
        plots = base_path + "/09_crispr_iii_CorA/CorA_negatives/graphs",
        CRISPR_tables = base_path + "/09_crispr_iii_CorA/CorA_negatives/tables",
        cctyper = base_path + "/07_cctyper"
    shell:
        """
        awk '($2=="True") && ($3=="True") {{print $1}}' {input} > {output}
        mkdir {params.plots}
        mkdir {params.CRISPR_tables}
        while read sample; do
            PLOT={params.cctyper}/$sample/plot.png
            INFO={params.cctyper}/$sample/CRISPR_Cas.tab
            if [ -s $PLOT ]; then
                cp $PLOT {params.plots}/$sample.png
                cp $INFO {params.CRISPR_tables}/$sample.tab
            fi
        done < {output}
        """

rule CorA_concatenate:
    input: aggregate_CorA_sequences
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
        if [ {protein_clustering} = "True" ]; then
            cd-hit -i {input} -o {output.proteins} -c 0.99 -n 5 -d 0 -M 16000 -T {threads}
        else
            cp {input} {output.proteins}
            touch {output.clusterinfo}
        fi
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


rule Cas10_concatenate:
    input: aggregate_cas10_sequences
    output: base_path + "/09_crispr_iii_CorA/cas10s.faa"
    shell:
        '''
        cat {input} > {output} 
        '''

rule Cas10_cluster:
    '''
    Sketching:
    Apparently it would be best to cluster cas10 proteins and then use representative sequences for the tree.
    These representatives would then carry alongside with them information on their respective host and 16S
    associations.
    '''
    input: 
        cas10 = rules.Cas10_concatenate.output,
        corA = rules.CorA_tree.output #just for chaining rules, not actually used as input
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

rule Cas5_concatenate:
    input: aggregate_cas5_sequences
    output: base_path + "/09_crispr_iii_CorA/cas5s.faa"
    shell:
        '''
        cat {input} > {output} 
        '''

rule Cas5_cluster:
    '''
    Sketching:
    Apparently it would be best to cluster cas5 proteins and then use representative sequences for the tree.
    These representatives would then carry alongside with them information on their respective host and 16S
    associations.
    '''
    input: 
        cas5 = rules.Cas5_concatenate.output,
        corA = rules.CorA_tree.output #just for chaining rules, not actually used as input
    output:
        proteins = base_path + "/22_Cas5_cluster/cas5_cluster.faa",
        clusterinfo = base_path + "/22_Cas5_cluster/cas5_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/22_Cas5_cluster/logs/cas5_align.out",
        err = base_path + "/22_Cas5_cluster/logs/cas5_align.err"
    threads: 40
    shell:
        '''
        echo {protein_clustering}
        if [ {protein_clustering} = "True" ]; then
            cd-hit -i {input.cas5} -o {output.proteins} -c 0.99 -n 5 -d 0 -M 16000 -T {threads}
        else
            cp {input.cas5} {output.proteins}
            touch {output.clusterinfo}
        fi
        '''


rule Cas5_align:
    input: rules.Cas5_cluster.output.proteins
    output: base_path + "/23_Cas5_align/cas5_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/23_Cas5_align/logs/cas5_align.out",
        err = base_path + "/23_Cas5_align/logs/cas5_align.err"
    threads: 40
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''



rule Cas5_tree:
    input: rules.Cas5_align.output
    output: base_path + "/24_Cas5_tree/cas5_tree.txt",
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''


rule Cas7_concatenate:
    input: aggregate_cas7_sequences
    output: base_path + "/09_crispr_iii_CorA/cas7s.faa"
    shell:
        '''
        cat {input} > {output} 
        '''

rule Cas7_cluster:
    '''
    Sketching:
    Apparently it would be best to cluster cas7 proteins and then use representative sequences for the tree.
    These representatives would then carry alongside with them information on their respective host and 16S
    associations.
    '''
    input: 
        cas7 = rules.Cas7_concatenate.output,
        corA = rules.CorA_tree.output #just for chaining rules, not actually used as input
    output:
        proteins = base_path + "/25_Cas7_cluster/cas7_cluster.faa",
        clusterinfo = base_path + "/25_Cas7_cluster/cas7_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/25_Cas7_cluster/logs/cas7_align.out",
        err = base_path + "/25_Cas7_cluster/logs/cas7_align.err"
    threads: 40
    shell:
        '''
        echo {protein_clustering}
        if [ {protein_clustering} = "True" ]; then
            cd-hit -i {input.cas7} -o {output.proteins} -c 0.99 -n 5 -d 0 -M 16000 -T {threads}
        else
            cp {input.cas7} {output.proteins}
            touch {output.clusterinfo}
        fi
        '''


rule Cas7_align:
    input: rules.Cas7_cluster.output.proteins
    output: base_path + "/26_Cas7_align/cas7_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/26_Cas7_align/logs/cas7_align.out",
        err = base_path + "/26_Cas7_align/logs/cas7_align.err"
    threads: 40
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule Cas7_tree:
    input: rules.Cas7_align.output
    output: base_path + "/27_Cas7_tree/cas7_tree.txt",
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''



rule align_16s:
    input: rules.concatenate_16s.output,
    output: base_path + "/20_16s_align/16s_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/20_16s_align/logs/cas10_align.out",
        err = base_path + "/20_16s_align/logs/cas10_align.err"
    threads: 40
    shell:
        '''
        mafft {input} > {output}
        '''


rule tree_16s:
    '''

    '''
    input: rules.align_16s.output
    output: base_path + "/21_16s_tree/16s_tree.txt",
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -nt -gtr -gamma {input} > {output}
        '''


checkpoint align_matcher:
    '''
    This rule checks all possible pairwise comparisons and reduces the alignments
    so that a tree will only be created out of those species that
    have both proteins in their genomes. This enables the creation
    of trees with the same number of nodes for each comparison.

    If thi
    '''
    input:
        cas10_alignment = rules.Cas10_align.output,
        cas7_alignment = rules.Cas7_align.output,
        cas5_alignment = rules.Cas5_align.output,
        cora_alignment = rules.CorA_align.output,
        rna16S_alignment = rules.align_16s.output,
    output: directory(base_path + "/40_alignmatcher/")
        # cas10_cora_cas10 = base_path + "/40_alignmatcher/cas10_cora/cas10.afa",
        # cas10_cora_cora = base_path + "/40_alignmatcher/cas10_cora/cora.afa",
        # cas10_cas7_cas10 = base_path + "/40_alignmatcher/cas10_cas7/cas10.afa",
        # cas10_cas7_cas7 = base_path + "/40_alignmatcher/cas10_cas7/cas7.afa",
        # cas10_cas5_cas10 = base_path + "/40_alignmatcher/cas10_cas5/cas10.afa",
        # cas10_cas5_cas5 = base_path + "/40_alignmatcher/cas10_cas5/cas5.afa",
        # cas10_16s = base_path + "/40_alignmatcher/cas10_16s/cas10.afa",
        # cas10_16s = base_path + "/40_alignmatcher/cas10_16s/16s.afa",
        # cas7_cas5_cas7 = base_path + "/40_alignmatcher/cas7_cas5/cas7.afa",
        # cas7_cas5_cas5 = base_path + "/40_alignmatcher/cas7_cas5/cas5.afa",
        # cas7_16s_cas7 = base_path + "/40_alignmatcher/cas7_16s/cas7.afa",
        # cas7_16s_16s = base_path + "/40_alignmatcher/cas7_16s/16s.afa",
        # cas7_cora_cas7 = base_path + "/40_alignmatcher/cas7_cora/cas7.afa",
        # cas7_cora_cora = base_path + "/40_alignmatcher/cas7_cora/cora.afa",
        # cas5_16s_cas5 = base_path + "/40_alignmatcher/cas5_16s/cas5.afa",
        # cas5_16s_16s = base_path + "/40_alignmatcher/cas5_16s/16s.afa",
        # cas5_cora_cas5 = base_path + "/40_alignmatcher/cas5_cora/cas5.afa",
        # cas5_cora_cora = base_path + "/40_alignmatcher/cas5_cora/cora.afa",
        # cora_16s_cora = base_path + "/40_alignmatcher/cora_16s/cora.afa",
        # cora_16s_16s = base_path + "/40_alignmatcher/cora_seq_16s/16s.afa",
    params:
        outfolder = base_path + "/40_alignmatcher"
    shell:
        '''
        mkdir {params.outfolder}
        python scripts/align_matcher.py --cas10 {input.cas10_alignment} --cas7 {input.cas7_alignment} --cas5 {input.cas5_alignment} --cora {input.cora_alignment} --rna16s {input.rna16S_alignment} --outfolder {params.outfolder}
        '''

rule tree_reduced_alignments:
    '''
    Runs FastTree for all alignments in the pairwise folder (two files).
    The wildcard is now the combination of two proteins (e.g. k = cas10_cora).
    The wildcard is formulated automatically in the python script of the rule align_matcher.
    Different runs of the pipeline may produce different naming schemes regarding
    which of the two components comes first (cas10_cora or cora_cas10), as they are constructed
    from a non-structured dictionary.
    '''
    input: base_path + "/40_alignmatcher/{k}.done"
        #first = base_path + "/40_alignmatcher/{i}/done",
        #second = base_path + "/40_alignmatcher/{i}/",
    output:
        done = base_path + "/41_matching_trees/{k}/done.txt"
    params:
        in_folder = base_path + "/40_alignmatcher/{k}",
        out_folder = base_path + "/41_matching_trees/{k}",
    conda: "envs/trees.yaml"
    shell:
        '''
        for input_file in {params.in_folder}/*.afa; do
            filename=$(basename "$input_file" .afa)
            FastTree -wag -gamma "$input_file" > "{params.out_folder}/$filename.tree"
        done
        touch {output}
        '''

rule matched_distance_matrices_comparison:
    '''
    Extracts distance matrices from the pairwise alignments (not trees!)
    '''
    input:
        #tree = rules.tree_reduced_alignments.output,
        alignments = rules.tree_reduced_alignments.output.done
        #first = base_path + "/40_alignmatcher/{i}/done",
        #second = base_path + "/40_alignmatcher/{i}/",
    output:
        #done = base_path + "/42_matching_trees_distance_matrices/{k}/done.txt,
        correlation = base_path + "/42_matching_trees_distance_matrices/{k}/{k}_correlation.tsv"
    params:
        #in_folder = base_path + "/41_matching_trees/{k}",
        in_folder = base_path + "/40_alignmatcher/{k}",
        out_folder = base_path + "/42_matching_trees_distance_matrices/{k}"
    run:
        from Bio import Phylo, AlignIO
        print(dir(Phylo))
        from Bio.Phylo import TreeConstruction
        from Bio.Phylo.TreeConstruction import DistanceCalculator
        import os, glob
        import numpy as np
        from scipy.stats import pearsonr
        import matplotlib.pyplot as plt
        from scipy.spatial.distance import squareform

        trees = {} #this dictionary contains two values: a dictionary for the tree and a dic for the distance matric
        names = []

        path = params.in_folder + '/*.tree'
        print(path)

        calculator_prot = DistanceCalculator('blosum62') #used to create prot dms
        calculator_nt = DistanceCalculator('identity') #used to create nt dms


        for filename in glob.glob(params.in_folder + '/*.afa'):
            print(filename) #full path
            simplename = str(os.path.basename(str(filename))).split(".")[0] #filename without extension, e.g. cas10_cas7_cas10
            print(simplename)
            protein = simplename.split("_")[2]
            names.append(simplename)
            with open(os.path.join(os.getcwd(), filename), 'r') as f: # open in readonly mode
                print("Extracting DM from " + simplename)
                trees[simplename] = {}
                trees[simplename]["alignment"] = AlignIO.read(filename, "fasta")
                if protein == "16s":
                    trees[simplename]["dm"] = calculator_nt.get_distance(trees[simplename]["alignment"])
                else:
                    trees[simplename]["dm"] = calculator_prot.get_distance(trees[simplename]["alignment"])
                #trees[simplename]["dm"] = DistanceMatrix.from_tree(trees[simplename]["tree"])
                #np.savetxt(params.out_folder + "/" + simplename + ".dm", trees[simplename]["dm"])
            
        name_0 = names[0].split("_")[2]
        name_1 = names[1].split("_")[2]
        print(name_0)
        print(name_1)
        dm0 = trees[names[0]]["dm"]
        dm1 = trees[names[1]]["dm"]
        corr_coef, p_value = pearsonr(squareform(dm0), squareform(dm1))

        print("Correlation coefficient for " + wildcards.k + ": " + str(corr_coef))
        print("P-value for " + wildcards.k + ": " + str(p_value))

        with open(params.out_folder + "/" + wildcards.k + "_correlation.tsv", 'w') as f:
            f.write(str(wildcards.k + "\t" + str(p_value) + "\t" + str(corr_coef) + "\n"))

        with open(params.out_folder + "/" + wildcards.k + "_dm1.tsv", 'w') as f:
            f.write(str(trees[names[0]]["dm"]))

        with open(params.out_folder + "/" + wildcards.k + "_dm2.tsv", 'w') as f:
            f.write(str(trees[names[1]]["dm"]))

        plt.figure(figsize=(8, 8))
        plt.scatter(dm0, dm1, s=75, alpha = 0.2)
        plt.title("DM (" + name_0 + "/" + name_1 + ")", fontsize=18, y=1.03)
        plt.text(1.2, 0.1, 'r = ' + str(round(float(corr_coef),2)), fontsize = 16)
        plt.xlabel(name_0)
        plt.ylabel(name_1)
        plt.xlim([0, 1.5])
        plt.ylim([0, 1.5])
        plt.savefig(params.out_folder + "/" + wildcards.k + "_dm_plot.png")
        
rule gather_dms:
    input: aggregate_dms
    output: base_path + "/42_matching_trees_distance_matrices/dm_stats.txt"
    shell:
        '''
        echo "Pair\tP-value\tPearson_correlation" > {output}
        cat {input} >> {output}
        mkdir {base_path}/42_matching_trees_distance_matrices/plots
        cp {base_path}/42_matching_trees_distance_matrices/*/*.png {base_path}/42_matching_trees_distance_matrices/plots
        '''


rule clusterLinker:
    '''
    EXPERIMENTAL (MAYBE NOT NEEDED)
    We might need to link original wildcards to each cluster for meaningful results.
    This rule produces a table that, for each protein,
    shows the corresponding cluster representative of that cluster.
    '''
    input:
        cas10 = rules.Cas10_cluster.output.clusterinfo,
        cas5 = rules.Cas5_cluster.output.clusterinfo,
        cas7 = rules.Cas7_cluster.output.clusterinfo,
        cora = rules.CorA_cluster.output.clusterinfo,
        #rna_16s = rules.cluster.output.clusterinfo
    output:
        cluster_table = base_path + "/30_clustertable/{j}/{j}.tsv"
    params:
        revCas10 = base_path + "/30_clustertable/{j}/{j}_rev.clstr"
    shell:
        '''
        if [ {protein_clustering} = "True" ]; then
            python3 scripts/clusterLinker.py --cas10_clusters {input.cas10} --cas7_clusters {input.cas7} --cas5_clusters {input.cas5} --cora_clusters {input.cora} --out {output.cluster_table} --sample {wildcards.j}
        else
            touch {output}
        fi
        '''

rule concatenate_clustertables:
    input: aggregate_clustertables
    output: base_path + "/clusters.tsv"
    shell:
        '''
        echo "Sample\tCas10_link\tCas7_link\tCas5_link\tCorA_link\t16S_link" > {output}
        cat {input} >> {output}
        '''

rule final:
    input:
        tree_16s = rules.tree_16s.output,
        tree_Cas10 = rules.Cas10_tree.output,
        tree_CorA = rules.CorA_tree.output,
        tree_Cas5 = rules.Cas5_tree.output,
        tree_Cas7 = rules.Cas7_tree.output,
        #clustertables = rules.concatenate_clustertables.output,
        matched_distance_matrices_comparison = rules.gather_dms.output,
        type_iii_info = rules.concatenate_type_iii_info.output,
        getCas10_without_CorA = rules.getCas10_without_CorA.output,
        concat_taxInfo = rules.concat_taxInfo.output
    output: base_path + "/done"
    shell:
        '''
        touch {output}
        '''


# rule merge:
#     '''
#     We need a final table like this to annotate the trees:

#     id              |   Cas10  |   CorA    |   cas10_dist  |   cora_dist   |   cas5_dist   |  
#     GCF_987234.1    |    true  |   true    |       2.1     |       0.2     |       1.2     |
#     GCF_854732.1    |    false |   true    |       NaN     |       2.2     |       4.2     |

#     '''
#     input:
#         tree_16s = rules.tree_16s.output,
#         tree_Cas10 = rules.Cas10_tree.output,
#         tree_CorA = rules.CorA_tree.output,
#     output:
#     shell:
#         '''

#         '''
# rule 16SBlast:
#     input:
#         db = rules.makeBlastDBnt.output.blastdb,
#     output: base_path + "/061_host_genomes_blastDB/genomes/blast.out"
#     threads: 40
#     shell:
#         '''
#         blastn -num_threads {threads} -db {input.db} -task blastn -query {output.fasta} -outfmt "6 qseqid qlen sseqid  stitle staxids saccver sstart send evalue length pident mismatch gapopen gaps sstrand" -evalue 0.01 > {output.blast}
#         '''


# rule align_cas10:
#     '''
#     Requires as input extracted Cas10 sequences from each sample i.
#     Metainfo 
#     '''
#     input:
#         cas10s = concatenate_cas10s
#     output: base_path + "/09_crispr_iii_CorA/{i}/{i}_crispr_iii_CorA.tsv"
#     params:
#         cctyper_folder = base_path + "/07_cctyper/{i}",
#         this_folder = base_path + "/09_crispr_iii_CorA/{i}"
#     conda: "envs/gff_utils.yaml"
#     log:
#         out = base_path + "/09_crispr_iii_CorA/logs/{i}.out",
#         err = base_path + "/09_crispr_iii_CorA/logs/{i}.err"
#     shell:
#         '''

#         '''


#     run:
#         import pandas as pd
#         if not os.path.exists(str(output)):
#             os.makedirs(str(output))
#         cas10files = pd.read_csv(str(input), sep = "\t", header = 0)
#         phages = phagefile["phage"]
#         for i in phages:
#             with open(str(output) + "/" + str(i) + '.txt', 'w') as f:
#                 f.write(str(output))



 #   shell:
 #       '''
 #       cd {base_path}/04_host_genome_files/
 #       datasets download genome accession $(cat {input}) --filename {output.zip} --include genome,protein,gff3 --assembly-level complete
 #       '''
#        ncbi-genome-download --accessions {wildcards.i} --formats fasta,protein-fasta,gff --assembly-levels complete  bacteria,archaea


