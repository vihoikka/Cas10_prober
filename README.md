# Cas10_prober
Snakemake pipeline for characterizing type III CRISPR-Cas loci and related CorAs. Generates phylogenetic trees with annotations. Chi et al. 2023, manuscript submitted.

## Requirements
- You need to have [Snakemake]([url](https://anaconda.org/bioconda/snakemake)), [Conda]([url](https://docs.conda.io/en/latest/index.html)) and [Hmmer]([url](https://anaconda.org/bioconda/hmmer)) installed. Other dependencies are installed automatically by the pipeline.
- Runs on Unix environments. Tested only on Ubuntu
- Data: you need a local database of genomes. You can use [NCBI Datasets command line tool].(https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=datasets-command-line-20221012) to download the genomes. Use the command ```datasets download genome taxon 2 --annotated --assembly-level complete --include genome,protein,gff3```

## How to run
1. Clone the repository
2. HMM profiles and related protein alignments are provided in the msa_050523 -folder. Use Hmmer to hmmpress the HMM database (```hmmpress 050523.hmm```) and modify paths in the pipeline to point to the database. Also modify the hmm_msa_folder variable in the script to point to the directory with the alignments.
3. Run using the following command

```
snakemake --snakefile cas10_prober.smk --use-conda --cores 40 --config protein_clustering="False" getGenomesBy="local" genome_mode="all" cas10_anchor="True" --rerun-triggers mtime
```

