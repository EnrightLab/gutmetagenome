# Gut Metagenomics Supporting Material

# Dalia Bornstein Part II Project

This is a Github supporting Dalia's Part II Project which contains computational methods and code snippets described in my part II Writeup.


### Illumina 16S Amplicon Analysis

### Nanopore Basecalling Command

|Flow cell type | Flow cell ID | Kit Used |
|---------------|--------------|----------|
|FLO-PRO114M    | PAU22393 | SQK-NBD114-24|

`High-accuracy model v4.3.0, 400 bps`

`GCF_000001405.40_GRCh38.p14_genomic.fna`

|Software | Version |
|---------|---------|
|MinKNOW | 24.06.14 |
|Bream | 8.0.12 | 
|Configuration | 6.0.15 |
|Dorado  | 7.4.13 | 
|MinKNOW | Core 6.0.11 |

Kit type


This command was run and it was version `guppy 4.0.32x 12/10/24`. `guppy-ont-gpu -in filex.fastq -out basecalled.fastq.gz -evalue 10`

`guppy-ont-gpu -in filex.fastq -out basecalled.fastq.gz -evalue 10`

### Metagenomic Mapping Tools Used

-   Kraken version 2.1.3
-   Bracken v2.9

#### Mapping Options

`kraken2 --use-names --threads 4 --confidence 0.5 --db kraken_db --report barcodeX.krakenreport.txt --gzip-compressed barcodeX.fastq.gz > barcodeX.kraken2.txt`

* --use-names (Print scientific names instead of just taxids)
* --threads 4 (Use 4 cpus per job)
* --confidence 0.5 (Confidence score threshold)
* --db kraken_db (Use the database, details below)
* --report (Save text report to this file of read mappings)
* --gzip-compressed (Read query reads from a gzip file, e.g. fastq.gz)

#### Bracken Abundance Analysis
`est_abundance.py -i barcodeX.kraken2.txt -k kraken_db/database75mers.kmer_distrib -l S -t 10 -o barcodeX.bracken.txt`

* -i barcodeX.kraken2.txt (input file from kraken)
* -k kraken_db/database75mers.kmer_distrib (k-mers file from database)
* -l S (species level calls)
* -t 10 (use 10 cpus)

#### Database Used

The metagenomic reference used was the *core_nt* from NCBI. This very large collection, includes of GenBank, RefSeq, TPA and PDB.
This database had a *k-mer* size of 35.

* Database Date: 12/28/2024 
* Database Details [Link](https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20241228.tar.gz) 

#### Database options

```
# Database options: nucleotide db, k = 35, l = 31
# Spaced mask = 11111111111111111111111111111111110011001100110011001100110011
# Toggle mask = 1110001101111110001010001100010000100111000110110101101000101101
# Total taxonomy nodes: 2085417
# Table size: 43826298697
# Table capacity: 62617751405
# Min clear hash value = 0
```
