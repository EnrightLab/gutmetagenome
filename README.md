# Gut Metagenomics Supporting Material

# Dalia Bornstein Part II Project

This is a Github supporting Dalia's Part II Project which contains computational methods and code snippets described in my part II Writeup.

### Nanopore Basecalling Command

This command was run and it was version `guppy 4.0.32x 12/10/24`. `guppy-ont-gpu -in filex.fastq -out basecalled.fastq.gz -evalue 10`

`guppy-ont-gpu -in filex.fastq -out basecalled.fastq.gz -evalue 10`

### Anton Adds something here

`some command`

### Metagenomic Mapping Used

-   Kraken version 2.1.3
-   Bracken v2.9

`kraken2 --use-names --threads 4 --confidence 0.5 --db kraken_db --report barcodeX.krakenreport.txt --gzip-compressed barcodeX.fastq.gz > barcodeX.kraken2.txt`

#### Database Used

The metagenomic reference used was the *core_nt* from NCBI. This very large collection, includes of GenBank, RefSeq, TPA and PDB.
This database had a *k-mer* size of 35.

* Database Date: 12/28/2024 
* Database Details [Link](https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20241228.tar.gz) 

#### Database options
`
# Database options: nucleotide db, k = 35, l = 31
# Spaced mask = 11111111111111111111111111111111110011001100110011001100110011
# Toggle mask = 1110001101111110001010001100010000100111000110110101101000101101
# Total taxonomy nodes: 2085417
# Table size: 43826298697
# Table capacity: 62617751405
# Min clear hash value = 0
`
