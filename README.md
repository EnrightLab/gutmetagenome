# Dalia Bornstein - Part II Project

This is a Github supporting Dalia Bornstein's Part II Project which contains computational methods and code as described in my part II report.

---

## Illumina 16S Amplicon Analysis
This section describes the analysis framework (mostly in R/Bioconductor) for processing and analysing
the 16S amplicon sequencing data. The raw data is from Illumina MiSeq and the runs were Paired-end 250nt.

* NEXTFLEX® 16S V1 - V3 rRNA Amplicon-Seq Library Prep Kit
This [kit](https://www.revvity.com/gb-en/product/nextflex-16s-v1-v3-12-24-rxn-nova-4202-02s#product-overview) is designed for the preparation of multiplexed amplicon libraries that span the hypervariable domains one through three (V1-V3) of microbial 16S ribosomal RNA (rRNA) genes. These libraries are compatible with paired-end sequencing on the Illumina® sequencing platforms.
* Sequencer Used Illumina MiSeq at Cambridge Genomic Services (Pathology).

### 16S Database Used

We used the [SILVA database](https://benjjneb.github.io/dada2/training.html).

* Training Database - SILVA non-redundant version v132 `silva_nr_v132_train_set.fa`
* Species Database - SILVA species assignment v132 `silva_species_assignment_v132.fa`

### DADA2 Filtering Options 

* `truncLen=c(221,141)` - Filter forward reads by removing 29nt and 110nt from the reverse
* `maxN=0` - Maximum number of 'N' bases before filtering
* `truncQ=2` - Truncate reads with Quality threshold 2
* `rm.phix=TRUE` - Remove reads which are PhiX contaminants
* `compress=TRUE` - Compress the output (gz)
* `multithread=TRUE` - Use maximum available CPUs

### R Packages Used

|package | version | type         |
|--------|---------|--------------|
|dada2   | 1.34.0  | bioconductor |
|phyloseq| 1.50.0  | bioconductor |
|microbiome| 1.28.0 | bioconductor |
|microViz | 0.12.6 | bioconductor |
|RColorBrewer|1.1.3|bioconductor |
|kableExtra|1.4.0|Cran|
|Biostrings|2.74.1|bioconductor|
|ggplot2|3.5.1|Cran|
|dplyr|1.1.4|Cran|
|gdata|3.0.1|Cran|

### Illumina 16S Analysis R/BioConductor Code

The full codebase in R/BioConductor and source data are in the [markdown provided](Illumina-16S-DADA2.md).

* [Illumina-16S-DADA2.md](Illumina-16S-DADA2.md) - Markdown For Github Viewing.
* [Illumina-16S-DADA2.Rmd](Illumina-16S-DADA2.Rmd) - RMarkdown master script.

---

## Nanopore Metagenomics

### Nanopore Details

* Sequencer - PromethION P24 - Cambridge Genomic Services (Pathology).
* Run Time  - 120 hours (both runs)
* Run 1 - Adaptive Sampling - Deplete for Human Genes
* Run 2 - No Adaptive Sampling - for comparison

|Run  |Flow cell type | Flow cell ID | Kit Used |
|-----|---------------|--------------|----------|
|Run 1| FLO-PRO114M   | PAU22393  |SQK-NBD114-24|
|Run 2| FLO-PRO114M   | PAS76652  |SQK-NBD114-24|


### Dorado Basecalling Model Used
`High-accuracy model v4.3.0, 400 bps`

### Adaptive Sampling Human Reference Used
`GCF_000001405.40_GRCh38.p14_genomic.fna`

### Software versions Used
|Software | Version |
|---------|---------|
|MinKNOW | 24.06.14 |
|Bream | 8.0.12 | 
|Configuration | 6.0.15 |
|Dorado  | 7.4.13 | 
|MinKNOW | Core 6.0.11 |

### Metagenomic Mapping Tools Used

-   Kraken version 2.1.3
-   Bracken v2.9

### Metagenomic Mapping Options

These were the command-lines used for analysis where `X` represents the barcode being analysed.

```
kraken2 --use-names --threads 4 --confidence 0.5 --db kraken_db --report barcodeX.krakenreport.txt --gzip-compressed barcodeX.fastq.gz > barcodeX.kraken2.txt
```

* `--use-names` (Print scientific names instead of just taxids)
* `--threads 4` (Use 4 cpus per job)
* `--confidence 0.5` (Confidence score threshold)
* `--db kraken_db` (Use the database, details below)
* `--report` (Save text report to this file of read mappings)
* `--gzip-compressed` (Read query reads from a gzip file, e.g. fastq.gz)

### Bracken Abundance Analysis
```
est_abundance.py -i barcodeX.kraken2.txt -k kraken_db/database75mers.kmer_distrib -l S -t 10 -o barcodeX.bracken.txt
```

* `-i barcodeX.kraken2.txt` (input file from kraken)
* `-k kraken_db/database75mers.kmer_distrib` (k-mers file from database)
* `-l S` (species level calls)
* `-t 10` (use 10 cpus)
* `-o barcodeX.bracken.txt` (output abundance information in new report file)

### Metagenomic Database Used

The metagenomic reference used was the *core_nt* from NCBI. This very large collection, includes of GenBank, RefSeq, TPA and PDB.
This database had a *k-mer* size of 35.

* Database Date: 12/28/2024 
* Database Details [Link](https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20241228.tar.gz) 

### Database options
```
# Database options: nucleotide db, k = 35, l = 31
# Spaced mask = 11111111111111111111111111111111110011001100110011001100110011
# Toggle mask = 1110001101111110001010001100010000100111000110110101101000101101
# Total taxonomy nodes: 2085417
# Table size: 43826298697
# Table capacity: 62617751405
# Min clear hash value = 0
```

---

## Metagenome Visualisation

### Krona Plots
Individual Bracken report files (or Kraken report files) can be turned into a Krona HTML visualisation.
For this we use the following command from `KronaTools v2.8.1`:

```
ktImportTaxonomy -t 5 -m 3 -o barcodeX.html barcodeX.bracken.report.txt
```

This uses 5 threads and generates a html report (barcodeX.html) from the input report.

### Circos Plots

* Circos Version Used: `circos 0.69-9`

For genomes of interest we download a reference in *Genbank* format, e.g. `Akkermansia.gb` is the Genbank entry (CP001071.1) for *Akkermansia sp*.
We then use a perl script [find_species_hits.pl](scripts/find_species_hits.pl) to interrogate all kraken mapping files to identify reads hitting the species of interest.
These reads are extracted into a new *FASTA* file. We also extract a genome *FASTA* file from the same genbank file.

To identify per genome mappings we use *blastn* to map the reads against the reference genome *FASTA*.

* `blast 2.16.0, build Jun 25 2024 12:36:54`

We build a BLAST database as follows:
`makeblastdb -in genome.fasta -dbtype="nucl"`

blastn is performed as follows:
```
blastn -query query.fasta -db genome.fasta -outfmt '7 qseqid sseqid pident mismatch gapopen qstart qend sstart send evalue bitscore length qlen qseq sseq'  -word_size 8 -evalue 1e-5 -num_threads 12
```

This is performed by an accessory script [blaster.pl](scripts/blaster.pl) which handles the *blastn* and filtering to one hit per metagenomic read.

This produces tabular output which is processed in perl to annotate hits like below:
```
b9e3f48f-dc8e-41de-af98-590048dbcb13    No Hits NA
13264101-5fe5-4c8f-8d91-8b8667594101    CP001071.1      95.039  96.90% [1532 / 1581]    2515534 2517058       74      1580
932f7c9c-e167-4c25-a687-8514ad674880    CP001071.1      96.647  71.63% [3788 / 5288]    2058879 2055109       75      3830
7eec27c8-4c8a-4135-84b3-5a6d48f4eb59    CP001071.1      97.650  91.63% [766 / 836]      94455   93691         77      836
9c5f7fe5-658a-4954-ac5f-a14e57f41191    CP001071.1      88.475  87.93% [590 / 671]      972071  971489        75      631
7545b3ec-6706-4c8e-8835-c12f2023b137    CP001071.1      97.102  96.10% [2588 / 2693]    188464  191038        69      2640
9b229828-2dc0-4c01-b21a-d7283041bae2    CP001071.1      96.458  72.35% [3783 / 5229]    2055109 2058879       1438    5181
8980fbd9-053b-49d6-b6fd-345f9a9b3023    CP001071.1      97.114  55.41% [2564 / 4627]    150877  148323        71      2613
459602a3-72e7-408d-93f2-65bb2bb80f92    CP001071.1      98.169  99.11% [5242 / 5289]    1803894 1798665       68      5284
a3c1161d-2d67-4a35-ac4f-66d18aa0dd87    CP001071.1      97.557  96.35% [2824 / 2931]    1834344 1831529       72      2881
1d47265a-ccd6-4e8f-bc16-b97bf067b45b    CP001071.1      99.571  77.15% [233 / 302]      74149   73917         71      302
```

### Genome Skew

The Python script gcskew.py (Jennifer Lu, jlu26\@jhmi.edu) is used to compute GC skew in 10kb bins from the *FASTA* file.
This is saved into a new GC skew file, e.g.:

```
Sequence        Index   GC Skew (20kb)
CP001071.1      0       -0.10094022
CP001071.1      20000   -0.07867751
CP001071.1      40000   -0.06978027
CP001071.1      60000   -0.08492865
CP001071.1      80000   -0.05014085
CP001071.1      100000  -0.09060134
CP001071.1      120000  -0.04795584
CP001071.1      140000  -0.08262656
CP001071.1      160000  -0.09920739
CP001071.1      180000  -0.06964381
```

### Genome Annotation and Coverage

A custom perl script interrogates selected genomes, processes the *BLAST* output and explores the genbank annotation.
* BLAST hits are summarised as coverage in bins (usually over 1kb) by taking an average count of all reads spanning each nt in that bin.
* Genbank annotation is stored with its location as a highlights track to display genes
* Genbank annotation is stored with locations as a labels track to highlight key genes
* The coverage is stored in a histogram track file normalised to the maximum coverage value as 1.0.
* The GC skew file (if computed) is also stored in a histogram track file for display.

These track files are:

| filename                | circos type | description.         |
|-------------------------|-------------|----------------------|
|`labels_forward.txt`     | labels | forward strand gene names |
|`labels_reverse.txt`     | labels | reverse strand gene names |
|`highlights_forward.txt` | highlight | forward strand genes   |
|`highlights_reverse.txt` | highlight | reverse strand genes   |
|`karyotype.txt`          | karyotype | visual karyotype       |
|`data_tile.txt`          | histogram | coverage of nanopore   |
|`data_skew.txt`          | histogram | GC skew                |

#### Command Line
This perl script [make_karyotype_gb.pl](scripts/make_karyotype_gb.pl) computes the key track files:

```
./make_karyotype_gb.pl akkermansia.gb akkermansia.hits akkermansia.skew --names; circos -conf metagenome.conf
```
This script makes use of [genbank_gtf.pl](scripts/genbank_gtf.pl) by Jiang Li (Vanderbilt Center for Quantitative Sciences).


This creates a *png* and an *svg* vector graphic for each species of interest.

All species are run with a simple shell script:
```
#!/bin/sh
./make_karyotype_gb.pl akkermansia.gb akkermansia.hits akkermansia.skew --names; circos -conf circos.conf; cp circos.svg akkermansia.svg; open circos.png
./make_karyotype_gb.pl herelleviridae.gb herelleviridae.hits --names;  circos -conf circos.conf; cp circos.svg herelleviridae.svg; open circos.png
./make_karyotype_gb.pl methano_smithii.gb methanobrevibacter.hits methanobrevibacter.skew; circos -conf circos.conf; cp circos.svg methano_smithii.svg; open circos.png
./make_karyotype_gb.pl methanomassiliicoccus.gb methanomassiliicoccus.hits methanomassiliicoccus.skew --names; circos -conf circos.conf; cp circos.svg methanomassiliicoccus.svg; open circos.png
./make_karyotype_gb.pl methanosphaera.gb methanosphaera.hits methanosphaera.skew --names;  circos -conf circos.conf; cp circos.svg methanosphaera.svg; open circos.png
./make_karyotype_gb.pl suterella.gb suterella.hits suterella.skew;  circos -conf circos.conf; cp circos.svg suterella.svg; open circos.png
./make_karyotype_chr.pl candida.gb candida.hits; circos -conf circos.conf; cp circos.svg candida.svg; open circos.png
```


### Circos configuration file

The configuration file [circos.conf](scripts/circos.conf) to show an outer karyotype, inner gene bodies and labels, coverage and finally genome skew is as follows:

```
<<include etc/colors_fonts_patterns.conf>>
<<include ideogram.conf>>
<<include ticks.conf>>

karyotype = karyotype.txt

<image>
<<include etc/image.conf>>
</image>

chromosomes_units           = 1000
chromosomes_display_default = yes

<plots>
show = yes

## Show Gene Labels
<plot>
type  = text
file  = labels_forward.txt
color = black
r1    = 0.98r
r0    = 0.88r
label_size = 8
label_font = light
padding    = 4p
rpadding   = 4p
show_links     = yes
link_dims      = 5p,4p,8p,4p,0p
link_thickness = 1p
link_color     = dgrey
label_snuggle        = yes
</plot>

<plot>
type  = text
file  = labels_reverse.txt
color = black
r1    = 0.85r
r0    = 0.75r
label_size = 8
label_font = light
padding    = 4p
rpadding   = 4p
show_links     = yes
link_dims      = 5p,4p,8p,4p,0p
link_thickness = 1p
link_color     = dgrey
label_snuggle        = yes
</plot>

## Show Coverage Histogram
<plot>
type        = histogram
file        = data_tile.txt
r1          = 0.68r
r0          = 0.58r
min         = 0
max         = 1
extend_bin  = no
fill_color  = dblue
color       = blue
thickness   = 0
orientation = out
<axes>
<axis>
color     = grey_a1
thickness = 1
spacing   = 0.10r
</axis>
</axes>
</plot> 

## Show GC Skew Histogram
<plot>
type        = histogram
file        = data_skew.txt
color       = black
thickness   = 2
r1          = 0.58r
r0          = 0.48r
max         = 0.25 
min         =-0.25 
orientation = out
bgy1 = 0.48
bgy2 = 0.58
bgc1 = red
bgc2 = blue
<rules>
<rule>
condition = var(value) > 0
color     = blue 
fill_color = lblue
</rule>
<rule>
condition = var(value) <= 0
color     = red
fill_color = lred
</rule>
</rules>
<backgrounds>
<background>
color = vvlconf(.,bgc2)
y1    = conf(.,bgy1)
y0    = 0
</background>
<background>
color = vvlconf(.,bgc1)
y1    = 0
y0    = -conf(.,bgy1)
</background>
</backgrounds>
<axes>
<axis>
color     = grey_a1
thickness = 2
spacing   = 0.25r
</axis>
</axes>
<axis>
color     = white
thickness = 5
position  = -conf(.,bgy2),-conf(.,bgy1),conf(.,bgy1),conf(.,bgy2)
</axis>
</plot>

</plots>

## Show Gene Bodies 
<highlights>
z = 0
fill_color = grey
<highlight>
file       = highlights_forward.txt
r0         = 0.85r
r1         = 0.88r
</highlight>
<highlight>
file       = highlights_reverse.txt
r0         = 0.73r
r1         = 0.75r
</highlight>
</highlights>

<<include etc/housekeeping.conf>>
data_out_of_range* = trim
```