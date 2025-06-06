# Dalia Bornstein - Part II Project

This is a Github supporting Dalia Bornstein's Part II Project which contains computational methods and code as described in my part II report.

---

## Illumina 16S Amplicon Analysis
This section describes the analysis framework (mostly in R/Bioconductor) for processing and analysing
the 16S amplicon sequencing data. The raw data is from Illumina MiSeq and the runs were Paired-end 250nt.

* NEXTFLEX® 16S V1 - V3 rRNA Amplicon-Seq Library Prep Kit
This [kit](https://www.revvity.com/gb-en/product/nextflex-16s-v1-v3-12-24-rxn-nova-4202-02s#product-overview) is designed for the preparation of multiplexed amplicon libraries that span the hypervariable domains one through three (V1-V3) of microbial 16S ribosomal RNA (rRNA) genes. These libraries are compatible with paired-end sequencing on the Illumina® sequencing platforms.
* Sequencer Used Illumina MiSeq at Cambridge Genomic Services (Pathology).

### DADA2 Filtering Options 

*DaDa2*[^1] was used to process the raw illumina sequencing data in order to detect amplicon sequence variants (ASVs) and to compare them to a reference databases.
The key options used are below but also available in the Rmarkdown (below). Additonal analysis of ASVs across taxa are performed using *phyloseq*[^3]. Again the Rmarkdown code for this is below.

* `truncLen=c(221,141)` - Filter forward reads by removing 29nt and 110nt from the reverse
* `maxN=0` - Maximum number of 'N' bases before filtering
* `truncQ=2` - Truncate reads with Quality threshold 2
* `rm.phix=TRUE` - Remove reads which are PhiX contaminants
* `compress=TRUE` - Compress the output (gz)
* `multithread=TRUE` - Use maximum available CPUs

### 16S Database Used

We used the [SILVA database](https://benjjneb.github.io/dada2/training.html)[^2] as a reference for *DaDa2*.

* Training Database - SILVA non-redundant version v132 `silva_nr_v132_train_set.fa`
* Species Database - SILVA species assignment v132 `silva_species_assignment_v132.fa`

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
|Bream   | 8.0.12 | 
|Configuration | 6.0.15 |
|Dorado  | 7.4.13 | 
|MinKNOW | Core 6.0.11 |

### Metagenomic Mapping Tools Used

*Kraken2* was used to employ fast *k-mer* indexing to map metagenome nanopore reads rapidly against the *core_nt* database. It produces *k-mer* counts on taxonomic nodes.
The *Bracken* tool aims to estimate actual abundance from the raw *Kraken2* counts.

-   Kraken version 2.1.3[^4].
-   Bracken v2.9[^5].

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

### Bracken Abundance Estimation
```
est_abundance.py -i barcodeX.kraken2.txt -k kraken_db/database75mers.kmer_distrib -l S -t 10 -o barcodeX.bracken.txt
```

* `-i barcodeX.kraken2.txt` (input file from kraken)
* `-k kraken_db/database75mers.kmer_distrib` (k-mers file from database)
* `-l S` (species level calls)
* `-t 10` (use 10 cpus)
* `-o barcodeX.bracken.txt` (output abundance information in new report file)

### Metagenomic Tables

The metagenomic tables from kraken and bracken are in the [tables folder](tables/kraken_bracken).
Some example output is shown below:

```
 49.11	3272876	3272876	U	0	unclassified
 50.89	3391970	1039443	R	1	root
 34.27	2283844	30236	R1	131567	  cellular organisms
 33.81	2253183	292082	D	2	    Bacteria
 19.27	1284120	24498	D1	1783272	      Terrabacteria group
 14.38	958583	46700	P	1239	        Bacillota
 13.42	894704	80116	C	186801	          Clostridia
  6.62	441032	36705	O	186802	            Eubacteriales
  5.76	383834	116272	F	216572	              Oscillospiraceae
  2.64	176082	6	F1	552397	                Oscillospiraceae incertae sedis
  2.64	176076	173949	S	39492	                  [Eubacterium] siraeum
  0.02	1471	1471	S1	717961	                    [Eubacterium] siraeum V10Sc8a
  0.01	656	656	S1	657319	                    [Eubacterium] siraeum 70/3
  0.85	56760	41265	G	216851	                Faecalibacterium
  0.18	12135	11856	S	853	                  Faecalibacterium prausnitzii
  0.00	264	264	S1	718252	                    Faecalibacterium prausnitzii L2-6
  0.00	15	15	S1	657322	                    Faecalibacterium prausnitzii SL3/3
  0.05	3267	1885	G1	2646395	                  unclassified Faecalibacterium
  0.01	630	630	S	2929493	                    Faecalibacterium sp. I3-3-89
  0.00	330	330	S	2929489	                    Faecalibacterium sp. IP-3-29
  0.00	144	144	S	2929495	                    Faecalibacterium sp. I4-3-84
  0.00	140	140	S	2929490	                    Faecalibacterium sp. I2-3-92
  0.00	68	68	S	3141184	                    Faecalibacterium sp. i21-0019-B1
  0.00	49	49	S	3141185	                    Faecalibacterium sp. i25-0019-C1
  0.00	12	12	S	2929494	                    Faecalibacterium sp. I4-1-79
  0.00	7	7	S	2929491	                    Faecalibacterium sp. HTF-F
  0.00	1	1	S	2929488	                    Faecalibacterium sp. IP-1-18
  0.00	1	1	S	2929492	                    Faecalibacterium sp. I3-3-33
  ...
  ```

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

### R Analysis of metagenomic Data
Some R scripts are used to process and analyse the nanopore metagenomic data from bracken and kraken.

You can see these scripts in the [nanopore_count_levels](nanopore_count_levels.md) files.

An example plot from R showing the number of mapped, unmapped and the different taxonomic levels reached is shown below:
![nanopore_counts](nanopore_count_levels_files/figure-gfm/unnamed-chunk-2-1.png)

---

## Metagenome Visualisation

### Krona Plots
Individual Bracken report files (or Kraken report files) can be turned into a *Krona*[^6] HTML visualisation.
For this we use the following command from `KronaTools v2.8.1`:

```
ktImportTaxonomy -t 5 -m 3 -o barcodeX.html barcodeX.bracken.report.txt
```

This uses 5 threads and generates a html report (barcodeX.html) from the input report.

The multi-sample plot of *bracken* abundances from *Krona* is here [dalia-krona.html](https://enrightlab.github.io/gutmetagenome/krona/dalia-krona.html)
An example image from *sample 1* is below:
![krona](krona/krona.png)


### Circos Plots

* Circos Version Used: `circos 0.69-9`[^7].

For genomes of interest we download a reference in *Genbank* format, e.g. `Akkermansia.gb` is the Genbank entry (CP001071.1) for *Akkermansia sp*.
We then use a perl script [find_species_hits.pl](scripts/find_species_hits.pl) to interrogate all kraken mapping files to identify reads hitting the species of interest.
These reads are extracted into a new *FASTA* file. We also extract a genome *FASTA* file from the same genbank file.

To identify per genome mappings we use *blastn*[^8] to map the reads against the reference genome *FASTA*.

* `blast 2.16.0, build Jun 25 2024 12:36:54`

### Genbank References Used
The genbank reference data use is in the folder provided [here](circos/genbank).

|Species (and strain)|Genbank Accession|Size|Type|Origin|Date|
|--------------------|-----------------|----|----|----|------|
|Akkermansia muciniphila ATCC BAA-835|[CP001071](circos/genbank/akkermansia.gb)|2664102 bp|DNA circular|BCT|16-AUG-2022|
|Candida parapsilosis strain CDC317 Chr 1|[NW_023503284](circos/genbank/candida.gb)|957321 bp|DNA linear|CON|28-OCT-2020|
|Candida parapsilosis strain CDC317 Chr 2|[NW_023503283](circos/genbank/candida.gb)|1789679 bp|DNA linear|CON|28-OCT-2020|
|Candida parapsilosis strain CDC317 Chr 3|[NW_023503282](circos/genbank/candida.gb)|962442 bp|DNA linear|CON|28-OCT-2020|
|Candida parapsilosis strain CDC317 Chr 4|[NW_023503281](circos/genbank/candida.gb)|3023470 bp|DNA linear|CON|28-OCT-2020|
|Candida parapsilosis strain CDC317 Chr 5|[NW_023503280](circos/genbank/candida.gb)|2091826 bp|DNA linear|CON|28-OCT-2020|
|Candida parapsilosis strain CDC317 Chr 6|[NW_023503279](circos/genbank/candida.gb)|1039767 bp|DNA linear|CON|28-OCT-2020|
|Candida parapsilosis strain CDC317 Chr 7|[NW_023503278](circos/genbank/candida.gb)|2235583 bp|DNA linear|CON|28-OCT-2020|
|Candida parapsilosis strain CDC317 Chr 8|[NW_023503277](circos/genbank/candida.gb)|898305 bp|DNA linear|CON|28-OCT-2020|
|Candida parapsilosis strain CDC317 Chr MT|[NC_005253](circos/genbank/candida.gb)|32745 bp|DNA linear|PLN|03-APR-2023|
|Herelleviridae sp. isolate ctX5e1|[BK045609](circos/genbank/herelleviridae.gb)|29362 bp|DNA linear|ENV|24-JUN-2021|
|Methanobrevibacter smithii ATCC 35061|[CP000678](circos/genbank/methano_smithii.gb)|1853160 bp|DNA circular|BCT|31-JAN-2014|
|Candidatus Methanomassiliicoccus intestinalis isolate 138|[LVVS01000017](circos/genbank/methanomassiliicoccus.gb)|461488 bp|DNA linear|ENV|12-JUL-2019|
|Methanosphaera stadtmanae DSM 3091|[CP000102](circos/genbank/methanosphaera.gb)|1767403 bp|DNA circular|BCT|31-JAN-2014|
|Sutterella wadsworthensis strain FDAARGOS_1159|[CP068055](circos/genbank/suterella.gb)|3026517 bp|DNA circular|BCT|20-JAN-2021|


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
...
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

#### Example labels_forward track
```
chr1    1       234     mobA
chr1    608     2142    rRNA 16S
chr1    5589    5703    rRNA
chr1    7442    9463    uvrB
chr1    12165   12572   iscU
chr1    12578   12901   iscA
chr1    12915   13409   hscB
chr1    13500   15386   hscA
chr1    15428   15760   fdx
...
```
#### Example data_tile track
```
chr1    0       999     0.275594405594406
chr1    1000    1999    0.41734965034965
chr1    2000    2999    0.646265734265734
chr1    3000    3999    0.727986013986014
chr1    4000    4999    0.677216783216783
chr1    5000    5999    0.655524475524475
chr1    6000    6999    0.655258741258741
...
```

#### Example highlight_forward track
```
chr1    608     2142    fill_color=30,30,30
chr1    2228    2304    fill_color=30,30,30
chr1    2597    5477    fill_color=30,30,30
chr1    5589    5703    fill_color=30,30,30
chr1    7442    9463    fill_color=94,252,130
chr1    9585    10124   fill_color=30,30,30
chr1    10343   10828   fill_color=30,30,30
chr1    10865   12076   fill_color=30,30,30
chr1    12165   12572   fill_color=142,173,119
chr1    12578   12901   fill_color=102,122,106
chr1    12915   13409   fill_color=189,151,146
...
```

#### Command Line
This perl script [make_karyotype_gb.pl](scripts/make_karyotype_gb.pl) computes the key track files:

```
./make_karyotype_gb.pl akkermansia.gb akkermansia.hits akkermansia.skew --names; circos -conf metagenome.conf
```
This script makes use of [genbank_gtf.pl](scripts/genbank_gtf.pl) by Jiang Li (Vanderbilt Center for Quantitative Sciences).


This creates a *png* and an *svg* vector graphic for each species of interest. The *Circos* plots are available [circos_plots](./circos).
![akkermansia](circos/akkermansia.png)
![candida](circos/candida.png)
![herelleviridae](circos/herelleviridae.png)
![methano_smithii](circos/methano_smithii.png)
![methanomassiliicoccus](circos/methanomassiliicoccus.png)
![methanosphaera](circos/methanosphaera.png)
![suterella](circos/suterella.png)


All species were run with the following simple shell script:
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

---

## Citations

[^1]:Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP; DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods. (2016) Jul;13(7):581-3. doi: https://10.1038/nmeth.3869.
[^2]:Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO; The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. (2013) Nucleic Acids Res.
[^3]:McMurdie and Holmes; phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE. (2013) 8(4):e61217.
[^4]:Lu J, Rincon N, Wood DE, Breitwieser FP, Pockrandt C, Langmead B, Salzberg SL, Steinegger M; Metagenome analysis using the Kraken software suite. Nat Protoc. (2022) Dec;17(12):2815-2839. doi: https://10.1038/s41596-022-00738-y. Epub 2022 Sep 28. Erratum in: Nat Protoc. 2024 Aug 29. doi: 10.1038/s41596-024-01064-1. PMID: 36171387; PMCID: PMC9725748.
[^5]:Jennifer Lu, Florian P Breitwieser, Peter Thielen, Steven L Salzberg; Bracken: Estimating species abundance in metagenomics data. bioRxiv (2016) 051813; doi: https://doi.org/10.1101/051813.
[^6]:Ondov, B.D., Bergman, N.H. & Phillippy, A.M; Interactive metagenomic visualization in a Web browser. BMC Bioinformatics 12, 385 (2011). https://doi.org/10.1186/1471-2105-12-385.
[^7]:Krzywinski, M. et al; Circos: an Information Aesthetic for Comparative Genomics. Genome Res (2009) 19:1639-1645.
[^8]:Altschul, S.F., Gish, W., Miller, W., Myers, E.W. and Lipman, D.J.;  Basic local alignment search tool. Journal of Molecular Biology, (1990) 215(3), pp.403-410.
