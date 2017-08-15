# scPipe


<img src=inst/scPipe.png height="200">

scPipe is an R package that allows barcode demultiplexing, transcript mapping and quality control of raw sequencing data generated by multiple 3 prime end sequencing protocols that include CEL-seq, MARS-seq and Drop-seq. scPipe produces a count matrix that is essential for downstream analysis along with a user-friendly HTML report that summarises data quality. These results can be used as input for downstream analyses including normalization, visualization and statistical testing.

The package is under active development, feel free to ask any question or submit a pull request.

## Installation

### from bioconductor

coming soon

### from github

```{r}
  install.packages("devtools")
  devtools::install_github("LuyiTian/scPipe")
```

## Getting started

the general workflow of scPipe:

<img src=inst/workflow.png height="800">

### data preprocessing



* `sc_trim_barcode` will reformat the read and put the cell barcode and UMI into fastq read names: `@ACGATCGA_TAGAGC#SIMULATE_SEQ::002::000::0000::0
AAGACGTCTAAGGGCGGTGTACACCCTTTTGAGCAATGATTGCACAACCTGCGATCACCTTATACAGAATTAT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA`

* after alignment. The `sc_exon_mapping` will put the cell barcode and UMI into bam file with different tags, together with gene information: `AAAGTCAA_AACTCA#SIMULATE_SEQ::007::000::0013::10        0       ERCC-00171      142     40      73M     *       0       0       GCCTCGGGAATAAGCTGACGGTGACAAGGTTTCCCCCTAATCGAGACGCTGCAATAACACAGGGGCATACAGT       AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA       HI:i:1  NH:i:1  NM:i:0  GE:Z:ERCC-00171 YC:Z:AAAGTCAA   YM:Z:AACTCA     YE:i:-364`. In this example the cell barcode is AAAGTCAA with tag `YC`, UMI is AACTCA with tag `YM` and the gene that this read mapped to is `ERCC-00171` within tag `GE`. This read located 364 bp upstream of transcription end site (TES), which is stored in `YE` tag.

* The `sc_demultiplex` will looking for the cell barcode in bam file (by default in `YC` tag). And compare it against the known cell barcode annotation file, which is a csv file consist of two columns. The first column is the cell name and second column is the cell barcode. For Drop-seq data we can run `sc_detect_bc` to find the barcode and generate the cell barcode annotation file before running `sc_demultiplex`. Here is an example barcode annotation file: `system.file("extdata", "barcode_anno.csv", package = "scPipe")`. the output of `sc_demultiplex` will be multiple csv files corresponding to each cell. The csv file will have three columns where the first column is gene id, second column is UMI sequence and third column is the relative location to TES. These files are used for `sc_gene_counting`.

for the example codes see vignettes

## TODO

* update UMI correction methods.
* use (singlecellExperiments)[https://github.com/drisso/SingleCellExperiment] package.

## Acknowledgements
This package is inspared by `scater` and `scran` packages. The idea to put cell barcode and UMI into BAM file comes from the (Drop-seq tools)[http://mccarrolllab.com/dropseq/]
