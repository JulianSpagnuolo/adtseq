
<!-- README.md is generated from README.Rmd. Please edit that file -->

# adtseq

adtseq provides simple Rcpp methods to identify Antibody-Derived-Tags in
single cell CITE/REAP-seq data. Currently the function is setup to
accept the 15 nt Total-Seq antibody barcodes provided and distributed by
BioLegend.

## Installation

Currently adtseq requires the Seqan headers provided by the RSeqan
package, the pacakge currently links to the Boost headers provided by
the BH package, thus it must also be installed - these are primarily
used in some utility diagnostic functions (in the dev branch). Please
refer to each of these packages for their installation requirements.

NOTE: the primary function will be able to accept BAM files only if the
Seqan headers are linked to zlib (should be present on most \*nix type
systems).

``` r
install.packages(c("Rcpp", "BH"))
BiocManager::install("RSeqan")
install_github("JulianSpagnuolo/adtseq")
```

### Citation

It is requested that any work published using this package or scripts
contained within cite adtseq

    #> 
    #> To cite package 'adtseq' in publications use:
    #> 
    #>   Julian Spagnuolo (2018). adtseq: Identification of
    #>   Antibody-Derived-Tags in single cell CITE/REAP-seq data. R
    #>   package version 0.99.3.
    #> 
    #> A BibTeX entry for LaTeX users is
    #> 
    #>   @Manual{,
    #>     title = {adtseq: Identification of Antibody-Derived-Tags in single cell
    #> CITE/REAP-seq data},
    #>     author = {Julian Spagnuolo},
    #>     year = {2018},
    #>     note = {R package version 0.99.3},
    #>   }
    #> 
    #> ATTENTION: This citation information has been auto-generated from
    #> the package DESCRIPTION file and may need manual editing, see
    #> 'help("citation")'.

## Example

The primary function `adtseq` will take two files as input; A) a fasta
file containing the antibody derived barcodes to identify in the
sequencing data, B) a bam or sam file that has been processed through
the dropseq-tools pipeline to tag read 2 reads with the cell barcodes
and UMIs, remove the SMART-PCR adapter sequences and (optionally) 3’
polyA stretches. It is intended that this function be used at command
line as it is currently extremely verbose and will output cell barcodes,
UMI, antibody-ID, hamming distance and map score directly to std output.

The function will tag each read with the antibody barcode name (form the
fasta file) using a combination of hamming distance (less than
max\_dist) and best Phred quality weighted mapping score calculated as
the sum of scores for each base, subtracting the score at mismatches or
adding at matches. This is performed similarly to the method used by
Bowtie2 to determine alignment scores.

Score = MN + ((MX - MN) \* (Phred/40))

At present Bowtie2’s defaults for MN (minimum score) and MX (maximum
score), 2 and 6, respectively, are used.

The function will place the name of the identified antibody in a new tag
“XA”, the mapping score under “AS”, and the hamming distance under “XN”.
Additionally, the function will output a second tab-delimited text file
containing the Cell Barcodes, UMI and identified antibody along with the
hamming distance and mapping scores for each bam/sam record.

The resulting BAM/SAM file can be directly used in latter stages of the
dropseq-tools pipeline which attempt to correct for synthesis errors in
the cell barcodes and UMIs. However, should users wish, they may
retrieve count data directly from the text file and omit the latter
steps (not
recommended).

``` r
adtseq(bamFileName = "ADT_processed_u.sam", bamOut = "out.sam", adtFasta = "total_seq_panel.fasta", max_dist = 4, sumoutput = "summary.txt", adt_panel = "A", bc_length = 15)
```

## NOTE:

*Experimental* support for TotalSeq panels B and C has been implemented
but not tested. In addition, support for custom barcodes has also been
implemented, this is controlled via the “adt\_panel” parameter. The
bc\_length parameter must be set for custom panels.

Default is to use the TotalSeq Panel “A” regex rules and a barcode
length (bc\_length) of 15.

``` r
# TotalSeq Panel A
adtseq(bamFileName = "ADT_processed_u.sam", bamOut = "out.sam", adtFasta = "total_seq_panel.fasta", max_dist = 4, sumoutput = "summary.txt", adt_panel = "A", bc_length = 15)
# TotalSeq Panel B
adtseq(bamFileName = "ADT_processed_u.sam", bamOut = "out.sam", adtFasta = "total_seq_panelB.fasta", max_dist = 4, sumoutput = "summary.txt", adt_panel = "B", bc_length = 15)
# TotalSeq Panel C
adtseq(bamFileName = "ADT_processed_u.sam", bamOut = "out.sam", adtFasta = "total_seq_panelC.fasta", max_dist = 4, sumoutput = "summary.txt", adt_panel = "C", bc_length = 15)
# Custom panel
adtseq(bamFileName = "ADT_processed_u.sam", bamOut = "out.sam", adtFasta = "custom_panel.fasta", max_dist = 4, sumoutput = "summary.txt", adt_panel = "custom", bc_length = 8)
```
