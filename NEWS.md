# adtseq 0.99.3

* Added support for TotalSeq Panel B, C and custom antibody barcodes (EXPERIMENTAL)

# adtseq 0.99.2

* Altered Makevar to enable ZLIB for RSeqan
* Added documentation for the functions.
* Added several functions (bcExtract, icgrEncode, complexityEntropy)

# adtseq 0.99.1

* Added Hamming Distance Function hdist.
* Added dynamic unload utility function.
* Updated pkg Rd file.

# adtseq 0.99.0

* Added a `NEWS.md` file to track changes to the package.
* Added a `README.Rmd` file to create and modify the github ReadMe providing some description of the pkg.
* Altered the Cell BC and UMI tag extraction methods to catch instances of XC/XM tags not being found.
* Rolled back version number to 0.99.0 as per Bioconductor recommendations for pre-release packages that are feature complete