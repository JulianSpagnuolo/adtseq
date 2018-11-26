\name{adtseq-package}
\alias{adtseq-package}
\alias{adtseq}
\docType{package}
\title{
  Identification of Antibody-Derived-Tags in single cell CITE/REAP-seq data
}
\description{
  adtseq provides simple Rcpp methods to identify Antibody-Derived-Tags in single cell CITE/REAP-seq data.
  Currently the function is setup to accept the 15 nt Total-Seq antibody barcodes provided and distributed by BioLegend.
}
\details{
  This section should provide a more detailed overview of how to use the
  package, including the most important functions.
}
\author{
Julian Spagnuolo

Maintainer: Julian Spagnuolo <julian.spagnuolo@unibas.ch>
}
\references{
  
}
\keyword{ package }
\seealso{
  
}
\examples{
  \dontrun{
     adtseq(bamFileName = "ADT_processed_u.sam", bamOut = "out.sam",
     adtFasta = "total_seq_panel.fasta", max_dist = 4, sumoutput = "summary.txt")
  }
}
