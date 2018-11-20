complexity.entropy <- function(object) {
  
  #' @importFrom Biostrings trinucleotideFrequency
  #' @importFrom BiocGenerics width
  scores = vector(mode="numeric", length=length(object))
  tfq = trinucleotideFrequency(object) # matrix with nreads rows and 64 cols
  rls = width(object)
  fac1 = tfq / (rls - 2)
  fac2 = ifelse(fac1 == 0, 0, log(fac1, base=ifelse(rls < 66, rls - 2, 64)))
  scores = (-100) * rowSums(fac1 * fac2)
  return(scores)
}