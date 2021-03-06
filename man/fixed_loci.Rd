% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixed_loci.R
\name{fixed_loci}
\alias{fixed_loci}
\title{Identify fixed loci}
\usage{
fixed_loci(X)
}
\arguments{
\item{X}{The \code{m}-by-\code{n} genotype matrix}
}
\value{
A length-\code{m} boolean vector where the \code{i} element is \code{TRUE} if locus \code{i} is fixed or completely missing, \code{FALSE} otherwise.
}
\description{
A locus is "fixed" if the non-missing sub-vector contains all 0's or all 2's (the locus is completely homozygous for one allele or completely homozygous for the other allele).
This function tests each locus, returning a vector that is \code{TRUE} for each fixed locus, \code{FALSE} otherwise.
Loci with only missing elements (\code{NA}) are treated as fixed.
Below \code{m} is the number of loci, and \code{n} is the number of individuals.
}
\examples{
# here's a toy genotype matrix
X <- matrix(
       data = c(
              2, 2, NA,  # fixed locus (with one missing element)
              0, NA, 0,  # another fixed locus, for opposite allele
              1, 1, 1,   # NOT fixed (heterozygotes are not considered fixed)
              0, 1, 2,   # a completely variable locus
              NA, NA, NA # completely missing locus (will be treated as fixed)
             ),
       ncol = 3, byrow = TRUE)

# test that we get the desired values
stopifnot(
  fixed_loci(X) == c(TRUE, TRUE, FALSE, FALSE, TRUE)
)

}
