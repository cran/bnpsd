% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coanc_admix.R
\name{coanc_admix}
\alias{coanc_admix}
\title{Construct the coancestry matrix of an admixture model}
\usage{
coanc_admix(admix_proportions, coanc_subpops)
}
\arguments{
\item{admix_proportions}{The \code{n}-by-\code{k} admixture proportion matrix}

\item{coanc_subpops}{The intermediate subpopulation coancestry, given either as a \code{k}-by-\code{k} matrix (for the complete admixture model), or the length-\code{k} vector of intermediate subpopulation FST values (for the BN-PSD model; implies zero coancestry between subpopulations), or a scalar FST value shared by all intermediate subpopulations (also implies zero coancestry between subpopulations).}
}
\value{
The \code{n}-by-\code{n} coancestry matrix.
}
\description{
The \code{n}-by-\code{n} coancestry matrix \code{Theta} of admixed individuals is determined by the \code{n}-by-\code{k} admixture proportion matrix \code{Q} and the \code{k}-by-\code{k} intermediate subpopulation coancestry matrix \code{Psi}, given by \code{Theta = Q \%*\% Psi \%*\% t(Q)}
In the more restricted BN-PSD model, \code{Psi} is a diagonal matrix (with FST values for the intermediate subpopulations along the diagonal, zero values off-diagonal).
}
\examples{
# a trivial case: unadmixed individuals from independent subpopulations
# number of individuals and subpops
n_ind <- 5
# unadmixed individuals
admix_proportions <- diag(rep.int(1, n_ind))
# equal Fst for all subpops
coanc_subpops <- 0.2
# diagonal coancestry matryx
coancestry <- coanc_admix(admix_proportions, coanc_subpops)

# a more complicated admixture model
# number of individuals
n_ind <- 5
# number of intermediate subpops
k_subpops <- 2
# non-trivial admixture proportions
admix_proportions <- admix_prop_1d_linear(n_ind, k_subpops, sigma = 1)
# different Fst for each of the k_subpops
coanc_subpops <- c(0.1, 0.3)
# non-trivial coancestry matrix
coancestry <- coanc_admix(admix_proportions, coanc_subpops)

}
