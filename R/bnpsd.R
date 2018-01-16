#' A package for modeling and simulating an admixed population
#'
#' The underlying model is called the BN-PSD admixture model, which combines the Balding-Nichols (BN) allele frequency model for the intermediate subpopulations with the Pritchard-Stephens-Donnelly (PSD) model of individual-specific admixture proportions.
#' The BN-PSD model enables the simulation of complex population structures, ideal for illustrating challenges in kinship coefficient and \eqn{F_{ST}}{FST} estimation.
#' Note that simulated loci are drawn independently (in linkage equilibrium).
#'
#' @examples
#' # dimensions of data/model
#' m <- 10 # number of loci
#' n <- 5 # number of individuals
#' k <- 2 # number of intermediate subpops
#'
#' # define population structure
#' F <- c(0.1, 0.3) # FST values for k=2 subpopulations
#' sigma <- 1 # dispersion parameter of intermediate subpops
#' Q <- q1d(n, k, sigma) # admixture proportions from 1D geography
#'
#' # get pop structure parameters of the admixed individuals
#' Theta <- coanc(Q,F) # the coancestry matrix
#' Fst <- fst(Q,F) # Fst
#'
#' # draw all random allele freqs and genotypes
#' out <- rbnpsd(Q, F, m)
#' X <- out$X # genotypes
#' P <- out$P # IAFs (individual-specific AFs)
#' B <- out$B # intermediate AFs
#' pAnc <- out$Pa # ancestral AFs
#'
#' # OR... draw each vector or matrix separately
#' # provided for additional flexibility
#' pAnc <- rpanc(m) # "anc"estral AFs
#' B <- rpint(pAnc, F) # "int"ermediate AFs
#' P <- rpiaf(B, Q) # "IAF"s (individual-specific AFs)
#' X <- rgeno(P) # "geno"types
#' 
#' @docType package
#' @name bnpsd
"_PACKAGE"