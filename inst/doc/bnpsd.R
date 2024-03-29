## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# this package
library(bnpsd)
# for nice colors
library(RColorBrewer)
# for visualizing coancestry matrix with plot_popkin, and other handy functions
library(popkin)

# dimensions of data/model
# number of individuals (NOTE this is 10x less than in publication!)
n_ind <- 100
# number of intermediate subpops
k_subpops <- 10

# define population structure
# subpopulation FST vector, up to a scalar
inbr_subpops <- 1 : k_subpops
# desired bias coefficient
bias_coeff <- 0.5
# desired FST for the admixed individuals
Fst <- 0.1
# admixture proportions from 1D geography
obj <- admix_prop_1d_linear(
    n_ind = n_ind,
    k_subpops = k_subpops,
    bias_coeff = bias_coeff,
    coanc_subpops = inbr_subpops,
    fst = Fst
)
admix_proportions <- obj$admix_proportions
# rescaled inbreeding vector for intermediate subpopulations
inbr_subpops <- obj$coanc_subpops

# get pop structure parameters of the admixed individuals
coancestry <- coanc_admix(admix_proportions, inbr_subpops)

# verify that we got the desired FST!
fst_admix(admix_proportions, inbr_subpops)

# the mean of the diagonal of the coancestry matrix also equals FST
inbr <- diag(coancestry)
fst(inbr) # `fst` is a function in the popkin package

# verify that we got the desired `bias_coeff` too!
mean(coancestry) / Fst

## ---- fig.width = 4, fig.height = 2, fig.align = 'center'---------------------
# visualize the per-subpopulation inbreeding coefficients (FSTs)
# shrink default margins
par(mar = c(4, 4, 0, 0) + 0.2)
# colors for independent subpopulations
col_subpops <- brewer.pal(k_subpops, "Paired")

barplot(inbr_subpops, col = col_subpops, names.arg = 1 : k_subpops, ylim = c(0, 1),
	xlab = 'Subpopulation', ylab = 'Inbreeding coeff.')

## ---- fig.width = 6, fig.height = 1.5, fig.align = 'center'-------------------
# visualize the admixture proportions
# shrink default margins
par(mar = c(1, 4, 0, 0) + 0.2)
barplot(
    t(admix_proportions),
    col = col_subpops,
    border = NA,
    space = 0,
    ylab = 'Admixture prop.'
)
mtext('Individuals', 1)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
# Visualize the coancestry matrix using "popkin"!
plot_popkin(
    coancestry,
    # zero inner margin (plus padding) because we have no labels
    mar = 0,
    leg_title = 'Coancestry'
)

## -----------------------------------------------------------------------------
# number of loci in simulation (NOTE this is 30x less than in publication!)
m_loci <- 10000
# draw all random Allele Freqs (AFs) and genotypes
# reuse the previous inbr_subpops, admix_proportions
out <- draw_all_admix(
    admix_proportions = admix_proportions,
    inbr_subpops = inbr_subpops,
    m_loci = m_loci,
    # NOTE by default p_subpops and p_ind are not returned, but here we will ask for them
    want_p_subpops = TRUE,
    # NOTE: asking for `p_ind` increases memory usage substantially,
    # so don't ask for it unless you're sure you want it!
    want_p_ind = TRUE
)
# genotypes
X <- out$X
# ancestral AFs
p_anc <- out$p_anc
# intermediate independent subpopulation AFs
p_subpops <- out$p_subpops
# individual-specific AFs
p_ind <- out$p_ind

## ---- fig.width = 4, fig.height = 2, fig.align = 'center'---------------------
# inspect distribution of ancestral AFs (~ Uniform(0.01, 0.5))
# shrink default margins for these figures
par(mar = c(4, 4, 0, 0) + 0.2)
hist(p_anc, xlab = 'Ancestral AF', main = '', xlim = c(0, 1))

# distribution of intermediate population AFs
# (all subpopulations combined)
# (will be more dispersed toward 0 and 1 than ancestral AFs)
hist(p_subpops, xlab = 'Intermediate Subpopulation AF', main = '', xlim = c(0, 1))

# distribution of individual-specific AFs (admixed individuals)
# (admixture reduces differentiation, so these resemble ancestral AFs a bit more)
hist(p_ind, xlab = 'Individual-specific AF', main = '', xlim = c(0, 1))

# genotype distribution of admixed individuals
barplot(table(X), xlab = 'Genotypes', ylab = 'Frequency', col = 'white')

## ---- fig.width = 6, fig.height = 2.8, fig.align = 'center'-------------------
# for best estimates, group individuals into subpopulations using the geography
# this averages more individuals in estimating the minimum kinship
subpops <- ceiling( ( 1 : n_ind ) / n_ind * k_subpops )
table(subpops) # got k_subpops = 10 with 100 individuals each
# now estimate kinship using popkin
kinship_estimate <- popkin(X, subpops)
# replace diagonal with inbreeding coeffs. to match coancestry matrix
coancestry_estimate <- inbr_diag(kinship_estimate)

# Visualize the coancestry matrix using "popkin"!
plot_popkin(
    list( coancestry, coancestry_estimate ),
    titles = c('Truth', 'Estimate'),
    # second value is top margin, for panel titles
    mar = c(0, 2.5),
    leg_title = 'Coancestry'
)

## -----------------------------------------------------------------------------
coancestry_subpops_estimate <- popkin_af( p_subpops )
coancestry_ind_estimate <- popkin_af( p_ind )

## ---- fig.width = 6, fig.height = 2.8, fig.align = 'center'-------------------
plot_popkin(
    list( diag( inbr_subpops ), coancestry_subpops_estimate ),
    titles = c('Truth', 'Estimate'),
    mar = c(0, 2.5),
    leg_title = 'Coancestry'
)

## ---- fig.width = 6, fig.height = 2.8, fig.align = 'center'-------------------
plot_popkin(
    list( coancestry, coancestry_ind_estimate ),
    titles = c('Truth', 'Estimate'),
    mar = c(0, 2.5),
    leg_title = 'Coancestry'
)

## -----------------------------------------------------------------------------
# reuse the previous m_loci, inbr_subpops, admix_proportions
# draw ancestral AFs
p_anc <- draw_p_anc(m_loci)
# draw intermediate independent subpopulations AFs
p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
# draw individual-specific AFs
p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
# draw genotypes
X <- draw_genotypes_admix(p_ind)

## ---- fig.width = 4, fig.height = 2, fig.align = 'center'---------------------
# this increases the proportion of rare alleles
p_anc_beta <- draw_p_anc(m_loci, beta = 0.1)
# shrink default margins for these figures
par(mar = c(4, 4, 0, 0) + 0.2)
hist(p_anc_beta, xlab = 'Ancestral AF', main = '', xlim = c(0, 1))

## ---- fig.width = 4, fig.height = 2, fig.align = 'center'---------------------
# this increases the proportion of rare alleles
p_anc_beta <- rbeta(m_loci, shape1 = 0.1, shape2 = 1)
# shrink default margins for these figures
par(mar = c(4, 4, 0, 0) + 0.2)
hist(p_anc_beta, xlab = 'Ancestral AF', main = '', xlim = c(0, 1))

## -----------------------------------------------------------------------------
# draw genotypes for one individual from the ancestral population
# use "cbind" to turn the vector p_anc into a column matrix
# ("draw_genotypes_admix" expects a matrix)
X_anc <- draw_genotypes_admix( cbind( p_anc ) )
# returns a column matrix:
dim(X_anc)

# draw genotypes from intermediate populations
# draws one individual per intermediate population
X_subpops <- draw_genotypes_admix(p_subpops)

## -----------------------------------------------------------------------------
# data dimensions
# number of individuals
n_ind <- 100
# number of intermediate subpops
k_subpops <- 10

# define population structure
# subpopulation FST vector, up to a scalar
inbr_subpops <- 1 : k_subpops
# desired bias coefficient
bias_coeff <- 0.5
# desired FST for the admixed individuals
Fst <- 0.1

# admixture proportions from *circular* 1D geography
obj <- admix_prop_1d_circular(
    n_ind = n_ind,
    k_subpops = k_subpops,
    bias_coeff = bias_coeff,
    coanc_subpops = inbr_subpops,
    fst = Fst
)
admix_proportions <- obj$admix_proportions
inbr_subpops <- obj$coanc_subpops

# get pop structure parameters of the admixed individuals
coancestry <- coanc_admix(admix_proportions, inbr_subpops)

# verify that we got the desired FST!
fst_admix(admix_proportions, inbr_subpops)

# verify that we got the desired bias_coeff too!
mean(coancestry) / Fst

## ---- fig.width = 4, fig.height = 1.2, fig.align = 'center'-------------------
# visualize the per-subpopulation inbreeding coefficients (FSTs)
# tweak margins/etc
par(mar = c(2.5, 2.5, 0.3, 0) + 0.2, lab = c(2, 1, 7), mgp = c(1.5, 0.5, 0))
# colors for independent subpopulations
col_subpops <- brewer.pal(k_subpops, "Paired")
barplot(inbr_subpops, col = col_subpops, names.arg = colnames(admix_proportions), xlab = 'Subpopulation', ylab = 'Inbr')

## ---- fig.width = 4, fig.height = 1, fig.align = 'center'---------------------
# visualize the admixture proportions
# tweak margins/etc
par(mar = c(1, 4, 0.4, 0) + 0.2, lab = c(2, 2, 7))
barplot(
    t(admix_proportions),
    col = col_subpops,
    border = NA,
    space = 0,
    ylab = 'Admix prop'
)
mtext('Individuals', 1)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
# Visualize the coancestry matrix using "popkin"!
plot_popkin(
    coancestry,
    leg_n = 3,
    mar = c(0, 0.4),
    leg_title = 'Coancestry'
)

## -----------------------------------------------------------------------------
# define population structure
# we'll have k_subpops = 3, each subpopulation with these sizes:
n1 <- 100
n2 <- 50
n3 <- 20

# here's the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
labs <- c(
    rep.int('S1', n1),
    rep.int('S2', n2),
    rep.int('S3', n3)
)
# data dimensions infered from labs:

# number of individuals "n_ind"
length(labs)

# number of subpopulations "k_subpops"
k_subpops <- length(unique(labs))
k_subpops

# desired admixture matrix
admix_proportions <- admix_prop_indep_subpops(labs)

# got a numeric matrix with a single 1 value per row
# (denoting the sole subpopulation from which each individual draws its ancestry)
head(admix_proportions, 2)

# construct the intermediate subpopulation FST vector
# the desired final FST
Fst <- 0.2
# subpopulation FST vector, unnormalized so far
inbr_subpops <- 1 : k_subpops
# normalized to have the desired FST
# NOTE fst is a function in the `popkin` package
inbr_subpops <- inbr_subpops / fst(inbr_subpops) * Fst
# verify FST for the intermediate subpopulations
fst(inbr_subpops)

# get coancestry of the admixed individuals
coancestry <- coanc_admix(admix_proportions, inbr_subpops)
# before getting FST for individuals, weigh then inversely proportional to subpop sizes
weights <- weights_subpops(labs) # function from `popkin` package
# verify FST for individuals (same as for intermediate subpops for this pop structure)
fst_admix(admix_proportions, inbr_subpops, weights)

## ---- fig.width = 4, fig.height = 1.2, fig.align = 'center'-------------------
# visualize the per-subpopulation inbreeding coefficients (FSTs)
# tweak margins/etc
par(mar = c(2.5, 2.5, 0, 0) + 0.2, lab = c(2, 1, 7), mgp = c(1.5, 0.5, 0))
# colors for independent subpopulations
col_subpops <- brewer.pal(k_subpops, "Paired")
barplot(inbr_subpops, col = col_subpops, names.arg = colnames(admix_proportions), xlab = 'Subpopulation', ylab = 'Inbr')

## ---- fig.width = 4, fig.height = 1, fig.align = 'center'---------------------
# visualize the admixture proportions
# tweak margins/etc
par(mar = c(1, 4, 0.4, 0) + 0.2, lab = c(2, 2, 7))
barplot(
    t(admix_proportions),
    col = col_subpops,
    border = NA,
    space = 0,
    ylab = 'Admix prop'
)
mtext('Individuals', 1)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
# Visualize the coancestry matrix using "popkin"!
plot_popkin(
    coancestry,
    leg_n = 3,
    mar = c(0, 0.4),
    leg_title = 'Coancestry'
)

## -----------------------------------------------------------------------------
tree_str <- '(S1:0.1,(S2:0.1,(S3:0.1,(S4:0.1,S5:0.1)N3:0.1)N2:0.1)N1:0.1)T;'

## ---- fig.width = 4, fig.height = 3.5, fig.align = 'center'-------------------
# base package for phylogenetics
library(ape)
# parses tree, creates object of class "phylo"
tree_subpops <- read.tree( text = tree_str )
# plot tree
par( mar = c(3,0,0,0) + 0.2 )
plot( tree_subpops, show.node.label = TRUE )
# add axis and label
axisPhylo( backward = FALSE )
mtext( 'Coancestry', side = 1, line = 2 )

## -----------------------------------------------------------------------------
p_subpops_tree <- draw_p_subpops_tree(
    p_anc = p_anc,
    tree_subpops = tree_subpops
)

## ---- fig.width = 6, fig.height = 2.8, fig.align = 'center'-------------------
# estimate coancestry
coancestry_subpops_tree_estimate <- popkin_af( p_subpops_tree )
# expected coancestry according to model
coancestry_subpops_tree <- coanc_tree( tree_subpops )
# plot
plot_popkin(
    list( coancestry_subpops_tree, coancestry_subpops_tree_estimate ),
    titles = c('Truth', 'Estimate'),
    mar = c(2, 2.5),
    leg_title = 'Coancestry',
    names = TRUE
)

## -----------------------------------------------------------------------------
1 - (1 - 0.1)^4

## -----------------------------------------------------------------------------
# number of admixed individuals
n_ind <- 100
# number of subpopulations, extracted from true coancestry
k_subpops <- nrow( coancestry_subpops_tree )

# define population structure
# for simplicity here we set spread parameter `sigma` directly
# admixture proportions from 1D geography
admix_proportions <- admix_prop_1d_linear(
    n_ind = n_ind,
    k_subpops = k_subpops,
    sigma = 0.5
)

# get pop structure parameters of the admixed individuals
coancestry <- coanc_admix(admix_proportions, coancestry_subpops_tree)

# calculate true FST of simulation
fst_admix(admix_proportions, coancestry_subpops_tree)

## ---- fig.width = 6, fig.height = 1.5, fig.align = 'center'-------------------
# visualize the admixture proportions
# colors for these independent subpopulations
col_subpops <- brewer.pal(k_subpops, "Paired")
# shrink default margins
par(mar = c(1, 4, 0, 0) + 0.2)
barplot(
    t(admix_proportions),
    col = col_subpops,
    border = NA,
    space = 0,
    ylab = 'Admixture prop.'
)
mtext('Individuals', 1)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
# Visualize the coancestry matrix using "popkin"!
plot_popkin(
    coancestry,
    # zero inner margin (plus padding) because we have no labels
    mar = 0,
    leg_title = 'Coancestry'
)

## -----------------------------------------------------------------------------
# number of loci in simulation (NOTE this is 30x less than in publication!)
m_loci <- 10000
# draw all random Allele Freqs (AFs) and genotypes
# reuse the previous inbr_subpops, admix_proportions
out <- draw_all_admix(
    admix_proportions = admix_proportions,
    tree_subpops = tree_subpops,
    m_loci = m_loci
)
# genotypes
X <- out$X

## ---- fig.width = 6, fig.height = 2.8, fig.align = 'center'-------------------
# estimate kinship using popkin
# here we didn't use subpopulation labels for simplicity
kinship_estimate <- popkin(X)
# replace diagonal with inbreeding coeffs. to match coancestry matrix
coancestry_estimate <- inbr_diag(kinship_estimate)

# Visualize the coancestry matrix using "popkin"!
plot_popkin(
    list( coancestry, coancestry_estimate ),
    titles = c('Truth', 'Estimate'),
    # second value is top margin, for panel titles
    mar = c(0, 2.5),
    leg_title = 'Coancestry'
)

