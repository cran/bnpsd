# Find the sigma that gives the desired bias coefficient
#
# This function uses a numerical solver to find the value of \code{sigma} in \code{func} that give the desired bias coefficient \eqn{s} given the other PSD parameters.
# The admixture scenarios guaranteed to work are the 1D linear and 1D circular geographies as provided in \code{func}.
# Therefore the structure of the admixture coefficients is not explicitly stated in the inputs.
# Assumes uniform weights in calculating the bias coefficient.
#
# @param bias_coeff The desired bias coefficient.
# @param coanc_subpops The coancestry matrix of the intermediate subpopulations (or equivalent vector or scalar), up to a scaling factor (which is irrelevant as it cancels out in the bias coefficient).
# @param n_ind The number of individuals.
# @param k_subpops The number of subpopulations.
# @param func A function that accepts \code{(n_ind, k_subpops, sigma, coord_ind_first = coord_ind_first, coord_ind_last = coord_ind_last)} as inputs and returns the admixture proportions matrix.
# @param coord_ind_first Location of first individual (to pass to func).
# @param coord_ind_last Location of last individual (to pass to func).
#
# @return The desired value of \code{sigma}
#
# @examples
# # number of intermediate subpops
# k <- 10
# # set bias coeff. of 1/2, coanc_subpops and n, implicitly assumes our 1D scenario
# sigma <- bias_coeff_admix_fit(
#     bias_coeff = 0.5,
#     coanc_subpops = (1:k)/k,
#     n_ind = 1000,
#     k_subpops = k,
#     func = admix_prop_1d_linear,
#     coord_ind_first = 0.5,
#     coord_ind_last = k + 0.5
# ) 
bias_coeff_admix_fit <- function(
                                 bias_coeff,
                                 coanc_subpops,
                                 n_ind,
                                 k_subpops,
                                 func,
                                 coord_ind_first,
                                 coord_ind_last
                                 ) {
    # finds the value of param that give a desired bias coefficient "s"
    # we use a very generic numeric method, since the function is complicated at best
    if (missing(bias_coeff))
        stop('Desired `bias_coeff` is required!')
    if (missing(coanc_subpops))
        stop('`coanc_subpops` is required!')
    if (missing(n_ind))
        stop('Number of individuals `n_ind` is required!')
    if (missing(k_subpops))
        stop('Number of subpopulations `k_subpops` is required!')
    if (missing(func))
        stop('Admixture proportion function `func` is required!')
    if (missing(coord_ind_first))
        stop('`coord_ind_first` is required!')
    if (missing(coord_ind_last))
        stop('`coord_ind_last` is required!')
    
    # validate range
    if (bias_coeff > 1)
        stop('Desired `bias_coeff` must be <= 1, passed ', bias_coeff)
    if (bias_coeff < 0)
        stop('Desired `bias_coeff` must be >= 0, passed ', bias_coeff)
    
    # actual minimum is greater than zero, but it seems that message should be different for non-negative values out of this range
    # set sigma = 0 here, in both cases that we have (admix_prop_1d_linear, admix_prop_1d_circular) it's the minimum bias (independent subpopulations)!
    admix_prop_bias_coeff_min <- func(
        n_ind,
        k_subpops,
        sigma = 0,
        coord_ind_first = coord_ind_first,
        coord_ind_last = coord_ind_last
    )
    bias_coeff_min <- bias_coeff_admix(admix_prop_bias_coeff_min, coanc_subpops)
    if (bias_coeff < bias_coeff_min)
        stop('Desired `bias_coeff` must be greater than ', bias_coeff_min, ' (the minimum achievable with `sigma = 0`), passed ', bias_coeff)
    
    # function whose zero we want!
    # "sigma" is the only parameter to optimize
    # need this function inside here as it depends on extra parameters
    bias_coeff_admix_objective <- function( x ) {
        # this is a transformation that helps us explore up to `sigma = Inf` in a compact space instead:
        sigma <- x / (1 - x)
        # x == 0 => sigma == 0
        # x == 1 => sigma == Inf
        # first, construct the admixture coefficients, which critically depend on sigma
        admix_proportions <- func(
            n_ind,
            k_subpops,
            sigma,
            coord_ind_first = coord_ind_first,
            coord_ind_last = coord_ind_last
        )
        # get bias coefficient, return difference from desired value!
        bias_coeff_admix(admix_proportions, coanc_subpops) - bias_coeff
    }
    
    # default tolerance of .Machine$double.eps^0.25 (~ 1e-4) was not good enough, reduced to ~ 2e-16
    obj <- stats::uniroot(
                      bias_coeff_admix_objective,
                      interval = c(0, 1),
                      extendInt = "no",
                      tol = .Machine$double.eps
                  )
    # this is the value in terms of `x`
    x_root <- obj$root
    # transform back to sigma, return this
    x_root / (1 - x_root)
}