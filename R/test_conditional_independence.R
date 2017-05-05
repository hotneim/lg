
# Functions for testing for conditional independence
# --------------------------------------------------

#' Interpolate a univariate conditional density function
#'
#' Estimates the conditional density function for one free variable on
#' a grid. Returns a function that interpolates between these grid points
#' so that it can be evaluated more quickly, without new optimizations.
#'
#' @param lg_object An object of type \code{lg}, as produced by the \code{lg}-function
#' @param condition A vector with conditions for the variables that we condition upon.
#'   Must have exactly one more element than there are columns in the data
#' @param nodes Either the number of equidistant nodes to generate, or a vector of nodes
#'   supplied by the user
#' @param extend How far to extend the grid beyond the extreme data points, in share of
#'   the range
#' @export
interpolate_conditional_density <- function(lg_object,
                                            condition,
                                            nodes,
                                            extend = .3) {

    # Extract some basic info
    n  <- nrow(lg_object$x)        # Sample size
    d  <- ncol(lg_object$x)        # Number of variables
    nc <- length(condition)        # Number of conditioning variables
    
    # Do some checks
    check_lg(lg_object)
    if((d - nc) != 1)
        stop("The number of free variables must be exactly 1")

    # Specify the nodes if not supplied
    # Create the grid
    if(length(nodes) == 1) {
        width <- diff(range(lg_object$x[,1]))
        grid <- matrix(seq(from = min(lg_object$x[,1]) - extend*width,
                           to   = max(lg_object$x[,1]) + extend*width,
                           length.out = nodes),
                       ncol = 1)
        delta <- diff(grid)[1]
    } else {
        grid <- matrix(nodes,
                       ncol = 1)
    }

    # Estimate the conditional density
    conditional_density_estimate <- clg(lg_object, grid = grid, condition = condition)

    # Create the interpolation function
    x <- c(grid)
    y <- conditional_density_estimate$f_est

    conditional_density <- function(t) {
        exp(akima::aspline(x[!(log(y) == Inf | log(y) == -Inf)],
                           log(y)[!(log(y) == Inf | log(y) == -Inf)],
                           t,
                           method = "improved",
                           degree = 2)$y)
    }

    return(conditional_density)
}


#' Test for conditional independence
#'
#' Test for conditional independence using local Gaussian correlations.
#'
#' This is the main function for performing a test for conditional independence
#' between two variables given a (set of) variable(s). The function takes in an
#' lg-object as produced by the main \code{lg}-function, and always assumes that
#' we want to test conditional independence between th variables represented by
#' the first two columns, given the rest.
#'
ci_test <- function(lg_object,
                    n_replicates <- 1000,
                    h = function(x) x^2,
                    return_replicated_test_statistics = TRUE,
                    return_replicated_data = FALSE,
                    return_auxiliary_info = FALSE) {
}


#' Generate replicates under the null hypothesis for conditional independence test
#'
generate_replicates <- function(x,
                                n_replicates,
                                h) {
}
