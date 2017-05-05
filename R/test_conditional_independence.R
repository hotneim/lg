
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

    # Mean and variance of the density estimate, useful for accept-reject
    if(length(nodes) == 1) {
        mean <- sum(x*y)*delta
        var <- sum((x - mean)^2*y)*delta
    } else {
        mean = NULL
        var = NULL
    }

    conditional_density <- function(t) {
        exp(akima::aspline(x[!(log(y) == Inf | log(y) == -Inf)],
                           log(y)[!(log(y) == Inf | log(y) == -Inf)],
                           t,
                           method = "improved",
                           degree = 2)$y)
    }

    return(list(conditional_density = conditional_density,
                mean = mean,
                var = var))
}

#' Generate sample from a conditional density estimate
#'
#' Generate a sample from a locally Gaussian conditional density estimate using
#' the accept-reject algorithm.
#'
#' @param lg_object An object of type \code{lg}, as produced by the \code{lg}-function
#' @param n_new The number of observations to generate
#' @param nodes Either the number of equidistant nodes to generate, or a vector of nodes
#'   supplied by the user
#' @param M The value for M in the accept-reject algorithm if already known
#' @param M_sim The number of replicates to simulate in order to find a value for M
#' @param M_corr Correction factor for M, to be on the safe side
#' @param n_corr Correction factor for n_new, so that we mostly will generate enough
#'   observations in the first go
#' @param return_just_M \code{TRUE} if we just want to find M, without actually generating
#'   any replications.
#' @param extend How far to extend the grid beyond the extreme data points when interpolating,
#'   in share of the range
#' @export
accept_reject <- function(lg_object,
                          condition,
                          n_new,
                          nodes,
                          M = NULL,
                          M_sim = 1500,
                          M_corr = 1.5,
                          n_corr = 1.2,
                          return_just_M = FALSE,
                          extend = .3) {

    # Get function for interpolation
    int_object <- interpolate_conditional_density(lg_object, condition, nodes, extend)
    
    # Do a quick check
    if(is.null(M) && is.null(int_object$mean))
        stop("Not enough information to calculate M")

    # Calculate M for use in the accept-reject algorithm
    if(is.null(M)){
        test_grid <- rnorm(M_sim,
                           mean = int_object$mean,
                           sd = sqrt(int_object$var))
        test_f <- int_object$conditional_density(test_grid)
        test_g <- dnorm(test_grid,
                        mean = int_object$mean,
                        sd = sqrt(int_object$var))
        M <- min(10,  M_corr*quantile(test_f/test_g, .95, type = 1))
    }

    # Function for generating a sample from the candidate, and retaining the accepts
    get <- function(k) {
        # Generate k samples from the candidate and the uniform
        Y <- rnorm(k,
                   mean = int_object$mean,
                   sd = sqrt(int_object$var))
        U <- runif(k)
        g <- dnorm(Y,
                   mean = int_object$mean,
                   sd = sqrt(int_object$var))
        f <- int_object$conditional_density(Y)
        Y[U < (f/(M*g))]
    }

    # Get the replicates, or just return M if that was specified
    if(!return_just_M) {
        generated <- get(round(n_corr*M*n_new))

        # Do we have enough samples?
        enough = length(generated) >= n_new
 
        # Get more samples if needed
        while(!enough) {
            n_missing <- n_new - length(generated)
            generated <- c(generated, get(round(n_corr*M*n_missing)))
            enough = length(generated) >= n_new        
        }

        return(list(sample = generated[1:n_new], M = M))

    } else {

        return(list(M = M))
    }
}

#' Test for conditional independence
#'
#' Test for conditional independence using local Gaussian correlations.
#'
#' This is the main function for performing a test for conditional independence
#' between two variables given a (set of) variable(s). The function takes in an
#' lg-object as produced by the main \code{lg}-function, 
#'
ci_test <- function(lg_object,
                    n_replicates = 1000,
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
