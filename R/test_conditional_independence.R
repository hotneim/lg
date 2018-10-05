
# Functions for testing for conditional independence
# --------------------------------------------------

#' Interpolate a univariate conditional density function
#'
#' Estimates the conditional density function for one free variable on a grid.
#' Returns a function that interpolates between these grid points so that it can
#' be evaluated more quickly, without new optimizations.
#'
#' @param lg_object An object of type \code{lg}, as produced by the
#'   \code{lg_main}-function
#' @param condition A vector with conditions for the variables that we condition
#'   upon. Must have exactly one more element than there are columns in the data
#' @param nodes Either the number of equidistant nodes to generate, or a vector
#'   of nodes supplied by the user
#' @param extend How far to extend the grid beyond the extreme data points, in
#'   share of the range
#' @param gaussian_scale Stay on the standard Gaussian scale, useful for the
#'   accept-reject algorithm
interpolate_conditional_density <- function(lg_object,
                                            condition,
                                            nodes,
                                            extend = .3,
                                            gaussian_scale = lg_object$transform_to_marginal_normality) {

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
    } else {
        grid <- matrix(nodes,
                       ncol = 1)
    }

    # Estimate the conditional density
    conditional_density_estimate <- clg(lg_object, grid = grid, condition = condition)

    # Create the interpolation function
    if(gaussian_scale) {
        x <- c(conditional_density_estimate$transformed_grid[,1])
        y <- conditional_density_estimate$f_est/conditional_density_estimate$normalizing_constants[,1]
    } else {
        x <- c(grid)
        y <- conditional_density_estimate$f_est
    }

    # Mean and variance of the density estimate, useful for accept-reject
    if(length(nodes) == 1) {
        delta <- diff(x)
        mean <- sum(x[-1]*y[-1]*delta, na.rm = TRUE)
        var <- sum((x[-1] - mean)^2*y[-1]*delta, na.rm = TRUE)
    } else {
        mean = NULL
        var = NULL
    }

    conditional_density <- function(t) {
      exp(stats::splinefun(x = x[!(log(y) == Inf | log(y) == -Inf)],
                           y = log(y)[!(log(y) == Inf | log(y) == -Inf)],
                           method = "natural")(t))
    }

    return(list(conditional_density = conditional_density,
                mean = mean,
                var = var))
}

#' Generate sample from a conditional density estimate
#'
#' Generate a sample from a locally Gaussian conditional density estimate using
#' the accept-reject algorithm. If the \code{transform_to_marginal_normality}-
#' component of the lg_object is \code{TRUE}, the replicates will be on the
#' standard normal scale.
#'
#' @param lg_object An object of type \code{lg}, as produced by the
#'   \code{lg_main}-function
#' @param condition The value of the conditioning variables
#' @param n_new The number of observations to generate
#' @param nodes Either the number of equidistant nodes to generate, or a vector
#'   of nodes supplied by the user
#' @param M The value for M in the accept-reject algorithm if already known
#' @param M_sim The number of replicates to simulate in order to find a value
#'   for M
#' @param M_corr Correction factor for M, to be on the safe side
#' @param n_corr Correction factor for n_new, so that we mostly will generate
#'   enough observations in the first go
#' @param return_just_M \code{TRUE} if we just want to find M, without actually
#'   generating any replications.
#' @param extend How far to extend the grid beyond the extreme data points when
#'   interpolating, in share of the range
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
        test_grid <- stats::rnorm(M_sim,
                           mean = int_object$mean,
                           sd = sqrt(int_object$var))
        test_f <- int_object$conditional_density(test_grid)
        test_g <- dnorm(test_grid,
                        mean = int_object$mean,
                        sd = sqrt(int_object$var))
        M <- min(10,  M_corr*quantile(test_f/test_g, .95, type = 1))
    }

    # For some very extreme observations, M can become very small here. For the time being,
    # we will return nothing if that happens
    if(M < 1) {
        return(0)
    } else {
        # Function for generating a sample from the candidate, and retaining the accepts
        get <- function(k) {
            # Generate k samples from the candidate and the uniform
            Y <- stats::rnorm(k,
                       mean = int_object$mean,
                       sd = sqrt(int_object$var))
            U <- stats::runif(k)
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
}

#' Bootstrap replication under the null hypothesis
#'
#' Generate bootstrap replicates under the null hypothesis that the first two
#' variables are conditionally independent given the rest of the variables.
#'
#' @param lg_object An object of type \code{lg}, as produced by the
#'   \code{lg_main}-function
#' @param n_rep  The number of replicated bootstrap samples
#' @param nodes Either the number of equidistant nodes to generate, or a vector
#'   of nodes supplied by the user
#' @param M The value for M in the accept-reject algorithm if already known
#' @param M_sim The number of replicates to simulate in order to find a value
#'   for M
#' @param M_corr Correction factor for M, to be on the safe side
#' @param n_corr Correction factor for n_new, so that we mostly will generate
#'   enough observations in the first go
#' @param extend How far to extend the grid beyond the extreme data points when
#'   interpolating, in share of the range
replicate_under_ci <- function(lg_object,
                               n_rep,
                               nodes,
                               M = NULL,
                               M_sim = 1500,
                               M_corr = 1.5,
                               n_corr = 1.2,
                               extend = .3) {

    # Extract some info
    n <- nrow(lg_object$x)
    d <- ncol(lg_object$x)

    # Initialize a big list of replicates. Fill the matrices of X1 and X2 by row, pick them
    # apart by column. Also, we want to return the Ms that were used for accept reject.
    X <- list(X1 = matrix(NA, nrow = n, ncol = n_rep),
              X2 = matrix(NA, nrow = n, ncol = n_rep))
    Mret <- matrix(NA, nrow = n, ncol = 2)

    # Go over two times, one for each variable
    for(i in 1:2) {

        # Data without the other variable
        temp_lg_object <- lg_main(x = lg_object$x[,-i],
                                  bw_method = lg_object$bw_method,
                                  est_method = lg_object$est_method,
                                  transform_to_marginal_normality = lg_object$transform_to_marginal_normality,
                                  plugin_constant_marginal = lg_object$plugin_constant_marginal,
                                  plugin_constant_joint = lg_object$plugin_constant_joint,
                                  plugin_exponent_marginal = lg_object$plugin_exponent_marginal,
                                  plugin_exponent_joint = lg_object$plugin_exponent_joint,
                                  tol_marginal = lg_object$tol_marginal,
                                  tol_joint = lg_object$tol_joint)

        # Fill each row with a new sample, one for each condition
        for(j in 1:n) {
            new_sample <- accept_reject(temp_lg_object,
                                        condition = temp_lg_object$x[j, 2:(d-1)],
                                        n_new = n_rep,
                                        nodes = nodes,
                                        M = NULL,
                                        M_sim = M_sim,
                                        M_corr = M_corr,
                                        n_corr = n_corr,
                                        extend = extend)
            if(!is.numeric(new_sample)) {
                X[[i]][j,] <- new_sample$sample
                Mret[j,i] <- new_sample$M
            } else {
                X[[i]][j,] <- NA
                Mret[j,i] <- NA
            }
        }
    }

    # Distribute the generated observations into a list of replicated data
    replicates <- list()
    for(i in 1:n_rep) {
        replicates[[i]] <- lg_object$transformed_data
        replicates[[i]][, 1] <- X[[2]][,i]
        replicates[[i]][, 2] <- X[[1]][,i]
    }

    return(list(replicates = replicates,
                Mret = Mret))
}

#' Calculate the local conditional covariance between two variables
#'
#' Wrapper for the \code{clg} function that extracts the local Gaussian conditional
#' covariance between two variables from an object that is produced by the clg-function.
#'
#' This function is a wrapper for the clag-function, and extracts the estimated local
#' conditional covariance between the first two variables in the data matrix, on the grid
#' specified to the clg-function.
#'
#' @param clg_object The object produced by the clg-function
#' @param coord The variables for which the conditional covariance should be extracted
local_conditional_covariance <- function(clg_object, coord = c(1, 2)) {

    unlist(lapply(clg_object$c_cov, `[`, coord[1], coord[2]))

}

#' Calculate the value of the test statistic for the conditional independence
#' test
#'
#' Calculate the test statistic in the test for conditional independence between
#' the first two variables in the data set, given the remaining variables.
#'
#' @param lg_object An object of type \code{lg}, as produced by the
#'   \code{lg_main}-function
#' @param h The \code{h}-function used in the calculation of the test statistic.
#'   The default value is \code{h(x) = x^2}.
ci_test_statistic <- function(lg_object, h = function(x) x^2) {

    # Calculate the conditional coveriance between the first two variables in the data points
    clg_object <- clg(lg_object, fixed_grid = lg_object$x)

    # Extract the conditional covariances and calculate the value of the test statistic
    mean(h(local_conditional_covariance(clg_object)))
}

#' Test for conditional independence
#'
#' Perform a test for conditional independence between the first two variables
#' in the data set, given the remaining variables.
#'
#' @param lg_object An object of type \code{lg}, as produced by the
#'   \code{lg_main}-function
#' @param h The \code{h}-function used in the calculation of the test statistic.
#'   The default value is \code{h(x) = x^2}.
#' @param n_rep The number of replicated bootstrap samples
#' @param nodes Either the number of equidistant nodes to generate, or a vector
#'   of nodes supplied by the user
#' @param M The value for M in the accept-reject algorithm if already known
#' @param M_sim The number of replicates to simulate in order to find a value
#'   for M
#' @param M_corr Correction factor for M, to be on the safe side
#' @param n_corr Correction factor for n_new, so that we mostly will generate
#'   enough observations in the first go
#' @param extend How far to extend the grid beyond the extreme data points when
#'   interpolating, in share of the range
#' @param return_time Measure how long the test takes to run, and return along
#'   with the test result
#' @export
ci_test <- function(lg_object, h = function(x) x^2, n_rep = 1000, nodes = 1000, M = NULL,
                    M_sim = 1500, M_corr = 1.5, n_corr = 1.2, extend = .3, return_time = TRUE) {

    ## Start the clock
    if(return_time) {
        start_time <- Sys.time()
    }

    ## Generate bootstrap samples under the null
    replicates <- replicate_under_ci(lg_object, n_rep = n_rep, nodes = nodes, M = M, M_sim = M_sim,
                                     M_corr = M_corr, n_corr = n_corr, extend = extend)

    ## The observed test functional
    observed <- ci_test_statistic(lg_object, h = h)

    ## Initialize the vector where we want to store the replicated test functionals
    replicated <- rep(NA, n_rep)

    ## Run over the replicated data sets in order to calculate the test statistics
    ## under the null hypothesis
    for(i in 1:n_rep) {

        # First, create the temporary lg-object. Parameters are the same as for the original data.
        temp_lg_object <- lg_main(replicates$replicates[[i]], bw_method = lg_object$bw_method,
                                  est_method = lg_object$est_method,
                                  transform_to_marginal_normality = lg_object$transform_to_marginal_normality,
                                  plugin_constant_marginal = lg_object$plugin_constant_marginal,
                                  plugin_constant_joint = lg_object$plugin_constant_joint,
                                  plugin_exponent_marginal = lg_object$plugin_exponent_marginal,
                                  plugin_exponent_joint = lg_object$plugin_exponent_joint,
                                  tol_marginal = lg_object$tol_marginal,
                                  tol_joint = lg_object$tol_joint)

        # Then calculate the test statistic
        replicated[i] <- ci_test_statistic(temp_lg_object, h = h)

    }

    ## Stop the clock
    if(return_time) {
        end_time <- Sys.time()
        duration <- as.numeric(end_time - start_time, unit = "mins")
    } else {
        duration <- NA
    }

    return(list(p_value = mean(observed < replicated),
                observed = observed,
                replicated = replicated,
                Mret = replicates$Mret,
                duration = duration))
}
