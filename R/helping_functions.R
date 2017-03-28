
# Various helping functions
# -------------------------

#' Wrapper for \code{dmvnorm} - single point
#'
#' Function that evaluates the bivariate normal in a single point
#'
#' @param x1 The first component of the evaluation point
#' @param x2 The second component of the evaluation point
#' @param mu_1 The first expectation
#' @param mu_2 The second expectation
#' @param sig_1 The first standard deviation
#' @param sig_2 The second standard deviation
#' @param rho The correlation
dmvnorm_wrapper_single <- function(x1, x2, mu_1, mu_2, sig_1, sig_2, rho) {
    mvtnorm::dmvnorm(x = c(x1, x2),
                     mean = c(mu_1, mu_2),
                     sigma = matrix(c(sig_1^2, rho*sig_1*sig_2, rho*sig_1*sig_2, sig_2^2),
                                    ncol = 2,
                                    byrow = TRUE))
}

#' Wrapper for \code{dmvnorm}
#'
#' \code{dmvnorm_wrapper} is a function that evaluates the bivariate normal distribution
#' in a matrix of evaluation points, with local parameters.
#'
#' This functions takes as arguments a matrix of grid points, and vectors of parameter values,
#' and returns the bivariate normal density at these points, with these parameter values.
#'
#' @param eval_points A \code{kx2} matrix with evaluation points
#' @param mu_1 The first expectation vector
#' @param mu_2 The second expectation vector
#' @param sig_1 The first standard deviation vector
#' @param sig_2 The second standard deviation vector
#' @param rho The correlation vector
#' @param run_checks Run sanity check for the arguments
dmvnorm_wrapper <- function(eval_points,
                            mu_1  = rep(0, nrow(eval_points)),
                            mu_2  = rep(0, nrow(eval_points)),
                            sig_1 = rep(1, nrow(eval_points)),
                            sig_2 = rep(1, nrow(eval_points)),
                            rho   = rep(0, nrow(eval_points)),
                            run_checks = TRUE) {

    # Check that the arguments have the correct format
    if(run_checks) {
        check_dmvnorm_arguments(eval_points, mu_1, mu_2, sig_1, sig_2, rho)
    }
    
    # Collect the arguments in one matrix, so that we can apply the single wrapper function
    arguments <- cbind(eval_points, mu_1, mu_2, sig_1, sig_2, rho)

    # Calculate the bivariate normal
    apply(arguments, 1, function(x) dmvnorm_wrapper_single(x1 = x[1],
                                                           x2 = x[2],
                                                           mu_1 = x[3],
                                                           mu_2 = x[4],
                                                           sig_1 = x[5],
                                                           sig_2 = x[6],
                                                           rho = x[7]))                        
                    
}
