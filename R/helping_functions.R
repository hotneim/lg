
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
#' \code{dmvnorm_wrapper} is a function that evaluates the bivariate normal
#' distribution in a matrix of evaluation points, with local parameters.
#'
#' This functions takes as arguments a matrix of grid points, and vectors of
#' parameter values, and returns the bivariate normal density at these points,
#' with these parameter values.
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

    # Collect the arguments in one matrix, so that we can apply the single
    # wrapper function
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

#' Evaluate the multivariate normal
#'
#' Function that evaluates the multivariate normal distribution with local
#' parameters
#'
#' Takes in a grid, where we want to evaluate the multivariate normal, and in
#' each grid point we have a new set of parameters.
#'
#' @param eval_points A matrix of grid points
#' @param loc_mean A matrix of local means, one row per grid point, one column
#'   per component
#' @param loc_sd A  matrix of local standard deviations, one row per grid point,
#'   one column per component
#' @param loc_cor A matrix of local correlations, one row per grid point, on
#'   column per pair of variables
#' @param pairs A data frame specifying the components that make up each pair,
mvnorm_eval <- function(eval_points,
                        loc_mean,
                        loc_sd,
                        loc_cor,
                        pairs) {

    # Evaluation in one point
    single_eval <- function(i) {
        mu_vec <- loc_mean[i,]
        sigma <- diag(loc_sd[i,]^2)
        for(j in 1:nrow(pairs)) {
            var1 <- pairs$x1[j]
            var2 <- pairs$x2[j]
            rho <- loc_cor[i, j]
            sigma[var1, var2] <- rho*sqrt(sigma[var1, var1]*sigma[var2, var2])
            sigma[var2, var1] <- sigma[var1, var2]
        }
        mvtnorm::dmvnorm(eval_points[i,],
                        mean = mu_vec,
                        sigma = sigma)
    }

    unlist(lapply(X = as.list(1:nrow(eval_points)),
                  FUN = single_eval))
}

#' Auxiliary function for calculating the local score function u
#'
#' Auxiliary function for calculating the local score function u
#'
#' This function is used to estimate the asymptotic variance of the estimates.
#'
#' @param z1 z1
#' @param z2 z2
#' @param rho rho

u <- function(z1, z2, rho) {
  (rho^3 - z1*z2*(1 + rho^2) + (z1^2 + z2^2 - 1)*rho)/((1 - rho^2)^2)
}

#' Auxiliary function for calculating the asymptotic standard deviations for the
#' local Gaussian correlations
#'
#' Auxiliary function for calculating the asymptotic standard deviations for the
#' local Gaussian correlations
#'
#' @param r r
#' @param pairs pairs
#' @param p p
make_C <- function(r, pairs, p) {

  R <- diag(p)

  for(j in 1:nrow(pairs)) { # For each pair of variables

    var1 <- pairs[j,1]      #| The variables
    var2 <- pairs[j,2]      #|

    R[var1, var2] <- R[var2, var1] <- r[j]
  }

  R22 <- R[3:ncol(R), 3:ncol(R)]
  R21 <- R[3:ncol(R), 1:2]

  solve(R22) %*% R21

}

#' Auxiliary function for calculating the asymptotic standard deviations for the
#' local Gaussian correlations
#'
#' Auxiliary function for calculating the asymptotic standard deviations for the
#' local Gaussian correlations
#'
#' @param sigma sigma
#' @param sigma_k sigma_k
gradient <- function(sigma, sigma_k) {
  (sigma[1,1]*sigma[2,2]*sigma_k[1,2] - sigma[1,2]*(sigma_k[1,1]*sigma[2,2] + sigma[1,1]*sigma_k[2,2]))/
    (sigma[1,1]*sigma[2,2])^2
}
