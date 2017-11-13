
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
#'   two colimns names 'x1' and 'x2', one row per pair of components
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



#' Helping function for concise looping in CV bandwidth selection
#'
#' Function that is used in order to accomodate parellell computing in bandwidth
#' selection. This function is only used internally.
#'
#' @param i Indicate pair number
#' @param x The bivariate data set
#' @param joint_bandwidths Pairwise bandwidths
#' @param tol_joint Tolerance in the cross validation of pairwise bandwidths
#' @param est_method The estimation method, must either be "1par" for estimation
#'   with just the local correlation, or "5par"  for a full locally Gaussian fit
#'   with all 5 parameters
#' @param marginal_bandwidths The marginal bandwidths, used in the optimization
bandwidth_selection_cv_loop_helpfunc <- function(i,
                                                 x,
                                                 joint_bandwidths,
                                                 tol_joint,
                                                 est_method,
                                                 marginal_bandwidths) {

  variables <- c(joint_bandwidths$x1[i], joint_bandwidths$x2[i])

  # Extract the pairs of variables
  bivariate_data <- x[, variables]


  result <- bw_select_cv_bivariate(x = bivariate_data,
                                   tol = tol_joint,
                                   est_method = est_method,
                                   bw_marginal = marginal_bandwidths[variables])

  if(result$convergence != 0)
    warning(paste("Cross valdidation for joint bandwidths",
                  as.character(variables[1]), "and",  as.character(variables[2]),
                  "did not converge properly"))

  return(c(result$bw[1],result$bw[2],result$convergence))
}

#' Help function for concise looping in sequential bivariate density estimation
#'
#' Function that is used in accomodating paralell computing of locally Gaussian
#' density estimation. This function is only used internally.
#'
#' @param i Indicate pair number
#' @param x The bivariate data set
#' @param x0 Evaluation points
#' @param lg_object The \code{lg}-object calculated for this data set
#' @param marginal_estimates The marginal density estimates
#' @param pairs Matrix indicating the pairs of variables.
dlg_bivariate_loop_helpfunc <- function(i,
                                        x,
                                        x0,
                                        lg_object,
                                        marginal_estimates,
                                        pairs){

  if(lg_object$est_method == "1par") {
    pairwise_marginal_estimates <- NA
  } else {
    pairwise_marginal_estimates <- list(marginal_estimates[[pairs$x1[i]]],
                                        marginal_estimates[[pairs$x2[i]]])
  }

  pairwise_estimate <-
    dlg_bivariate(x = x[, c(pairs$x1[i], pairs$x2[i])],
                  eval_points = x0[, c(pairs$x1[i], pairs$x2[i])],
                  bw = c(lg_object$bw$joint[i, "bw1"],
                         lg_object$bw$joint[i, "bw2"]),
                  est_method = lg_object$est_method,
                  marginal_estimates = pairwise_marginal_estimates)

  return(pairwise_estimate$par_est[, "rho"])
}

#' 1-parameter likelihood function for bivariate data
#'
#' Used for optimization, only used internally.
#'
#' @param rho The local Gaussian correlation, with respect to which we optimize
#' @param m1 Empirical mean 1, used in the optimization
#' @param m2 Empirical mean 2, used in the optimization
#' @param m3 Empirical mean 3, used in the optimization
#' @param m4 Empirical mean 4, used in the optimization
#' @param x1_0 Grid point 1
#' @param x2_0 Grid point 2
#' @param h1 Bandwidth 1
#' @param h2 Bandwidth 2
#'
#' @return The likelihood for the one-parameter case.

lik_1par <- function(rho,m1,m2,m3,m4,x1_0,x2_0,h1,h2) {
  - log(2*pi*sqrt(1 - rho^2))*m1 - m2/(2*(1 - rho^2)) - m3/(2*(1 - rho^2)) +
    rho*m4/(1 - rho^2) - 1/2*exp(-1/2*(x2_0^2*h1^2 + x1_0^2 + x1_0^2*h2^2 -
                                         2*x1_0*rho*x2_0 + x2_0^2)/(-rho^2 + h2^2 + 1 + h1^2 + h1^2*h2^2))/
    (pi*(-rho^2 + h2^2 + 1 + h1^2 + h1^2*h2^2)^(1/2))
}



#' 1-parameter likelihood maximization function for bivariate data
#'
#' Optimization of the local likelihood function.
#'
#' @param grid_point Grid point
#' @param x1 The first data vector
#' @param x2 The second data vector
#' @param h1 The first bandwidth
#' @param h2 The second bandwidth
#' @param tol The tolerance used by the \code{optim}-function
#'
#' @return The maximized likelihood for the one-parameter case.
#'
#' @export

maximize_likelihood_1par = function(grid_point,
                                    x1, x2,
                                    h1, h2,
                                    tol) {

  x1_0 <- grid_point[1]
  x2_0 <- grid_point[2]

  # We need weights and some empirical moments in this grid point
  W <- dnorm(x1, mean = x1_0, sd = h1)*dnorm(x2, mean = x2_0, sd = h2)

  m1 <- mean(W)
  m2 <- mean(W*x1^2)
  m3 <- mean(W*x2^2)
  m4 <- mean(W*x1*x2)


  # Return the maximum of the likelihood and the density estimate
  opt <- try(optimise(lik_1par,
                      lower = -1,
                      upper = 1,
                      maximum = TRUE,
                      tol = tol,
                      m1 = m1,
                      m2 = m2,
                      m3 = m3,
                      m4 = m4,
                      x1_0 = x1_0,
                      x2_0 = x2_0,
                      h1 = h1,
                      h2 = h2),
             silent = TRUE)

  # Store the result if the optimization went alright. Return NA if not.
  if(class(opt) != "try-error") {
    return(c(opt$maximum,
             mvtnorm::dmvnorm(c(x1_0, x2_0), mean = c(0,0),
                              sigma = matrix(c(1, opt$maximum, opt$maximum, 1), 2))))
  } else {
    return(c(NA, NA))
  }
}


