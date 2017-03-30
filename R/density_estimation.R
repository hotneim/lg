
# Density estimation functions
# ----------------------------

#' Bivariate density estimation
#'
#' \code{dlg_bivariate} returns the locally Gaussian density estimate of a bivariate
#' distribution on a given grid.
#'
#' This function serves as the backbone in the body of methods concerning local
#' Gaussian correlation. It takes a bivariate data set, \code{x}, and a bivariate
#' set of grid points \code{eval_points}, and returns the bivariate, localy Gaussian
#' density estimate in these points. We also need a vector of bandwidths, \code{bw},
#' with two elements, and an estimation method \code{est_method}
#'
#' @param x The data matrix (or data frame). Must have exactly 2 columns
#' @param eval_points The grid where the density should be estimated. Must have exactly
#'   2 columns
#' @param grid_size If \code{eval_points} is not supplied, then the function will create
#'   a suitable grid diagonally through the data, with this many grid points
#' @param bw The two bandwidths, a numeric vector of length 2
#' @param est_method The estimation method, must either be "1par" for estimation with just
#'   the local correlation, or "5par"  for a full locally Gaussian fit with all 5 parameters
#' @param tol The numerical tolerance to be used in the optimization. Opnly applicable in
#' the 1-parameter optimization
#' @param run_checks Logical. Should sanity checks be run on the arguments? Useful to disable
#'   this when doing cross-validation for example.
#' @param marginal_estimates Provide the marginal estimates here if estimation method is 
#'   "\code{5par_marginals_fixed}", and the marginal estimates have already been found. Useful
#'   for cross-valdiation. List with two elements as returned by \code{dlg_marginal_wrapper}
#' @param bw_marginal Vector of bandwidths used to estimate the marginal distributions
#' @export
dlg_bivariate <- function(x,
                          eval_points = NA,
                          grid_size = 15,
                          bw = c(1, 1),
                          est_method = "1par",
                          tol = .Machine$double.eps^0.25/10^4,
                          run_checks = TRUE,
                          marginal_estimates = NA,
                          bw_marginal = NA) {

    # Check that everything is the way it should be -----------------------------
    if(run_checks) {
        x <- check_data(x, dim_check = 2, type = "data")
        check_bw_bivariate(bw = bw)
        check_est_method(est_method)
        if(!identical(eval_points, NA))
            eval_points <- check_data(eval_points, dim_check = 2, type = "grid")
        if(!identical(bw_marginal, NA))
             check_bw_bivariate(bw = bw_marginal)
    }

    # If eval_points not supplied, create a suitable grid
    if(identical(eval_points, NA)) {
        eval_points <- apply(x,
                             2,
                             function(x) seq(quantile(x,.001),quantile(x,.999),
                                             length.out = grid_size))
    }

    # Preliminary housekeeping - mostly so that we can re-use code from the previous
    # package that is working alright. We split the data matrix and the bandwidth vector
    # into their individual components:
    x1 <- x[,1]
    x2 <- x[,2]

    h1 <- bw[1]
    h2 <- bw[2]

    # OPTION 1: THE ONE-PARAMETER OPTIMIZATION

    if((est_method == "1par")) {

        # We declare a function that maximizes the local likelihood in one grid point.
        maximize_likelihood = function(grid_point) {
      
            x1_0 <- grid_point[1]
            x2_0 <- grid_point[2]
            
            # We need weights and some empirical moments in this grid point
            W <- dnorm(x1, mean = x1_0, sd = h1)*dnorm(x2, mean = x2_0, sd = h2)

            m1 <- mean(W)
            m2 <- mean(W*x1^2)
            m3 <- mean(W*x2^2)
            m4 <- mean(W*x1*x2)
            
            # The likelihood function
            lik <- function(rho) {
                - log(2*pi*sqrt(1 - rho^2))*m1 - m2/(2*(1 - rho^2)) - m3/(2*(1 - rho^2)) +
                    rho*m4/(1 - rho^2) - 1/2*exp(-1/2*(x2_0^2*h1^2 + x1_0^2 + x1_0^2*h2^2 -
                    2*x1_0*rho*x2_0 + x2_0^2)/(-rho^2 + h2^2 + 1 + h1^2 + h1^2*h2^2))/
                    (pi*(-rho^2 + h2^2 + 1 + h1^2 + h1^2*h2^2)^(1/2))
            }
            
            # Return the maximum of the likelihood and the density estimate
            opt <- try(optimise(lik,
                                lower = -1,
                                upper = 1,
                                maximum = TRUE,
                                tol = tol),
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

        ## Send the grid points to 'maximize_likelihood'
        est <- cbind(do.call(rbind,
                             lapply(X = split(eval_points, row(eval_points)),
                                    FUN = maximize_likelihood)))
        par_est <- matrix(est[,1])
        f_est = as.vector(est[,2])
        colnames(par_est) <- c('rho')
    }

    # OPTION 2: THE FIVE PARAMETER OPTIMIZATION

    if(est_method == "5par") {

        # Use the localgauss-package to find the local parameter estimates
        localgauss_result <- localgauss::localgauss(x = x1, y = x2,
                                                    b1 = h1, b2 = h2,
                                                    xy.mat = eval_points)

        # Collect the results
        par_est <- localgauss_result$par.est
        f_est <- dmvnorm_wrapper(eval_points = eval_points,
                                 mu_1 = par_est[, "mu_1"],
                                 mu_2 = par_est[, "mu_2"],
                                 sig_1 = par_est[, "sig_1"],
                                 sig_2 = par_est[, "sig_2"],
                                 rho = par_est[, "rho"],
                                 run_checks = run_checks)
    }

    # OPTION 3: FIVE PARAMETERS WITH MARGINALS FIXED
    if(est_method == "5par_marginals_fixed") {

        # Have the marginals already been estimated? If not, do it
        if(identical(marginal_estimates, NA)) {
            marginal_estimates <- dlg_marginal_wrapper(data_matrix = x,
                                                       eval_matrix = eval_points,
                                                       bw_vector = bw_marginal)
        }

        # We declare a function that maximizes the local likelihood in one grid point,
        # with fixed marginal parameter estimates.
        maximize_likelihood_fixed_marginals =
            function(grid_point_and_marginal_parameters) {
      
            x1_0 <- grid_point_and_marginal_parameters[1]
            x2_0 <- grid_point_and_marginal_parameters[2]
            mu_1 <- grid_point_and_marginal_parameters[3]
            sig_1 <- grid_point_and_marginal_parameters[4]
            mu_2 <- grid_point_and_marginal_parameters[5]
            sig_2 <- grid_point_and_marginal_parameters[6]
            
            # We need weights and some empirical moments in this grid point
            W <- dnorm(x1, mean = x1_0, sd = h1)*dnorm(x2, mean = x2_0, sd = h2)

            m1 <- mean(W)
            m2 <- mean(W*((x1-mu_1)/sig_1)^2)
            m3 <- mean(W*((x2-mu_2)/sig_2)^2)
            m4 <- mean(W*((x1-mu_1)/sig_1)*((x2-mu_2)/sig_2))
            
            # The likelihood function
            lik <- function(rho) {
                - log(2*pi*sqrt(1 - rho^2)*sig_1*sig_2)*m1 - m2/(2*(1 - rho^2)) - m3/(2*(1 - rho^2)) +
                    rho*m4/(1 - rho^2) -
                        mvtnorm::dmvnorm(c(x1_0, x2_0),
                                         mean = c(mu_1, mu_2),
                                         sigma = matrix(c(sig_1^2 + h1^2,
                                                          sig_1*sig_2*rho,
                                                          sig_1*sig_2*rho,
                                                          sig_2^2 + h2^2),
                                                        byrow = TRUE,
                                                        ncol = 2))
            }
            
            # Return the maximum of the likelihood and the density estimate
            opt <- try(optimise(lik,
                                lower = -1,
                                upper = 1,
                                maximum = TRUE,
                                tol = tol),
                       silent = TRUE)

            # Store the result if the optimization went alright. Return NA if not. 
            if(class(opt) != "try-error") {
                return(c(opt$maximum,
                         mvtnorm::dmvnorm(c(x1_0, x2_0), mean = c(mu_1,mu_2),
                                          sigma = matrix(c(sig_1^2,
                                                           sig_1*sig_2*opt$maximum,
                                                           sig_1*sig_2*opt$maximum,
                                                           sig_2^2), 2))))
            } else {
                return(c(NA, NA))
            }
        }

        ## Send the grid points to 'maximize_likelihood_fixed_marginals'
        est <- cbind(do.call(rbind,
                             lapply(X = split(cbind(eval_points,
                                                    marginal_estimates[[1]], marginal_estimates[[2]]),
                                              row(eval_points)),
                                    FUN = maximize_likelihood_fixed_marginals)))
        par_est <- cbind(marginal_estimates[[1]][,1],
                         marginal_estimates[[2]][,1],
                         marginal_estimates[[1]][,2],
                         marginal_estimates[[2]][,2],
                         matrix(est[,1]))
        f_est = as.vector(est[,2])
        colnames(par_est) <- c("mu_1", "mu_2", "sig_1", "sig_2", "rho")

        
    }
    
    return(list(x = x,
                eval_points = eval_points,
                bw = bw,
                par_est = par_est,
                f_est = f_est))
       
}

#' Marginal density estimation
#'
#' Function that estimates a univariate density estimation by local Gaussian
#' approximations
#'
#' This functio is mainly mean to be used as a tool in multivariate analysis
#' as away to obtain the estimate of a univariate (marginal) density function,
#' but it can of course be used in general to estimate univariate densities.
#'
#' @param x The data vector.
#' @param bw The bandwidth (a single number).
#' @param grid_size Number of grid points if grid is not provided.
#' @param eval_points The grid where we want to evaluate the density. Chosen
#'   suitably if not provided, with length equal to grid_size. 
#' @export
dlg_marginal <- function(x,
                        bw = 1,
                        eval_points = seq(quantile(x, .01), quantile(x, .99), length.out = grid_size),
                        grid_size = 15) {

    # Function that maximizes the local Gaussian likelihood in one grid point
    maximize_likelihood_univariate <- function(eval_point) {

        W <- dnorm(x, mean = eval_point, sd = bw)

        m1 <- mean(W)
        m2 <- mean(W*x)
        m3 <- mean(W*(x^2))

        # The negative likelihood
        lik <- function(par) {
            mu <- par[1]
            sig <- par[2]
            - m2*mu/sig^2 + m3/(2*sig^2) + m1*(mu^2/(2*sig^2) + .5*log(2*pi*sig^2)) +
                dnorm((eval_point - mu)/sqrt(sig^2 + bw^2))/sqrt(sig^2 + bw^2) 
        }

        # Do the optimization
        opt <- try(optim(c(0, 1),
                         lik))

        # Return the result
        if(class(opt) != "try-error") {
            return(c(opt$par,
                     dnorm(eval_point,
                           mean = opt$par[1],
                           sd = opt$par[2])))
        } else {
            return(c(NA, NA))
        }        
    }

    # Send the grid points to 'maximize_likelihood_univariate'
    est <- cbind(do.call(rbind,
                         lapply(X = as.list(eval_points),
                                FUN = maximize_likelihood_univariate)))

    # Collect and return the results
    if(nrow(est) == 1) {
        par_est <- matrix(est[, 1:2], ncol = 2)
    } else {
        par_est <- est[, 1:2]
    }
    
    colnames(par_est) <- c("mu", "sig")
    f_est <- as.vector(est[, 3])

    return(list(x = x,
                eval_points = eval_points,
                bw = bw,
                par_est = par_est,
                f_est = f_est))
}

#' Marginal estimates for multivariate data
#'
#' Estimates the marginal locally Gaussian parameter estimates for a multivariate
#' data set
#'
#' This function takes in a matrix of observations, a matrix of evaluation points
#' and a vector of bandwidths, and does the marginal estimation pairwise. This
#' function assumes that the data and evaluation points are organized column-wise
#' in matrices, and that the bandwidth is found in the corresponding element in
#' the bandwidth matrix.
#'
#' @param data_matrix The matrix of data points. One column constitutes an
#'   observation vector.
#' @param eval_matrix The matrix of evaluation points. One column constitutes
#'   a vector of grid points.
#' @param bw_vector The vector of bandwidths, one element per component.
dlg_marginal_wrapper <- function(data_matrix, eval_matrix, bw_vector){

    # Return the parameters in a list
    ret <- list()

    # estimate the parameters for each marginal component
    for(i in 1:ncol(data_matrix)) {
        ret[[i]] <- dlg_marginal(x = data_matrix[,i],
                                 bw = bw_vector[i],
                                 eval_points = eval_matrix[,i])$par_est
    }

    return(ret)
}
