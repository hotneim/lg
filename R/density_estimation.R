
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

#' The locally Gaussian density estimator (LGDE)
#'
#' Estimate a mulivariate density function using locally Gaussian approximations
#'
#' This function does multivariate density estimation using the locally Gaussian
#' density estimator (LGDE), that was introduced by Otneim & Tjøstheim (2016)
#' [Proper citation needed!]. The function takes as arguments an \code{lg}-object
#' as produced by the main \code{lg}-function, and a grid of points where the density
#' estimate should be estimated.
#' @param lg_object An object of type \code{lg}, as produced by the \code{lg}-function
#' @param grid A matrix of grid points, where we want to evaluate the density estimate
#' @export
dlg <- function(lg_object, grid = NULL) {

    # Do some checks first
    check_lg(lg_object)
    if(!is.null(grid)) {
        grid <- check_data(grid,
                           dim_check = ncol(lg_object$transformed_data),
                           type = "grid")
    }

    # Sample size and number of variables
    n <- nrow(lg_object$x)
    d <- ncol(lg_object$x)

    # Extract the data that we will us in the optimization, transform the grid if needed
    if(lg_object$transform_to_marginal_normality) {
        x <- lg_object$transformed_data
        transformed <- lg_object$trans_new(grid)
        x0 <- transformed$trans
        normalizing_constants <- transformed$normalizing_constants
    } else {
        x <- lg_object$x
        x0 <- grid
        normalizing_constants <- rep(1, nrow(grid))
    }

    # Extract the pairs from the list of bandwidths
    pairs <- lg_object$bw$joint[, c("x1", "x2")]

    # Initialize the matrices where we eventually are going to store the parameter estimates
    loc_mean <- matrix(0, ncol = d, nrow = nrow(grid))
    loc_sd <- loc_mean + 1
    loc_cor <- matrix(0, ncol = nrow(pairs), nrow = nrow(grid))

    # If method "5par_marginals_fixed" we estimate the marginals now.
    if(lg_object$est_method == "5par_marginals_fixed") {
        marginal_estimates <- dlg_marginal_wrapper(data_matrix = x,
                                                   eval_matrix = x0,
                                                   bw_vector = lg_object$bw$marginal)
        for(i in 1:d) {
            loc_mean[,i] <- marginal_estimates[[i]][, "mu"]
            loc_sd[,i] <- marginal_estimates[[i]][, "sig"]
        }
    }

    # If method is "5par", then we have only one pair, and we can do all of the estimation now
    if(lg_object$est_method == "5par") {
        estimate <- dlg_bivariate(x = x,
                                  eval_points = x0,
                                  bw = c(lg_object$bw$joint[1, "bw1"],
                                         lg_object$bw$joint[1, "bw2"]),
                                  est_method = lg_object$est_method)

        loc_mean <- cbind(estimate$par_est[, "mu_1"],
                          estimate$par_est[, "mu_2"])
        loc_sd <- cbind(estimate$par_est[, "sig_1"],
                        estimate$par_est[, "sig_2"])
        loc_cor[,1] <- estimate$par_est[, "rho"]        
    } else {  # And if not, we estimate the local correlation pairwise now
        for(i in 1:nrow(pairs)) {
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
            
            loc_cor[,i] <- pairwise_estimate$par_est[, "rho"]
        }
    }

    # Evaluate the density estimate
    f_est <- mvnorm_eval(eval_points = x0,
                         loc_mean = loc_mean,
                         loc_sd = loc_sd,
                         loc_cor = loc_cor,
                         pairs = pairs)*
        apply(normalizing_constants, 1, prod)

    # Return
    return(list(f_est = f_est,
                loc_mean = loc_mean,
                loc_sd = loc_sd,
                loc_cor = loc_cor,
                x = lg_object$x,
                transformed_data = lg_object$transformed_data,
                normalizing_constants = normalizing_constants,
                grid = grid,
                transformed_grid = x0))
    
}

#' The locally Gaussian conditional density estimator
#'
#' Estimate a conditional density function using locally Gaussian approximations.
#'
#' This function is the coditional version of the locally Gaussian density estimator
#' (LGDE), described in Otneim & Tjøstheim (2017) [Proper citation needed!]. The
#' function takes as arguments an \code{lg}-object as produced by the main \code{lg}-
#' function, a grid of points where the density estimate should be estimated, and a
#' set of conditions.
#'
#' The variables must be sorted befor they are supplied to this function. It will always
#' assume that the free variables come before the conditioning variables.
#' @param lg_object An object of type \code{lg}, as produced by the \code{lg}-function
#' @param grid A matrix of grid points, where we want to evaluate the density estimate
#' @param condition A vector with conditions for the variables that we condition upon
#' @export
clg <- function(lg_object, grid = NULL, condition = NULL) {

    # Extract some basic info
    n  <- nrow(lg_object$x)        # Sample size
    d  <- ncol(lg_object$x)        # Number of variables
    nc <- length(condition)        # Number of conditioning variables
    m  <- nrow(grid)               # Number of grid points
    
    # Do some checks first
    check_lg(lg_object)
    if(!is.null(grid)) {
        grid <- check_data(grid,
                           dim_check = d - nc,     # The grid must have d - c columns,
                           type = "grid")         # that is: the number of free variables
    }
    
    
    # Create the estimation grid, where we need the local correlation.
    estimation_grid <- matrix(NA, ncol = d, nrow = m)
    estimation_grid[, 1:(d-nc)] <- grid
    for(i in 1:nc) {
        estimation_grid[,d-nc+i] <- rep(condition[i], m)
    }

    # Estimate the local correlation in these points using the density function
    density_object <- dlg(lg_object, grid = estimation_grid)

    # In each grid point, we need to calculate the conditional mean vector and
    # covariance matrix in order to calculate the conditional density estimate.

    # Some useful quantities:
    pairs   <- lg_object$bw$joint[, c("x1", "x2")]                         # All pairs of variables
    n_pairs <- nrow(pairs)                                                 # The number of pairs
    z2      <- matrix(density_object$transformed_grid[1, (d - nc + 1):d],
                      ncol = 1) # The transformed location 
    
    # This is a function that does this task in grid point number i.
    f_eval <- function(i) {

        # This grid point
        grid_point <- density_object$transformed_grid[i, 1:(d - nc)]

        # First, this is the conditional mean at that point
        mu <- density_object$loc_mean[i,]
        # First, we build up the unconditional covariance matrix
        sigma <- diag(density_object$loc_sd[i,]^2)
        
        for(j in 1:n_pairs) { # For each pair of variables

            var1 <- pairs[j,1]      #| The variables
            var2 <- pairs[j,2]      #|

            rho <- density_object$loc_cor[i, j]     #| The local correlation in *this*
                                                    #| grid point, and for *this* variable
                                                    #| pair.

            # Finally, fill in the covariances
            sigma[var1, var2] <- rho*sqrt(sigma[var1, var1]*sigma[var2, var2])
            sigma[var2, var1] <- sigma[var1, var2]
        }

        # Next, we partition the mean vector and covariance matrix to produce the conditional
        # parameters
        mu1 <- matrix(mu[1:(d - nc)], ncol = 1)
        mu2 <- matrix(mu[(nc+1):d], ncol = 1)
        
        S11 <- as.matrix(sigma[1:(d - nc), 1:(d - nc)])
        S22 <- as.matrix(sigma[(d - nc + 1):d, (d - nc + 1):d])
        S12 <- as.matrix(sigma[1:(d - nc),  (d - nc + 1):d])
        dim(S12) <- c(d - nc, nc)
        S21 <- as.matrix(sigma[(d - nc + 1):d, 1:(d - nc)])
        dim(S21) <- rev(dim(S12))

        # The formulas for conditional mean and conditional covariance matrix:
        c_mean <- mu1 + S12%*%solve(S22)%*%(z2 - mu2)
        c_cov  <- S11 - S12%*%solve(S22)%*%S21

        f_est <- mvtnorm::dmvnorm(x = grid_point,
                                  mean = c_mean,
                                  sigma = c_cov)*
            prod(density_object$normalizing_constants[i, 1:(d - nc)])

        return(list(f_est = f_est,
                    c_mean = c_mean,
                    c_cov = c_cov))
    }

    # Apply the f_eval function to all the grid points
    estimate <- lapply(as.list(1:m), f_eval)

    # Extract the conditional density estimates into a vector
    f_est <- unlist(lapply(estimate, "[[", 1))

    # Extract the conditional means and covariance matrices into lists
    c_mean <- lapply(estimate, "[[", 2)
    c_cov  <- lapply(estimate, "[[", 3)

    # Return the result
    return(list(f_est = f_est,
                c_mean = c_mean,
                c_cov = c_cov,
                x = lg_object$x,
                transformed_data = lg_object$transformed_data,
                normalizing_constants = density_object$normalizing_constants,
                grid = grid,
                transformed_grid = density_object$transformed_grid))
}   
