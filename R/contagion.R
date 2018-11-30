#' Test for financial contagion
#'
#' Test for financial contagion by means of the local Gaussian correlation.
#'
#' This function is an implementation of the test for financial contagion
#' developed by Støve, Tjøstheim and Hufthammer (2013). They test whether the
#' local correlations between two financial time series are different before and
#' during crisis times. The distinction between crisis and non-crisis times must
#' be made by the user.
#'
#' @param lg_object_nc An object of type \code{lg}, as produced by the
#'   \code{lg_main}-function for the observations covering the non-crisis
#'   period. The data must be two dimensional.
#' @param lg_object_c An object of type \code{lg}, as produced by the
#'   \code{lg_main}-function for the observations covering the crisis period.
#'   The data must be two dimensional.
#' @param grid_range This test measures the local correlations a long the
#'   diagonal specified by this vector of length two.
#' @param grid_length The number of grid points.
#' @param n_rep The number of bootstrap replicates.
#' @param weight Weight function
#'
#' @return A list containing the test result as well as various parameters. The
#'   elements are:
#'
#'   \itemize{
#'     \item \code{observed} The observed value of the test statistic.
#'     \item \code{replicated} The replicated values of the test statistic.
#'     \item \code{p_value} The p-value of the test.
#'     \item \code{local_correlations} The local correlations measured along the
#'      diagonal, for the non-crisis and crisis periods respectively.
#'   }
#'
#' @examples
#'
#'    # Run the test on some built-in stock data
#'    data(EuStockMarkets)
#'    x <- apply(EuStockMarkets, 2, function(x) diff(log(x)))[, 1:2]
#'
#'    # Define the crisis and non-crisis periods (arbitrarily for this simple
#'    # example)
#'    non_crisis <- x[1:100, ]
#'    crisis     <- x[101:200, ]
#'
#'    # Create the lg-objects, with parameters that match the applications in the
#'    # original publication describibg the test
#'    lg_object_nc <- lg_main(non_crisis, est_method = "5par",
#'                            transform_to_marginal_normality = FALSE)
#'    lg_object_c  <- lg_main(crisis, est_method = "5par",
#'                            transform_to_marginal_normality = FALSE)
#'
#'    \dontrun{
#'    # Run the test (with very few resamples for illustration)
#'    test_result <- cont_test(lg_object_nc, lg_object_c,
#'                             n_rep = 10)
#'    }
#'
#' @references
#'
#'    Støve, Bård, Dag Tjøstheim, and Karl Ove Hufthammer. "Using local Gaussian
#'    correlation in a nonlinear re-examination of financial contagion." Journal
#'    of Empirical Finance 25 (2014): 62-82.
#'
#' @export

cont_test <- function(lg_object_nc, lg_object_c,
                      grid_range = quantile(rbind(lg_object_nc$x,
                                                  lg_object_c$x),
                                            c(.05, .95)),
                      grid_length = 30, n_rep = 1000,
                      weight = function(y) {rep(1, nrow(y))}) {

    # Do sanity checks of the arguments
    check_lg(lg_object_nc)
    check_lg(lg_object_c)

    # Both data sets must be two-dimensional
    if((ncol(lg_object_nc$x) != 2) | (ncol(lg_object_c$x) != 2)) {
      stop("Both data sets must have two columns.")
    }

    # Grid range must be a numeric vector of length two.
    if(!is.numeric(grid_range)) {
      stop("grid_range must be a numeric vector of length two.")
    } else if(length(grid_range) != 2) {
      stop("grid_range must be a numeric vector of length two.")
    }

    # Grid length must be a single number
    if(!is.numeric(grid_length)) {
      stop("grid_length must be a single number.")
    } else if(length(grid_length) != 1) {
      stop("grid_length must be a single number.")
    }

    # The number of replicates must also be a single number
    if(!is.numeric(n_rep)) {
      stop("grid_length must be a single number.")
    } else if(length(n_rep) != 1) {
      stop("grid_length must be a single number.")
    }

    # Some initialization
    replicated <- rep(NA, length = n_rep)
    grid <- cbind(seq(from = grid_range[1],
                      to = grid_range[2],
                      length.out = grid_length),
                  seq(from = grid_range[1],
                      to = grid_range[2],
                      length.out = grid_length))
    weights <- weight(grid)

    # Calculate the local correlation for the observed data
    dlg_object_nc <- dlg(lg_object_nc,
                                 grid = grid)
    dlg_object_c <- dlg(lg_object_c,
                             grid = grid)

    local_correlations <-
        data.frame(x = grid[,1],
                   non_crisis = dlg_object_nc$loc_cor,
                   crisis = dlg_object_c$loc_cor)

    observed <- mean((local_correlations$crisis -
                        local_correlations$non_crisis)*weights)

    # Do the bootstrapping
    complete_sample <- rbind(lg_object_nc$x,
                             lg_object_c$x)

    for(i in 1:n_rep) {

        resampled <-
            complete_sample[sample(1:nrow(complete_sample), replace = TRUE),]

        # Split the resampled data into non-crisis and crisis periods
        resampled_non_crisis <- resampled[1:nrow(lg_object_nc$x),]
        resampled_crisis     <- resampled[nrow((lg_object_nc$x) + 1):
                                              nrow(complete_sample),]

        # Create new lg-objects with the same parameters as for the observed
        # data
        lg_object_nc_resampled <-
            lg_main(x = resampled_non_crisis,
                    bw = lg_object_nc$bw,
                    est_method = lg_object_nc$est_method,
                    transform_to_marginal_normality =
                        lg_object_nc$transform_to_marginal_normality)

        lg_object_c_resampled <-
            lg_main(x = resampled_crisis,
                    bw = lg_object_c$bw,
                    est_method = lg_object_c$est_method,
                    transform_to_marginal_normality =
                        lg_object_c$transform_to_marginal_normality)

        # Calculate the local correlation for the resampled data
        dlg_object_nc_resampled <- dlg(lg_object_nc_resampled,
                                               grid = grid)
        dlg_object_c_resampled <- dlg(lg_object_c_resampled,
                                           grid = grid)

        # The resampled value of the test statistic
        replicated[i] <- mean((dlg_object_c_resampled$loc_cor -
            dlg_object_nc_resampled$loc_cor)*weights)
    }

    ret <- list(observed = observed,
                replicated = replicated,
                p_value = mean(replicated >= observed),
                local_correlations = local_correlations)

    return(ret)
}
