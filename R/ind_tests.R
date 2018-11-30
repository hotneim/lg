
#' Independence tests
#'
#' Independence tests based on the local Gaussian correlation
#'
#' Implementation of three independence tests: For iid data (Berentsen et al., 2014),
#' for serial dependence within a time series (Lacal and Tjøstheim, 2017a), and
#' for serial cross-dependence between two time series (Lacal and Tjøstheim,
#' 2017b). The first test has a different theoretical foundation than the latter
#' two, but the implementations are similar and differ only in the bootstrap
#' procedure. For the time series applications, the user must lag the series to
#' his/her convenience before making the \code{lg}_object and calling this
#' function.
#'
#' @param lg_object An object of type \code{lg}, as produced by the
#'   \code{lg_main}-function. The data must be two dimensional.
#' @param h The \code{h}-function used in the calculation of the test statistic.
#'   The default value is \code{h(x) = x^2}.
#' @param S The integration area for the test statistic. Must be a logical
#'   function that accepts an n x 2 matrix and returns TRUE if a row is in S.
#' @param bootstrap_type The bootstrap method. Choose "plain" for the ordinary
#'   nonparametric bootstrap valid for independence test for iid data and for
#'   serial dependence within a time series. Choose "stationary" or "block" for
#'   a test for cross dependence between two time series.
#' @param block_length Block length if using block bootstrap for the cross
#'   dependence test. Calculated by \code{np::b.star()} if not supplied.
#' @param n_rep Number of bootstrap replications.
#'
#' @return A list containing the test result as well as various parameters. The
#'   elements are:
#'
#'   \itemize{
#'     \item \code{lg_object} The lg-object supplied by the user.
#'     \item \code{observed} The observed value of the test statistic.
#'     \item \code{replicated} The replicated values of the test statistic.
#'     \item \code{bootstrap_type} The bootstrap type.
#'     \item \code{block_length} The block length used for the block bootstrap.
#'     \item \code{p_value} The p-value of the test.
#'   }
#'
#' @examples
#'
#'     # Remember to increase the number of bootstrap samplesin preactical
#'     # implementations.
#'
#'     \dontrun{
#'
#'     # Test for independence between two vectors, iid data.
#'     x1 <- cbind(rnorm(100), rnorm(100))
#'     lg_object1 <- lg_main(x1)
#'     test_result1 = ind_test(lg_object1,
#'                             bootstrap_type = "plain",
#'                             n_rep = 20)
#'
#'     # Test for serial dependence in time series, lag 1
#'     data(EuStockMarkets)
#'     logreturns <- apply(EuStockMarkets, 2, function(x) diff(log(x)))
#'     x2 <- cbind(logreturns[1:100,1], logreturns[2:101, 1])
#'     lg_object2 <- lg_main(x2)
#'     test_result2 = ind_test(lg_object2,
#'                             bootstrap_type = "plain",
#'                             n_rep = 20)
#'
#'     # Test for cross-dependence, lag 1
#'     x3 <- cbind(logreturns[1:100,1], logreturns[2:101, 2])
#'     lg_object3 <- lg_main(x3)
#'     test_result3 = ind_test(lg_object3,
#'                             bootstrap_type = "block",
#'                             n_rep = 20)
#'     }
#'
#' @references
#'
#'   Berentsen, Geir Drage, and Dag Tjøstheim. "Recognizing and visualizing
#'   departures from independence in bivariate data using local Gaussian
#'   correlation." Statistics and Computing 24.5 (2014): 785-801.
#'
#'   Lacal, Virginia, and Dag Tjøstheim. "Local Gaussian autocorrelation and
#'   tests for serial independence." Journal of Time Series Analysis 38.1
#'   (2017a): 51-71.
#'
#'   Lacal, Virginia, and Dag Tjøstheim. "Estimating and testing nonlinear local
#'   dependence between two time series." Journal of Business & Economic
#'   Statistics just-accepted (2017b).
#'
#' @export


ind_test <- function(lg_object, h = function(x) x^2,
                     S = function(y) as.logical(rep(1, nrow(y))),
                     bootstrap_type = "plain",
                     block_length = NULL,
                     n_rep = 1000) {

  # Check that the data is two-dimensional
  if(ncol(lg_object$x) != 2) {
      stop("The data must be two-dimensional.")
  }

  # Check the weight function S
  weight <- try(S(lg_object$x), silent = TRUE)

  if(!is.logical(weight)) {
      stop("The function S must return a logical vector when supplied with the data.")
  } else if(!all(weight %in% c(TRUE, FALSE))) {
      stop("The function S must return a logical vector when supplied with the data.")
  } else if(class(weight) == "try-error") {
      stop("The call S(x) fails. The function S must return a logical vector when supplied with the data.")
  }

  # Check the function h
  test_h <- try(h(c(-0.5, 0, 0.5)), silent = TRUE)

  if(class(test_h) == "try-error") {
      stop(paste("The call h(c(-0.5, 0, 0.5)) fails with error '",
                 test_h[1],
                 "': h must be a scalar function [-1, 1] -> R.",
                 sep = ""))
  } else if(length(test_h) != 3) {
      stop("h must be a scalar function [-1, 1] -> R.")
  }

  # Check that the supplied bootstrap type is allowed
  if(!(bootstrap_type %in% c("plain", "stationary", "block"))) {
      stop("bootstrap_type must be one of 'plain', 'stationary' or 'block'.")
  }

  # Calculate test statistic for the observed data
  dlg_object <- dlg(lg_object, grid = lg_object$x[weight, ])
  observed <- mean(h(c(dlg_object$loc_cor, rep (0, sum(!weight)))))

  # Store the replicated values here
  replicated <- rep(NA, n_rep)

  # Bootstrap
  if(bootstrap_type == "plain") {
    for(i in 1:n_rep) {

      # Generate a bootstrap sample
      x_replicated <- cbind(sample(lg_object$x[,1],
                                   size = nrow(lg_object$x),
                                   replace = TRUE),
                            sample(lg_object$x[,2],
                                   size = nrow(lg_object$x),
                                   replace = TRUE))

      replicated[i] <- ind_teststat(x_replicated, lg_object, S, h)
    }
  } else {

    # Determine the block length
    if(is.null(block_length)) {
      block_length <- np::b.star(lg_object$x)
      block_length[block_length < 2] <- 2
    }

    if(bootstrap_type == "stationary") {
      for(i in 1:n_rep) {

        # Generate the bootstrap-sample
        x_replicated <- cbind(tseries::tsbootstrap(lg_object$x[,1],
                                                   nb = 1,
                                                   type = "stationary",
                                                   b = block_length[1]),
                              tseries::tsbootstrap(lg_object$x[,2],
                                                   nb = 1,
                                                   type = "stationary",
                                                   b = block_length[2]))

        replicated[i] <- ind_teststat(x_replicated, lg_object, S, h)
      }
    } else {
        for(i in 1:n_rep) {

            # Generate the bootstrap-sample
            x_replicated <- cbind(tseries::tsbootstrap(lg_object$x[,1],
                                                       nb = 1,
                                                       type = "block",
                                                       b = block_length[3]),
                                  tseries::tsbootstrap(lg_object$x[,2],
                                                       nb = 1,
                                                       type = "block",
                                                       b = block_length[4]))

            replicated[i] <- ind_teststat(x_replicated, lg_object, S, h)
        }
    }
  }

  ret <- list(lg_object = lg_object,
              observed = observed,
              replicated = replicated,
              bootstrap_type = bootstrap_type,
              block_length = block_length,
              p_value = mean(replicated > observed))

  return(ret)
}

#' Function that calculates the test statistic in the independence tests.
#'
#' This is an auxiliary function used by the independence tests.
#'
#' @param x_replicated A sample.
#' @param lg_object An lg-object.
#' @param S Integration area, see \code{?ind_test}.
#' @param h h-function for test statistic, see \code{?ind_test}.
#'
ind_teststat <- function(x_replicated, lg_object, S, h) {
  weight_replicated <- S(x_replicated)

  # Create lg-object for the replicated data
  lg_object_replicated <-
    lg_main(x = x_replicated,
            bw_method = lg_object$bw_method,
            est_method = lg_object$est_method,
            transform_to_marginal_normality =
              lg_object$transform_to_marginal_normality,
            bw = lg_object$bw)

  # Calculate test statistic for the replicated data
  dlg_object_replicated <-
    dlg(lg_object_replicated,
        grid = lg_object_replicated$x[weight_replicated, ])

  mean(h(c(dlg_object_replicated$loc_cor,
           rep (0, sum(!weight_replicated)))))

}
