

#' 1-parameter likelihood function for bivariate data
#'
#' Did not bother to document this function properly
#'
#' @param rho rho
#' @param m1 m1
#' @param m2 m2
#' @param m3 m3
#' @param m4 m4
#' @param x1_0 x1_0
#' @param x2_0 x2_0
#' @param h1 h1
#' @param h2 h2
#'
#' @return The likelihood for the one-parameter case.
#'
#' @export

lik_1par <- function(rho,m1,m2,m3,m4,x1_0,x2_0,h1,h2) {
  - log(2*pi*sqrt(1 - rho^2))*m1 - m2/(2*(1 - rho^2)) - m3/(2*(1 - rho^2)) +
    rho*m4/(1 - rho^2) - 1/2*exp(-1/2*(x2_0^2*h1^2 + x1_0^2 + x1_0^2*h2^2 -
                                         2*x1_0*rho*x2_0 + x2_0^2)/(-rho^2 + h2^2 + 1 + h1^2 + h1^2*h2^2))/
    (pi*(-rho^2 + h2^2 + 1 + h1^2 + h1^2*h2^2)^(1/2))
}



#' 1-parameter likelihood maximization function for bivariate data
#'
#' Did not bother to document this function properly
#'
#' @param rho grid_point
#' @param x1 x1
#' @param x2 x2
#' @param h1 h1
#' @param h2 h2
#' @param tol tol
#'
#' @return The maximized likelihood for the one-parameter case.
#'
#' @export

maximize_likelihood_1par = function(grid_point,
                                    x1,x2,
                                    h1,h2,
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
