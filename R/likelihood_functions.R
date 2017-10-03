

# The likelihood function is now defined outside this function
#' 1-parameter likelihood function for bivariate data
#'
#' \code{lik_1par} returns the likelihood function
#'
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
