utils::globalVariables(c("constants"))
#' Inverse Laplace transform in 1D
#'
#' @importFrom parallel mcmapply
#' @importFrom utils globalVariables
#'
#' @param t A vector of values for which to evaluate the inverse Laplace transform
#' @param h.star FUN Laplace transform to invert
#' @param ... Additional arguments to h.star
#' @param n Accuracy parameter
#' @param n.cores Number of cores to use defaults to 1
#'
#' @return Inverse Laplace transform as a vector
#' @export
#'
#' @examples inverseLaplace1d(1, function(s){s/(s^2+1)})
inverseLaplace1d <- function(t, h.star, ..., n = 15, n.cores = 1){
  pars <- constants[[n]]
  eta = c(pars$c * pars$mu1, (pars$a * pars$mu1 + 1i * pars$b * pars$mu1) /
            2)
  eta <- c(eta[1], rbind(c(eta[-1], Conj(eta[-1]))))
  beta = c(1, 1 + 1i * (1:pars$n) * pars$omega) * pars$mu1
  beta <- c(beta[1], rbind(c(beta[-1], Conj(beta[-1]))))
  m <- length(eta)
  mcmapply(function(t1, t2, ...) {
    h4 <- 0
    for (i in 1:m) {
      h4 = h4 + eta[i] * h.star(beta[i] / t, ...)
    }
    Re(h4/t)
  },
  t,
  MoreArgs = list(...),
  mc.cores = n.cores)
}
