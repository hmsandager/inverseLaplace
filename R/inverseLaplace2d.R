utils::globalVariables(c("constants"))
#' Inverse Laplace transform in 2D
#'
#' @importFrom parallel mcmapply
#' @importFrom utils globalVariables
#'
#' @param t1 First argument to Laplace transform
#' @param t2 Second argument to Laplace transform
#' @param n Accuracy parameter
#' @param h.star FUN Laplace transform to be inverted
#' @param ... Additional arguments to h.star
#' @param n.cores Number of cores to use defaults to 1
#'
#' @return Inverse Laplace transform
#' @export
#'
#' @examples inverseLaplace2d(c(1,2),c(1,2),function(s1,s2){1/((s1+1)*(s2+1))})
inverseLaplace2d <-
  function(t1,
           t2,
           h.star,
           ...,
           n = 15,
           n.cores = 1) {
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
        for (j in 1:m) {
          h4 = h4 + eta[i] * eta[j] * h.star(beta[i] / t1, beta[j] / t2, ...)
        }
      }
      Re(h4/(t1*t2))
    },
    t1,
    t2,
    MoreArgs = list(...),
    mc.cores = n.cores)
  }




