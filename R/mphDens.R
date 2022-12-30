utils::globalVariables(c("constants"))
#' mphDens
#'
#' @importFrom parallel mcmapply
#' @importFrom utils globalVariables
#' @importFrom inverseLaplace mphDensity
#'
#' @param t1 First vector argument to inverse Laplace transform
#' @param t2 Second vector argument to inverse Laplace transform
#' @param alpha Initial probability vector
#' @param U Inverse generator
#' @param R Reward matrix
#' @param n.cores Number of cores to be used
#' @param n Approximation accuary
#'
#' @return A vector of density evaluations in the points (t1,t2)
#' @export
#'
#' @examples S <- matrix(c(-1.5, 0, 0,1.5, -1, 0,0, 1, -0.5), ncol = 3)
#' U <- solve(-S)
#' alpha <- c(0.2,0.4,0.4)
#' R <- matrix(c(0.2,0.5,0.8,0.4,0.1,0.9),nrow = 3)
#' mphDens(1, 1, alpha, U, R)
mphDens <- function(t1, t2, alpha, U, R, n.cores = 1, n = 15){
  pars <- constants[[n]]
  eta = c(pars$c * pars$mu1, (pars$a * pars$mu1 + 1i * pars$b * pars$mu1) /
            2)
  eta <- c(eta[1], rbind(c(eta[-1], Conj(eta[-1]))))
  beta = c(1, 1 + 1i * (1:pars$n) * pars$omega) * pars$mu1
  beta <- c(beta[1], rbind(c(beta[-1], Conj(beta[-1]))))

  mcmapply(mphDensity, t1, t2, MoreArgs = list(eta = eta, beta = beta, alpha = alpha, U = U, R = R), mc.cores = n.cores)
}
