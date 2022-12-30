#' Inverse Laplace transform in 2D
#'
#' @param t1 First argument to Laplace transform
#' @param t2 Second argument to Laplace transform
#' @param n Accuracy parameter
#' @param h.star Laplace transform to be inverted
#' @param ... Additional arguments to h.star
#'
#' @return Inverse Laplace transform
#' @export
#'
#' @examples inverseLaplace2d(1,1,function(s1,s2){1/((s1+1)*(s2+1))})
inverseLaplace2d <- function(t1, t2, h.star, ..., n = 15){
  #data('coef')
  coef$beta <- coef$beta[1:n]
  coef$eta <- coef$eta[1:n]
  m <- length(coef$eta)
  h4 <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      h4 = h4 + coef$eta[i]*coef$eta[j]*h.star(coef$beta[i]/t1, coef$beta[j]/t2, ...)
    }
  }
  Re(h4)
}




