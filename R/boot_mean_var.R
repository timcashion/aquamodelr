#' Returns mean and variance of new sample
#'
#' \code{boot_mean_var} This is used to compute a mean and variance within a bootstrap function
#' @param x Numeric vector
#' @param d None
#' @return mean and variance of x
#' @examples
#' boot_mean_var(x= c(1,2,2,2,2,3))
#' @export

boot_mean_var <- function(x, d) {
  return(c(mean(x[d]),var(x[d])))
}
