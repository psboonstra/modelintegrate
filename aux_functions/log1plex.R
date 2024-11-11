#' Log1plex
#' 
#' Safe function for accurately calculating `log(1 + exp(x))` even for very
#' large x's
#'
#' @param x A numeric vector
#'
#' @return `log(1 + exp(x))` for each element of x
#'
#' @examples
#' log1plex(1e5) # correct answer
#' log(1+exp(1e5)) # overflows
log1plex = function(x) {
  - plogis(-x, log.p = TRUE)
}
