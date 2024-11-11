#' Get prior number of non-zero parameters
#' 
#' Starts with piecewise linear growth but becomes O(x^(1/3)) for x > 20. So if
#' you have a small number of variables, then you think that most of the
#' corresponding parameters will be non-zero, but if you have a large number of
#' variables, you think that most of the corresponding parameters will be zero.
#'
#' @param x An integer
#'
#' @return A number less than `x` that is used as the expected number of
#'   non-zero parameters for prior elicitation purposes.
#'
#' @examples
#' get_num_eff(3) # 2.85/3 = 95% non-zero
#' get_num_eff(300) # 19.1/300 = 6% non-zero
get_num_eff = function(x) {
  (x <= 3) * 0.95 * x + 
    (x > 3) * (x <= 10) * (3 * 0.95 + (x - 3) * 0.75) + 
    (x > 10) * (x <= 20) * (3 * 0.95 + (10 - 3) * 0.75 + (x - 10) * 0.5) +
    (x > 20) * (3 * 0.95 + (10 - 3) * 0.75 + (20 - 10) * 0.5 + (pmax(20, x) - 20)^(1/3) - 0.5) 
}

