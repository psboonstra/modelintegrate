# function that starts with piecewise linear growth but becomes O(x^(1/3)) for x > 20
get_num_eff = function(x) {
  (x <= 3) * 0.95 * x + 
    (x > 3) * (x <= 10) * (3 * 0.95 + (x - 3) * 0.75) + 
    (x > 10) * (x <= 20) * (3 * 0.95 + (10 - 3) * 0.75 + (x - 10) * 0.5) +
    (x > 20) * (3 * 0.95 + (10 - 3) * 0.75 + (20 - 10) * 0.5 + (pmax(20, x) - 20)^(1/3) - 0.5) 
}

