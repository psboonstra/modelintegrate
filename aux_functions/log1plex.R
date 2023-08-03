# Safe function for accurately calculating log(1 + exp(x)) even for very large
# and small x's
log1plex = function(x) {
  - plogis(-x, log.p = TRUE)
  #ifelse(x > 0,
  #       x + log1p(exp(-x)),
  #       log1p(exp(x)))
}
