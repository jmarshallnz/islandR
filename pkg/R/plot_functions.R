get_var = function(post, variable) {
  ncol = 1;
  if (!is.null(dim(post[[1]][[variable]]))) {
    ncol = ncol(post[[1]][[variable]])
    lapply(1:ncol, function(col) { do.call(rbind, lapply(post, function(x) { x[[variable]][,col] })) })
  } else {
    list(do.call(rbind, lapply(post, function(x) { as.numeric(x[[variable]]) })))
  }
}

#' Extract a variable from the posterior of an attribution fit
#' @export
#' @param x an attribution object, fit using attribution.
#' @param variable a character string specifying the variable to extract.
#' @return a list of posterior values for the given variable
get_variable = function(x, variable) {
  get_var(x$posterior, variable)
}
