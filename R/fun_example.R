
#' Title
#'
#' @param x 
#' @param option 
#'
#' @return
#' @export
#'
#' @examples
some_function = function(x, option = c("r", "eigen")) {
  if (option == "r") {
    message("using base r")
    return(x^2)
  } else {
    message("using rcpp eigen")
    return(square_rcppeigen_internal(x))
  }
}