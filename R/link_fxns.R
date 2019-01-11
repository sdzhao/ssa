#' Logit function
#'
#' @param x value between 0 and 1
#'
#' @return log(x/(1-x))
#'
#' @examples
#' x <- 0.3;
#' logit(x);
#'
#' @export

logit <- function(x){
    log(x/(1-x));
}

#' Expit function
#'
#' @param x value between -infinity and infinity
#'
#' @return 1/(1+exp(-x))
#'
#' @examples
#' x <- 1;
#' expit(x);
#'
#' @export

expit <- function(x){
    1/(1+exp(-x));
}
