#####################################################
# ~ cutter: common utilities ~
#
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

# function strNthSplit ----------------------------------------------------

#' strNthSplit
#'
#' @param stri
#' @param split
#' @param nth
#'
#' @returns
#' @export
#'
#' @examples
strNthSplit <- function(stri,
                        split,
                        nth) {


  # confirm we are given string(s)
  stri <- as.character(unlist(stri))

  as.character(sapply(strsplit(stri, split=split),
                      function(s) {
                        s[nth]
                      }))
}


# function sem ------------------------------------------------------------

#' sem
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
sem <- function(x) {
  sd(x)/sqrt(length(x))
}
