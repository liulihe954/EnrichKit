#' Parse/split a tibble into a list according to one column.
#' Ideally one column has categorical variable thus can be eatracted, each level would become a new element name of the list.
#'
#' @param tibble a tibble object that has at least two columns.
#' @param column String. The column used as index, preferably contains categorical variable.
#' @param keep String. The single or multiple column names to keep.
#'
#' @return A list. The selected column with preferably differetn levels would be used as index and the rest or selected columns would be kept.
#' The tibble would become a list.
#' @export
#' @import tibble dplyr
#'
#' @examples
#' df = tibble(A = c(1:5),B = c(rep('M',3),rep('T',2)))
#' split_tibble(df,'B','A')
#'
#'
split_tibble <- function(tibble,column,keep) {
  tibble %>%
    split(., .[,column]) %>%
    lapply(., function(x) c(unlist(x[,keep],use.names = F)))
}

