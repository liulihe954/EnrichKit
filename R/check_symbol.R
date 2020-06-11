#' Suspicious gene symbol checking
#' @description Check whether the current gene symbols are the approved HGNC suggested gene symbols, associate the latter if not
#' @param notsure Single or multiple suspicious gene symbols to be checked. multiple elements should be in an vector.
#'
#' @return Single or multiple "suggested" gene symbols. If the original symbols are approved or invalid, return themselves as so, otherwise return the HGNC approved suggested symbols.
#'   Refer to \code{\link{checkGeneSymbols}} from package \code{\link{HGNChelper}}
#'
#' @export
#' @import HGNChelper dplyr
#'
#' @examples check_symbol(c("TNK2","ACR"))
#'
check_symbol = function(notsure){
  tmp = checkGeneSymbols(notsure) %>%
    mutate(Suggested.Symbol.Merge = ifelse(is.na(Suggested.Symbol),x,Suggested.Symbol))
  out = unlist(tmp$Suggested.Symbol.Merge,use.names = F)
  return(out)
}



