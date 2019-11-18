#' Parse results from over-representation test, create a new data frame contains all significant hits
#'
#' @param Results_List A Returning list of single or multiple element.
#' @param keyword A note just to remind user of the database name, will be included in the prompts.
#'
#' @return A data frame contains all enriched records under pre-defined significance value
#' @import stats recount limma sva WGCNA stats ggplot2 dplyr tidyr readr tibble stringr edgeR biomaRt readxl GOSemSim corrplot org.Bt.eg.db MeSH.db MeSH.Bta.eg.db meshr magrittr AnnotationDbi openxlsx
#' @importFrom grDevices pdf png dev.off
#' @export
#'

Parse_Results = function(Results_List,keyword = "Which D.B"){
  all_enrich = data.frame()
  for (i in 1:length(Results_List)){
    len = dim(data.frame(Results_List[i]))[1]
    if (len> 0){
      tmp = data.frame(Results_List[i])
      names(tmp) = names(Results_List[[1]])
      all_enrich = rbind(all_enrich,tmp)
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich)[1]
  total_modules = length(Results_List)
  print(paste("In database: ",keyword,"-",total_hits,"hits found in",total_modules,"tested modules: ",names(Results_List)))
  return(ParseResults = all_enrich)
}
