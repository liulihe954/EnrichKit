#' A funcition to convert Ensembl ID to Entrez ID.
#'
#' @param bg_gene A vector of genes with Ensembl id as format
#' @param TestingSubsetNames A vector that has the names of each module, indicates the number of test will be run.
#' @param TestingModAssign A vector of module assignment of each gene id
#' @param keyword A keyword that will be included in all the outputs, helps user to recognize.
#' @import stats recount limma sva WGCNA stats ggplot2 dplyr tidyr readr tibble stringr edgeR biomaRt readxl GOSemSim corrplot org.Bt.eg.db MeSH.db MeSH.Bta.eg.db meshr magrittr AnnotationDbi
#' @importFrom grDevices pdf png dev.off
#' @importFrom grDevices dev.off pdf png
#' @return Nothing, but will output a .RData file will all the enrichment records.
#' @export
#'

ConvertNformat = function(bg_gene,
                          TestingSubsetNames,
                          TestingModAssign,
                          keyword = "Ensembl2Entrez_Convert"){
  # Get match information
  key.symbol = AnnotationDbi::keys(org.Bt.eg.db,  keytype = c("ENSEMBL"))
  entrezUniverse = AnnotationDbi::select(org.Bt.eg.db, as.character(key.symbol),
                                         columns = c("ENTREZID"),keytype = "ENSEMBL") %>%
    dplyr::distinct(ENSEMBL,.keep_all= TRUE)
  Gather_all = data.frame(ENSEMBL  =  bg_gene,
                          assign = TestingModAssign) %>%
    dplyr::left_join(entrezUniverse, by  = c("ENSEMBL" = "ENSEMBL"))
  names( Gather_all)[3] = "ENTREZID_final"
  #
  Sig_list_out = list();Total_list_out = list()
  Sig_list_out_entrez = list();Total_list_out_entrez = list()
  Sig_list_out_ens = list();Total_list_out_ens = list()
  for (i in seq_along(TestingSubsetNames)){
    #
    target = TestingSubsetNames[i]
    tmp01 = dplyr::filter(Gather_all,assign == target)
    names(tmp01)[3] = "ENTREZID_final"
    Sig_list_out[[i]] = tmp01;names(Sig_list_out)[i] = TestingSubsetNames[i]
    tmp1 = dplyr::select(tmp01,ENTREZID_final) %>% dplyr::distinct() %>% na.omit();attributes(tmp1) = NULL
    Sig_list_out_entrez[[i]] = tmp1
    names(Sig_list_out_entrez)[i] = TestingSubsetNames[i]
    tmp2 = dplyr::select(tmp01,ENSEMBL) %>% dplyr::distinct() %>% na.omit();attributes(tmp1) = NULL
    Sig_list_out_ens[[i]] = tmp2;names(Sig_list_out_ens)[i] = TestingSubsetNames[i]
    names(Sig_list_out_entrez)[i] = TestingSubsetNames[i]

  }
  Total_list_out_tmp = list(unique(Gather_all$ENTREZID_final));names(Total_list_out_tmp) = "total_genes_entrez"
  Total_list_out_entrez = rep(Total_list_out_tmp, length(Sig_list_out_entrez))
  Total_list_out_tmp2 = list(unique(Gather_all$ENSEMBL));names(Total_list_out_tmp2) = "total_genes_ens"
  Total_list_out_ens = rep(Total_list_out_tmp2, length(Sig_list_out_ens))

  save(Sig_list_out,
       Sig_list_out_entrez,Total_list_out_entrez,
       Sig_list_out_ens,Total_list_out_ens,
       file = paste(trimws(keyword),".RData",sep = ""))
  message("Nice! Conversion finished")
}
