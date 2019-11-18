#' A function to retrieve all the pair-wise semantic similarity of each enriched GO record and visualize the relationship in the format of correlation plots.
#'
#' @param Enrich_Out A vector of significantly enriched GO terms.
#' @param GOthres The threshold used to plot.
#' @param label_sizeCC The font size of category CC
#' @param label_sizeBP The font size of category BP
#' @param label_sizeMF The font size of category MF
#' @param Database The database used to retrieve the semantic similarity measurement
#' @param measure Measurement of semantic similarity, consistant with \code{GOSemSim} package.
#' @param combine NULL
#' @param Dataset_Name A keyword that will be included in the output names, helps users to recognize.
#' @import stats recount limma sva WGCNA stats ggplot2 dplyr tidyr readr tibble stringr edgeR biomaRt readxl GOSemSim corrplot org.Bt.eg.db MeSH.db MeSH.Bta.eg.db meshr magrittr AnnotationDbi openxlsx
#' @importFrom grDevices pdf png dev.off
#' @return Nothing, but will generate PDF plots and the corresponding correlation matrix.
#' @export
#'

ReduceDim_GO_Plot = function(Enrich_Out,
                             GOthres = 0.001,
                             label_sizeCC = 0.4,
                             label_sizeBP = 0.4,
                             label_sizeMF = 0.4,
                             Database = "org.Bt.eg.db",
                             measure="Jiang",combine=NULL,
                             Dataset_Name){
  # load libraries + download ref database
  do.call(library,list(Database))
  semData_BP <- godata(paste(Database), ont="BP", computeIC=T)
  semData_MF <- godata(paste(Database), ont="MF", computeIC=T)
  semData_CC <- godata(paste(Database), ont="CC", computeIC=T)
  # selection + formating: for each category we have one vector containing all the sig GO terms
  BP_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "biological_process") %>%
    dplyr::select(go_id) %>% unlist();attributes(BP_List) = NULL # name is an attribute and we dont them, so set null
  CC_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "cellular_component") %>%
    dplyr::select(go_id) %>% unlist();attributes(CC_List) = NULL
  MF_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "molecular_function") %>%
    dplyr::select(go_id) %>% unlist();attributes(MF_List) = NULL
  ### Now we are trying to get all similarity matrix ready. N x N, symetric, diag = 1
  # For BP

  goSimMatrix_BP = GOSemSim::mgoSim(BP_List,
                                    BP_List,
                                    semData=semData_BP,measure=measure,combine = combine)
  suspectID_BP = rownames(goSimMatrix_BP)[is.na(goSimMatrix_BP[,1])]
  if (length(suspectID_BP) != 0){BP_List_new = setdiff(BP_List,suspectID_BP)
  message(length(suspectID_BP)," invalid ID captured in BP: ",suspectID_BP,", thus been removed!")
  } else {BP_List_new = BP_List;message("Nice! All IDs are valid in BP!")}
  goSimMatrix_BP_new = GOSemSim::mgoSim(BP_List_new,
                                        BP_List_new,
                                        semData=semData_BP,measure=measure,combine = combine)
  colnames(goSimMatrix_BP_new) = paste(BP_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% BP_List_new)])
  rownames(goSimMatrix_BP_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% BP_List_new)],BP_List_new)
  # For CC
  goSimMatrix_CC = GOSemSim::mgoSim(CC_List,
                                    CC_List,
                                    semData=semData_CC,measure=measure,combine = combine)
  suspectID_CC = rownames(goSimMatrix_CC)[is.na(goSimMatrix_CC[,1])]
  if (length(suspectID_CC) != 0){CC_List_new = setdiff(CC_List,suspectID_CC)
  message(length(suspectID_CC)," invalid ID captured in CC: ",suspectID_CC,", thus been removed!")
  } else {CC_List_new = CC_List;message("Nice! All IDs are valid in CC!")}
  goSimMatrix_CC_new = GOSemSim::mgoSim(CC_List_new,
                                        CC_List_new,
                                        semData=semData_CC,measure=measure,combine =combine)
  colnames(goSimMatrix_CC_new) = paste(CC_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% CC_List_new)])
  rownames(goSimMatrix_CC_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% CC_List_new)],CC_List_new)
  # For MF
  goSimMatrix_MF = GOSemSim::mgoSim(MF_List,
                                    MF_List,
                                    semData=semData_MF,measure=measure,combine = combine)
  suspectID_MF = rownames(goSimMatrix_MF)[is.na(goSimMatrix_MF[,1])]
  if (length(suspectID_MF) != 0){MF_List_new = setdiff(MF_List,suspectID_MF)
  message(length(suspectID_MF)," invalid ID captured in MF: ",suspectID_MF,", thus been removed!")
  } else {MF_List_new = MF_List;message("Nice! All IDs are valid in MF!")}
  goSimMatrix_MF_new = GOSemSim::mgoSim(MF_List_new,
                                        MF_List_new,
                                        semData=semData_MF,measure=measure,combine = combine)
  colnames(goSimMatrix_MF_new) = paste(MF_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% MF_List_new)])
  rownames(goSimMatrix_MF_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% MF_List_new)],MF_List_new)
  # Now we take the results and plot
  pdf(paste("Semantic_Similarity_Measure_",Dataset_Name,"_",formatC(GOthres, format = "e", digits = 0),".pdf",sep = ""))
  corrplot(goSimMatrix_CC_new,title = "Semantic_Similarity_Measure_CC",
           tl.col = "black", tl.cex = label_sizeCC,
           method = "shade", order = "hclust",
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  corrplot(goSimMatrix_BP_new,title = "Semantic_Similarity_Measure_BP",
           tl.col = "black", tl.cex = label_sizeBP,
           method = "shade", order = "hclust",
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  corrplot(goSimMatrix_MF_new,title = "Semantic_Similarity_Measure_MF",
           tl.col = "black", tl.cex = label_sizeMF,
           method = "shade", order = "hclust",
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  dev.off()
  message(dim(goSimMatrix_CC_new)[1],",",
          dim(goSimMatrix_BP_new)[1],",",
          dim(goSimMatrix_MF_new)[1]," GOs ploted in CC, BP and MF, respectively",
          ", cutting at, ",GOthres)
  CorMatrix <- list("CorMat_BP" = data.frame(goSimMatrix_BP),
                    "CorMat_CC" = data.frame(goSimMatrix_CC),
                    "CorMat_MF" = data.frame(goSimMatrix_MF))
  write.xlsx(CorMatrix, row.names=TRUE,
             file = paste("Semantic_Similarity_Measure_",
                          Dataset_Name,"_",
                          formatC(GOthres, format = "e", digits = 0),".xlsx",sep = ""))
  save(goSimMatrix_CC,goSimMatrix_BP,goSimMatrix_MF,
       file = paste("Semantic_Similarity_Measure_",Dataset_Name,"_",formatC(GOthres, format = "e", digits = 0),".RData",sep = ""))
  message("Nice! Excels, Plots exported and RData saved!")
}
