#' RNA-seq data preparation for correlation network analysis - With pc correction
#'
#' @param networkData A dataframe contains expression data, samples in rows and genes in columns.
#' @param cousin Percentage of reads that need to be removed. Only a few genes have many reads thus we remove them.
#' @param n1 Number of samples in referece group.
#' @param n2 Number of samples in test group.
#' @param thres_rmzero  Lowest count considered as not expressed, default 5 reads per sample.
#' @param count_rmzero threshold of individuals to remove a gene, usually 1/2 sample size.
#' @param perct Percentage of genes need to be removed based on overall expression variance. Only genes with high variance will be included in network analysis.
#' @param correct  Whether pc-based correction should be used or not.
#' @return A dataframe with corrected measurement of expression.
#'
#' @import stats recount limma sva WGCNA stats ggplot2 dplyr tidyr readr tibble stringr edgeR biomaRt readxl GOSemSim corrplot org.Bt.eg.db MeSH.db MeSH.Bta.eg.db meshr magrittr AnnotationDbi
#' @importFrom grDevices pdf png dev.off
#'
#' @export
#'
DataPre_C = function(networkData, cousin = 0.4, n1, n2, perct, correct = T,thres_rmzero = 5,count_rmzero){
  check_zero = function(networkData,thres_rmzero,count_rmzero){
    cow_count_index = rep("ok",length(rownames(networkData)))
    for (i in seq_along(rownames(networkData))){
      tmp_count = sum(networkData[i,] <= thres_rmzero)
      if (tmp_count >= count_rmzero){cow_count_index[i] = "out"}
    }
    return(cow_count_index)
  }
  remove_filter = function(networkData,thres){
    ID_meanexpr1 = data.frame(names = rownames(networkData), mean = apply(networkData, MARGIN = 1,mean))
    ID_meanexpr2 = cbind(ID_meanexpr1,percent = ID_meanexpr1$mean/sum(ID_meanexpr1$mean))
    ID_meanexpr3 = ID_meanexpr2[order(ID_meanexpr2$mean,decreasing = T),]
    accumulative = numeric(nrow(ID_meanexpr3))
    for (i in c(1:nrow(ID_meanexpr3))){
      accum = sum(ID_meanexpr3$percent[1:i])
      accumulative[i] = accum
    }
    remove_pos = (length(which(accumulative <= thres))+1)
    remove_index = ID_meanexpr3$names[1:remove_pos]
    networkData_filter = networkData[!(rownames(networkData)%in%remove_index),]
    Results = list(remove_index=remove_index,networkData_filter = networkData_filter)
    return(Results)
  }

  q_normalize <- function(dat){
    n = nrow(dat)
    p = ncol(dat)
    rank.dat =  dat # matrix for ranking
    for (i in 1:p){
      rank.dat[,i] = rank(dat[,i])
    }
    U = rank.dat/(n+1)
    qnorm(U)
  }

  Correct_pca = function(rse_raw,method){
    rse_raw <- t(rse_raw)# transpose data so that rows are samples and columns are gene expression measurements
    mod=matrix(1,nrow=dim(rse_raw)[1],ncol=1)
    colnames(mod)="Intercept"
    ## num.sv requires data matrix with features(genes) in the rows and samples in the column
    nsv=num.sv(t(rse_raw), mod, method = method)
    print(paste("Number of PCs estimated to be removed:", nsv))
    ## PC residualization of gene expression data using sva_network. Rows correspond to samples, Columns correspond to features
    exprs_corrected = sva_network(rse_raw, nsv)
    ## Quantile normalize the corrected gene expression measurements such that each expression of every gene follows a gaussian distribution
    exprs_corrected_norm <- q_normalize(exprs_corrected)
    return(list(exprs_corrected_norm = t(data.frame(exprs_corrected_norm))))
  }
  # step 0 - rm too many zeros
  zero_cm_label = check_zero(networkData,thres_rmzero = 5,count_rmzero = 9)
  networkData_nozero = data_expr_all_with0[zero_cm_label=="ok",]
  # step 1 - filter out top 40% counts
  ## filter out top 40% counts # function established for future use
  networkData_filter = remove_filter(networkData,cousin)$networkData_filter
  # step 2 - normalization (0s out and normalization)
  remove_index = which(rowSums(networkData_filter) == 0)#;length(remove_index14)
  networkData_nm1 = networkData_filter[-remove_index,]#;dim(networkData_nm1)
  networkData_nmList = DGEList(counts = networkData_nm1,group  = c(rep("ref",n1),rep("test",n2)))
  networkData_nm2 = calcNormFactors(networkData_nmList)
  networkData_normalized_normfactors = networkData_nm2$samples
  networkData_normalized = data.frame(networkData_nm2$counts)
  #dim(networkData_normalized)
  # step 3 - log 2 trans
  # log trans
  networkData_log2 = log2(networkData_normalized+2)
  # step 4 - filter out bottom xx% variation
  # select most var
  networkData_log2$variance = apply(networkData_log2,1,var)
  networkData_log2_50var = networkData_log2[networkData_log2$variance >= quantile(networkData_log2$variance,c(perct)),]
  networkData_log2_50var$variance <- NULL
  #dim(networkData14_log2_50var)
  # step 5 - pca correction
  networkData_correction = Correct_pca(networkData_log2_50var,"leek")
  networkData_final = data.frame(networkData_correction$exprs_corrected_norm);
  names(networkData_final) = names(networkData_log2_50var)
  save(networkData_final,
       networkData_log2_50var,
       networkData_normalized_normfactors,
       networkData_normalized,
       file = paste(deparse(substitute(networkData)),"prepare with corrections","_top",100*(1-perct),".RData",sep = ""))
  if (correct == T) {return(list(Corrected_log2_PC = networkData_final))}
  else {return(list(Networkdata_norm = networkData_normalized))}
}
