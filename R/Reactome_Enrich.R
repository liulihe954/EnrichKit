#' A function systematically performs hypergeomatric-based over-representation analysis in Reactome Database.
#'
#' @param total_genes_all A list of single or multiple elements, contains all back ground genes for each test, should have the same length with sig_genes_all.
#' @param sig_genes_all A list of single or multiple elements, contains all significant genes for each test, should have the same length with total_genes_all.
#' @param TestingSubsetNames A vector that has the names of each module, indicates the number of test will be run, should have the same length with total_genes_all/sig_genes_all.
#' @param InputSource One of the three database input category, all pathways, all reactions, lowest pathways
#' @param Sig_list_out A list contains all identifier conversion information, used to print out the significant genes in each enrichment records, in order to have the same type of identifier with imput.
#' @param Reacthres P-value threshold to determain significance of each records.
#' @param keyword A keyword that will be included in all the outputs, helps user to recognize.
#' @import stats recount limma sva WGCNA stats ggplot2 dplyr tidyr readr tibble stringr edgeR biomaRt readxl GOSemSim corrplot org.Bt.eg.db MeSH.db MeSH.Bta.eg.db meshr magrittr AnnotationDbi openxlsx
#' @importFrom grDevices pdf png dev.off
#'
#' @return Nothing, but will output a .RData file will all the enrichment records.
#'
#' @export
#'

Reactome_Enrich = function(total_genes_all,
                           sig_genes_all,
                           TestingSubsetNames,
                           InputSource,
                           Sig_list_out,
                           Reacthres = 0.05,
                           #biomart="ensembl",
                           #dataset="btaurus_gene_ensembl",
                           #host="http://www.ensembl.org",
                           #attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                           keyword = "Reactome_Enrichment"){
  total_enrich = 0
  raw_pvalue_all = numeric()
  Reactome_results_b = list()
  Reactome_results_b_raw = list()
  DB_List = list()
  Reactome_gene =   unique(InputSource[,c("EntrezID")])
  ReactomeRecords = dplyr::select(InputSource,ReactomeID,Reactome_Description) %>% dplyr::arrange(ReactomeID) %>% distinct()
  #ReactomeRecords = unique(InputSource[,c("ReactomeID","Reactome_Description")]) %>% arrange(ReactomeID) #
  ReactomeID = na.omit(ReactomeRecords$ReactomeID)
  ReactomeName = na.omit(ReactomeRecords$Reactome_Description)
  for ( p in seq_along(ReactomeID)){
    IDindex = ReactomeID[p]
    tmp = subset(InputSource, ReactomeID == IDindex)$EntrezID
    DB_List[[p]] = tmp #
    names(DB_List)[p]  <- paste(ReactomeID[p],"-",ReactomeName[p])
  }
  #ReactomeID = ReactomeID[1:300]
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Reactome to check: ",length(ReactomeID)," with total number of names: ",length(ReactomeName))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    #i = 1
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    #length(sig.genes)
    #length(total.genes)
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% Reactome_gene])
    S = length(sig.genes[sig.genes %in% Reactome_gene]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(ReactomeID=character(),
                     ReactomeTerm=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(ReactomeID)){
      # j = 101
      if (j%%100 == 0) {message("tryingd on Reactome ",j," - ",ReactomeID[j]," - ",ReactomeName[j])}
      #target = ReactomeID[j]
      #gENEs = unique(subset(InputSource, ReactomeID == target)$EntrezID)
      gENEs = DB_List[[j]]
      m = length(total.genes[total.genes %in% gENEs])
      findG = sig.genes[sig.genes %in% gENEs]
      s = length(findG)
      orig_list = data.frame(Sig_list_out[[i]]) %>% dplyr::filter(ENTREZID_final %in% findG)
      PastefindG = paste(orig_list[,1], collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(ReactomeID = ReactomeID[j],
                       ReactomeName = ReactomeName[j],
                       totalG = m,
                       sigG = s,
                       Pvalue = Pval,
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig,
                       findG = PastefindG)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("ReactomeID","ReactomeName", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Reactome_results_b_raw[[i]] = final_raw; names(Reactome_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Reactomeid raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= Reacthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("ReactomeID","ReactomeName", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    Reactome_results_b[[i]] = final;names(Reactome_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched ReactomeID")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    ##
    #   print(final %>%
    #           top_n(dim(final)[1], wt= -pvalue)%>%
    #           ggplot(final, aes( x = hitsPerc,
    #                      y = GO_Name,
    #                      colour = pvalue,
    #                      size = Significant_Genes)) +
    #           geom_point() +
    #           theme_gray()+
    #          labs(title= paste("GO Enrichment in module",
    #                              TestingSubsetNames[i])),
    #                              x="Hits (%)", y="GO term",
    #                              colour="p value", size="Count")+
    #      theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
  }
  #  dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_Reactome = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(Reactome_results_b, Reactome_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant Reactome domains found within ",
          length(TestingSubsetNames)," modules/subsets",
          " at the significance level of ",Reacthres)
  message("Nice! - Reactome enrichment finished and data saved")}
