#' A function systematically performs hypergeomatric-based over-representation analysis in Interpro Database.
#'
#' @param total_genes_all A list of single or multiple elements, contains all back ground genes for each test, should have the same length with sig_genes_all.
#' @param sig_genes_all A list of single or multiple elements, contains all significant genes for each test, should have the same length with total_genes_all.
#' @param TestingSubsetNames A vector that has the names of each module, indicates the number of test will be run, should have the same length with total_genes_all/sig_genes_all.
#' @param IPthres P-value threshold to determain significance of each records.
#' @param biomart Name of annotation database. for now, only have "ensembl".
#' @param dataset Name of annotation database. e.g. "btaurus_gene_ensembl".
#' @param Identifier Whether "enternal_gene_name" or "ensembl_gene_id" will be used as the identifier in the imput dataset.
#' @param attributes A vector contains all attributes to be retrieved from the database, e.g. \code{c("ensembl_gene_id","go_id","name_1006")}
#' @param keyword A keyword that will be included in all the outputs, helps user to recognize.
#' @import stats recount limma sva WGCNA stats ggplot2 dplyr tidyr readr tibble stringr edgeR biomaRt readxl GOSemSim corrplot org.Bt.eg.db MeSH.db MeSH.Bta.eg.db meshr magrittr AnnotationDbi openxlsx
#' @importFrom grDevices pdf png dev.off
#' @return Nothing, but will output a .RData file will all the enrichment records.
#' @export
#'
InterPro_Enrich = function(total_genes_all,
                           sig_genes_all,
                           TestingSubsetNames,
                           IPthres = 0.05,
                           biomart="ensembl",
                           dataset="btaurus_gene_ensembl",
                           Identifier = "external_gene_name",
                           attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                           keyword = "Interpro_Enrichment"){
  total_enrich = 0
  raw_pvalue_all = numeric()
  Interpro_results_b = list()
  Interpro_results_b_raw = list()
  DB_List = list()
  database = useMart(biomart)
  genome = useDataset(dataset, mart = database)
  gene = getBM(attributes,mart = genome)
  ##
  InterproName = unique(gene[,c("interpro","interpro_description")]) %>% arrange(interpro)
  Interpro = na.omit(InterproName$interpro)[-1]
  Name = na.omit(InterproName$interpro_description)[-1]
  #
  if (Identifier == "ensembl_gene_id"){
    genesInterpro = unique(subset(gene,interpro != "")$ensembl_gene_id)
    for ( p in seq_along(Interpro)){
      tmp = subset(gene, interpro == Interpro[p])$ensembl_gene_id
      DB_List[[p]] = tmp #
      names(DB_List)[p]  <- paste(Interpro[p],"-",Name[p])
    }
  } else if (Identifier == "external_gene_name") {
    genesInterpro= unique(subset(gene,interpro != "")$external_gene_name);
    genesInterpro = genesInterpro[-1]
    for ( p in seq_along(Interpro)){
      tmp = subset(gene, interpro == Interpro[p])$external_gene_name
      DB_List[[p]] = tmp #
      names(DB_List)[p]  <- paste(Interpro[p],"-",Name[p])
    }
  } else {message("Sorry, we only have ensembel and names available as identifier, please use one of the followings:
                  ensembl_gene_id OR external_gene_name.")}

  #length(genesGO)
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Interpro domains to check: ",length(Interpro)," with total number of names: ",length(Name))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesInterpro])
    S = length(sig.genes[sig.genes %in% genesInterpro]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(Interpro=character(),
                     Name=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(Interpro)){
      if (j%%100 == 0) {message("tryingd on Interpro ",j," - ",Interpro[j]," - ",Name[j])}
      gENEs = DB_List[[j]]
      m = length(total.genes[total.genes %in% gENEs]) # genes from target interpro and in our dataset
      findG = sig.genes[sig.genes %in% gENEs]
      s = length(findG) # # genes from target interpro also in the non-preserved module
      PastefindG = paste(findG, collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(Interpro = Interpro[j],
                       Name = Name[j],
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
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("InterproID","Interpro_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Interpro_results_b_raw[[i]] = final_raw; names(Interpro_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Interpro raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= IPthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("InterproID","Interpro_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    Interpro_results_b[[i]] = final;names(Interpro_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched Interpro")
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
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_Interpro = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(Interpro_results_b, Interpro_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant Interpro domains found within ",
          length(TestingSubsetNames)," modules/subsets",
          " at the significance level of ",IPthres)
  message("Nice! - Interpro enrichment finished and data saved")}
