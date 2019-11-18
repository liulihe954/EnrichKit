#' A function systematically performs hypergeomatric-based over-representation analysis in KEGG Database.
#'
#' @param sig_genes_all A list of single or multiple elements, contains all back ground genes for each test, should have the same length with sig_genes_all.
#' @param total_genes_all A list of single or multiple elements, contains all significant genes for each test, should have the same length with total_genes_all.
#' @param TestingSubsetNames A vector that has the names of each module, indicates the number of test will be run, should have the same length with total_genes_all/sig_genes_all.
#' @param KEGGthres P-value threshold to determain significance of each records.
#' @param species species = "bta"
#' @param id.type id.type = "kegg"
#' @param Sig_list_out A list contains all identifier conversion information, used to print out the significant genes in each enrichment records, in order to have the same type of identifier with imput.
#' @param keyword A keyword that will be included in all the outputs, helps user to recognize.
#' @import stats recount limma sva WGCNA stats ggplot2 dplyr tidyr readr tibble stringr edgeR biomaRt readxl GOSemSim corrplot org.Bt.eg.db MeSH.db MeSH.Bta.eg.db meshr magrittr AnnotationDbi openxlsx
#' @importFrom grDevices pdf png dev.off
#'
#' @return Nothing, but will output a .RData file will all the enrichment records.
#' @export
#'

Kegg_Enrich_Plot = function(sig_genes_all,
                            total_genes_all,# all genes in your dataset( vector - format - vector - ensembl iD)
                            TestingSubsetNames, # "sig" module names
                            KEGGthres = 0.05, # significant level (default - 0.05)
                            species = "bta",
                            id.type = "kegg",
                            Sig_list_out =Sig_list_out,
                            #biomart="ENSEMBL_MART_ENSEMBL",
                            #dataset="btaurus_gene_ensembl",
                            #host="http://www.ensembl.org",
                            #attributes = c("ensembl_gene_id","entrezgene_id"), # the items you need to retrive from the database
                            #filters="ensembl_gene_id", # with which keywords we match
                            keyword){
  #
  total_enrich = 0
  raw_pvalue_all = numeric()
  KEGG_results_b = list()
  KEGG_results_b_raw = list()
  sdb = kegg.gsets(species = species, id.type = id.type, check.new = F) # get database1
  kegg.gs = sdb$kg.sets[sdb$sigmet.id] # organize database
  # get all genes in database
  # for splitting
  rexp <- "^(\\w+)\\s?(.*)$"
  genesKEGG = character()
  KEGGID = character()
  KEGGTERM = character()
  for (i in seq_along(names(kegg.gs))){
    tmp_gene = unlist(kegg.gs[i]);attributes(tmp_gene) = NULL
    KEGGID[i] = sub(rexp,"\\1",names(kegg.gs)[i])
    KEGGTERM[i] = sub(rexp,"\\2",names(kegg.gs)[i])
    genesKEGG = append(genesKEGG,tmp_gene,length(genesKEGG))
  }
  genesKEGG = unique(genesKEGG)
  # the output is a pdf and every single page will be the point plot of the enriched item of a specific module.
  # pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    if (i%%1==0){message("Now digging in module #",i)} # can change the
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesKEGG])
    S = length(sig.genes[sig.genes %in% genesKEGG]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(ID=character(),
                     Term=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    # Double loop: trying to go through every single KEGG.db, so extract each one first
    for (j in 1:length(names(kegg.gs))){
      #KEGG_Index = unlist(str_split(names(kegg.gs)[j]," ",2))[1] # split to get the GO-index
      #KEGG_Name = unlist(str_split(names(kegg.gs)[j]," ",2))[2] # split to get the Go name
      all_ENTER_temp = (as.vector(unlist(kegg.gs[j]))) #
      if (j%%100==0){message("checking on KEGG #",j,"-",KEGGID[j],"-",KEGGTERM[j])}
      # Calculate and overlap
      m = length(total.genes[total.genes %in% all_ENTER_temp]) # genes from target GO and in our dataset
      findG = sig.genes[sig.genes %in% all_ENTER_temp]
      s = length(findG) # # genes from target GO also in the non-preserved module
      orig_list = data.frame(Sig_list_out[[i]]) %>% dplyr::filter(ENTREZID_final %in% findG)
      PastefindG = paste(orig_list[,1], collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value, digits = 100)
      tmp = data.frame(ID= KEGGID[j],
                       Term = KEGGTERM[j],
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
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("KEGGID","KEGGTERM", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    KEGG_results_b_raw[[i]] = final_raw; names(KEGG_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched kegg raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= KEGGthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("KEGGID","KEGGTERM", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    KEGG_results_b[[i]] = final;names(KEGG_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched KEGG")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    # plotting
    #print(final %>%
    #        top_n(dim(final)[1], wt= -pvalue) %>%
    #        mutate(hitsPerc=Significant_Genes*100/Total_Genes)%>%
    #        ggplot(aes(x=hitsPerc,
    #                   y=KEGG_Name,
    #                   colour=pvalue,
    #                   size=Significant_Genes)) +
    #        xlim(0,max(final$hitsPerc)+5)+
    #        geom_point() +
    #        theme_gray()+
    #        labs(title= paste("KEGG Enrichment in module",module_name), x="Hits (%)", y="Kegg term", colour="p value", size="Count")+
    #       theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #        theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5)))
  }
  #dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_GO = raw_pvalue_sum)
  save(KEGG_results_b,
       KEGG_results_b_raw,
       raw_pvalue_distribution,
       file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significantly pathways found within ",
          length(TestingSubsetNames)," modules/subsets",
          " at the significance level of ",KEGGthres)
  message("Nice! - KEGG enrichment finished and data saved")}
