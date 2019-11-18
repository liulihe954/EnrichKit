#' A function systematically performs hypergeomatric-based over-representation analysis in MeSH Database.
#'
#' @param total_genes_all A list of single or multiple elements, contains all back ground genes for each test, should have the same length with sig_genes_all.
#' @param sig_genes_all A list of single or multiple elements, contains all significant genes for each test, should have the same length with total_genes_all.
#' @param TestingSubsetNames A vector that has the names of each module, indicates the number of test will be run, should have the same length with total_genes_all/sig_genes_all.
#' @param Meshthres P-value threshold to determain significance of each records.
#' @param Sig_list_out A list contains all identifier conversion information, used to print out the significant genes in each enrichment records, in order to have the same type of identifier with imput.
#' @param MeshCate Mesh category. for now, only "D" and "G"
#' @param dataset Species specific. e.g. Bos Taurus - "MeSH.Bta.eg.db"
#' @param keyword A keyword that will be included in all the outputs, helps user to recognize.
#' @import stats recount limma sva WGCNA stats ggplot2 dplyr tidyr readr tibble stringr edgeR biomaRt readxl GOSemSim corrplot org.Bt.eg.db MeSH.db MeSH.Bta.eg.db meshr magrittr AnnotationDbi openxlsx
#' @importFrom grDevices pdf png dev.off
#' @return Nothing, but will output a .RData file will all the enrichment records.
#' @export

MeSH_Enrich = function(total_genes_all,
                       sig_genes_all,
                       TestingSubsetNames,
                       Meshthres = 0.05,
                       Sig_list_out,
                       MeshCate = c("D","G"),
                       #biomart="ensembl",
                       dataset="MeSH.Bta.eg.db",
                       #dataset= "btaurus_gene_ensembl",
                       #Identifier = "external_gene_name",
                       #attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"),
                       keyword = "MESH_Enrichment"){
  #total.genes = Total_list_out_entrez_test
  #sig.genes = Sig_list_out_entrez_test
  #TestingSubsetNames = TestingSubsetNames_test
  total_enrich = 0
  raw_pvalue_all = numeric()
  Mesh_results_b = list()
  Mesh_results_b_raw = list()
  DB_List = list()
  ### Three ways to get meshdb
  # 1 download from github: we are gonna use
  #githubURL <- "https://github.com/liulihe954/Repro_Estrous_0918/raw/master/MeshDB.RData"
  #githubURL <- "https://github.com/liulihe954/HeatStress0708/raw/master/MeshDB_new.RData"
  #load(url(githubURL))
  #if (all(MeshCate%in%c("G","D"))){list_Bta = dplyr::filter(list_Bta, MESHCATEGORY %in% MeshCate)}
  #else {
  #  message("Sorry, we only have G and D")
  #  message("Now reload the category you need, it will take a while...")
  #  KEY = keys(MeSH.db, keytype = "MESHID")
  #  List = select(MeSH.db, keys = KEY, columns = columns(MeSH.db), keytype = "MESHID")
  #  Match_List = dplyr::select(List, MESHID, MESHTERM)
  #  key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
  #  list_Bta = MeSHDbi::select(MeSH.Bta.eg.db, keys = key_Bta, columns = columns(MeSH.Bta.eg.db)[-4], keytype = "MESHID") %>%
  #    dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>% dplyr::filter(MESHCATEGORY %in% MeshCate) %>%
  #    dplyr::left_join(Match_List,by= c("MESHID" = "MESHID"))}
  # 2. match from the very begining (will take an hour or so)
  #KEY = keys(MeSH.db, keytype = "MESHID")
  #List = select(MeSH.db, keys = KEY, columns = columns(MeSH.db), keytype = "MESHID")
  #List = select(MeSH.db, keys = KEY[1:3], columns = columns(MeSH.db), keytype = "MESHID")
  #Match_List = dplyr::select(List, MESHID, MESHTERM)
  ##head(Match_List)
  ### Prepare Bta database
  #key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
  #list_Bta = MeSHDbi::select(MeSH.Bta.eg.db, keys = key_Bta, columns = columns(MeSH.Bta.eg.db)[-4], keytype = "MESHID") %>%
  #  dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>% dplyr::filter(MESHCATEGORY %in% MeshCate) %>%
  #  dplyr::left_join(Match_List,by= c("MESHID" = "MESHID"))

  # 3. alternatively, if you have them in your environment
  keyword_outer = "MeshDB"
  DB = paste(keyword_outer,".RData",sep = "")
  load(DB)
  #Sig_list_out_entrez_test2
  #Total_list_out_entrez_test2
  # Get index
  list_Bta = list_Bta[which(list_Bta$MESHCATEGORY %in% MeshCate),]
  #list_Bta = dplyr::filter(list_Bta,MESHCATEGORY%in%MeshCate)
  genesMesh = unique(list_Bta$GENEID)
  MeshRecords = unique(list_Bta[,c("MESHID","MESHTERM")]) %>% arrange(MESHID)
  MeshID = na.omit(MeshRecords$MESHID)
  MeshTerm = na.omit(MeshRecords$MESHTERM)
  for ( p in seq_along(MeshID)){
    tmp = subset(list_Bta, MESHID == MESHID[p])$GENEID
    DB_List[[p]] = tmp #
    names(DB_List)[p]  <- paste(MeshID[p],"-",MeshTerm[p])
  }
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Mesh to check: ",length(MeshID)," with total number of names: ",length(MeshTerm))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesMesh])
    S = length(sig.genes[sig.genes %in% genesMesh]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(MeshID=character(),
                     MeshTerm=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in c(1:length(MeshID))){
      if (j%%100 == 0) {message("tryingd on MeshID ",j," - ",MeshID[j]," - ",MeshTerm[j])}
      #target = MeshID[j]
      #gENEs = unique(subset(list_Bta, MESHID == target)$GENEID)
      gENEs = DB_List[[j]]
      m = length(total.genes[total.genes %in% gENEs]) # genes from target  and in our dataset
      findG = sig.genes[sig.genes %in% gENEs]
      s = length(findG)
      orig_list = data.frame(Sig_list_out[[i]]) %>% dplyr::filter(ENTREZID_final %in% findG)
      PastefindG = paste(orig_list[,1], collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      #length(gENEs);Pval
      tmp = data.frame(MeshID= MeshID[j],
                       MeshTerm = MeshTerm[j],
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
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("MeshID","MeshTerm", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Mesh_results_b_raw[[i]] = final_raw; names(Mesh_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Mesh raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= Meshthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("MeshID","MeshTerm", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    Mesh_results_b[[i]] = final;names(Mesh_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched Mesh")
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
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_Mesh = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(Mesh_results_b, Mesh_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant MeshIDs found within ",
          length(TestingSubsetNames)," modules/subsets",
          " at the significance level of ",Meshthres)
  message("Nice! - Mesh enrichment finished and data saved")}
