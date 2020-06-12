#
#' Over-representation test process
#'
#' @param GeneSet A list object; Ideal input is object resulted from \code{\link{convertNformatID}}
#' Each element is a gene set to be tested, should have different identifiers in columns, 1st column - ensembl id; 2nd - entrezid,
#' Also should have a column "Sig" indicating significant or not, 1 - significant; 0 insignificant
#'
#' @param Database Indicate one of the dataset to be used; choose from c("go","kegg","interpro","mesh","msig","reactome")
#' @param minOverlap Minimum overlap number of total genes and genes of a target pathway, default value is \code{4}
#' @param pvalue_thres Pvalue threshold of the Fisher's exact test. Refer to \code{\link{fisher.test}}
#' @param adj_pvalue_thres Adjusted value threshold for multiple testing control.
#' @param padj_method Adjusted pvalues; choose methods from c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"). Refer to \code{\link{p.adjust}}
#' @param NewDB If new database will be used, please do \code{NewDB = T}. Please make sure xxx.rda was generated using XXX_DB_Update() and is currently in the working directory.
#'
#' @return
#' @export
#' @import dplyr
#'
#' @examples
#'
HyperGEnrich = function(GeneSet = sigGene_All,
                        Database = '', #'go','kegg,'interpro','mesh','msig','reactome'
                        minOverlap = 4,
                        pvalue_thres = 0.05,
                        adj_pvalue_thres = 1,
                        padj_method = "BH", #c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
                        NewDB = F){
  # get db
  TestingSubsetNames = names(GeneSet)
  message("Total Number of subsets/module to check: ",length(TestingSubsetNames))
  if (NewDB){
    if (Database == 'go'){
      DB_List = get(load('./GO_DB.rda'));IDtype = 1
    } else if (Database == 'kegg'){
      DB_List = get(load('./KEGG_DB.rda'));IDtype = 2
    } else if (Database == 'interpro'){
      DB_List = get(load('./Interpro_DB.rda'));IDtype = 1
    } else if (Database == 'mesh'){
      DB_List = get(load('./MeSH_DB.rda'));IDtype = 2
    } else if (Database == 'msig'){
      DB_List = get(load('./Msig_DB.rda'));IDtype = 2
    } else if (Database == 'reactome'){
      DB_List = get(load('./Reactome_DB.rda'));IDtype = 2
    } else {
      print('Please choose database: go,kegg,interpro,mesh,msig,reactome')
    }
  } else {
    if (Database == 'go'){
      DB_List = get(data(GO_DB));IDtype = 1
    } else if (Database == 'kegg'){
      DB_List = get(data(KEGG_DB));IDtype = 2
    } else if (Database == 'interpro'){
      DB_List = get(data(Interpro_DB));IDtype = 1
    } else if (Database == 'mesh'){
      DB_List = get(data(MeSH_DB));IDtype = 2
    } else if (Database == 'msig'){
      DB_List = get(data(Msig_DB));IDtype = 2
    } else if (Database == 'reactome'){
      DB_List = get(data(Reactome_DB));IDtype = 2
    } else {
      print('Please choose database: go,kegg,interpro,mesh,msig,reactome')
    }
  }
  #
  GeneInDB = unique(unlist(DB_List,use.names = F))
  message('Database selected: ', Database,'.',
          '\nTotal number of terms to check ', length(DB_List),'.',
          '\nTotal number of genes showed up in any record ',length(GeneInDB),'.')
  #
  total_enrich = 0
  raw_pvalue_all = numeric()
  results = list()
  results_raw = list()

  #
  for (i in c(1:(length(TestingSubsetNames)))){
    Genelist_tmp = GeneSet[[i]]
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    # get sig gene
    sig.genes_tmp = Genelist_tmp %>% data.frame() %>%
      dplyr::filter(Sig == '1') %>%
      dplyr::select(names(GeneSet[[i]])[IDtype]) %>%
      unlist(use.names = F)
    sig.genes = sig.genes_tmp %>% na.omit()
    # get total gene
    total.genes_tmp = Genelist_tmp %>% data.frame() %>%
      dplyr::filter(Sig %in% c('1','0')) %>%
      dplyr::select(names(GeneSet[[i]])[IDtype]) %>%
      unlist(use.names = F)
    total.genes = total.genes_tmp %>% na.omit()

    message('Total genes: ',length(total.genes),'. with ',length(total.genes_tmp)-length(total.genes),' loss due to identifier match',
            '\nSig Genes: ', length(sig.genes),'. with ',length(sig.genes_tmp)-length(sig.genes) ,' loss due to identifier match')

    # overlap with the database
    N = length(total.genes[total.genes %in% GeneInDB])
    S = length(sig.genes[sig.genes %in% GeneInDB]) #

    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    # formatting out put
    out = data.frame(Term=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    #
    totalterm <- length(DB_List)
    pb <- txtProgressBar(min = 0, # create progress bar
                         max = totalterm,
                         style = 3)

    for(j in seq_along(DB_List)){
      setTxtProgressBar(pb, j)
      #if (j%%1000 == 0) {message("tryingd on term ",j," - ",names(DB_List)[j])}
      gENEs = DB_List[[j]] # all gene in target GO #### note
      m = length(total.genes[total.genes %in% gENEs]) # genes from target GO and in our dataset
      findG = sig.genes[sig.genes %in% gENEs]
      s = length(findG)
      PastefindG = paste(findG,collapse="/")
      matrix(c(1:4),byrow = 2, nrow = 2)
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(Term= names(DB_List)[j],
                       totalG = m,
                       sigG = s,
                       pvalue = Pval,
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig,
                       findG = PastefindG)
      out = rbind(out,tmp)
    }
    close(pb)
    names(out) = c('Term','totalG','sigG','pvalue','ExternalLoss_total','ExternalLoss_sig','findG')

    # raw
    "/-" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y))
    final_raw = out %>%
      arrange(pvalue) %>%
      dplyr::mutate(hitsPerc = sigG*100 /- totalG) %>%
      mutate(adj.pvalue = p.adjust(pvalue, method = padj_method))

    results_raw[[i]] = final_raw
    names(results_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"terms tested - raw")

    # selection
    final = final_raw %>%
      dplyr::filter(totalG >= minOverlap) %>%
      dplyr::filter(pvalue <= pvalue_thres) %>%
      dplyr::filter(`adj.pvalue` <= adj_pvalue_thres)

    results[[i]] = final
    names(results)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched terms - selected")

    # selection ends
    message("Significant Enrichment Hits: ",nrow(final))
    total_enrich = total_enrich + nrow(final)
  }
  #
  save(results_raw,results,
       file = paste0(Database,'-Enrichment-',minOverlap,'-',pvalue_thres,'-',adj_pvalue_thres,'.rda'))

  message(total_enrich," significant terms found within ",
          length(TestingSubsetNames)," modules/subsets ",
          "with pvalue and adj.pvalue set at ", pvalue_thres,' and ',adj_pvalue_thres)
}
