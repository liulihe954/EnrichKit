#' Match different identifiers to gene sets provided
#'
#' @param GeneSetNames Sigle or multiple gene set names, put into a vector if multiple. e.g. c('Lactation1','Lactation2','Lactation3')
#' @param SigGene_list An list object; where each element contains signifcant genes. list dimension must match \code{TotalGene_list}
#' @param TotalGene_list An list object; where each element contains total genes. list dimension must match \code{SigGene_list}
#' @param IDtype Gene Identifier type, choose from \code{c("ens","entrez","symbol")}
#'
#' @return An list object with same dimension as SigGene_list and Totalgene_list, where each element has comprehensive id information.
#' @export
#' @import dplyr HGNChelper
#' @examples
#'
convertNformatID = function(GeneSetNames,
                            SigGene_list,
                            TotalGene_list,
                            IDtype = "ens"){
  ## error handling
  if(length(SigGene_list) != length(TotalGene_list)){
    stop("SigGene_list must have same length as TotalGene_list. Try re-organize.")
  }
  #
  if (length(GeneSetNames) > length(SigGene_list)){
    warning("Extra gene set names provided, last ",length(GeneSetNames)-length(SigGene_list)," ones will be dropped!")
  } else if (length(GeneSetNames) < length(SigGene_list)){
    warning("Insufficient gene set names provided, thus UNKNOWN will be assigned to the last ",length(SigGene_list)-length(GeneSetNames)," sets!")
    GeneSetNames[(length(GeneSetNames)+1):length(SigGene_list)] = paste0('UNKNOWN',seq(1:(length(SigGene_list)-length(GeneSetNames))))
  }
  #
  overlap_tester = c()
  for (i in seq_along(GeneSetNames)){
    overlap_tester[i] = all(SigGene_list[[i]] %in% TotalGene_list[[i]])
  }
  if(all(overlap_tester) != TRUE){
    errorloc = which(overlap_tester == F)
    stop("Pair ",errorloc," has non-overlap sig-total gene set\n")
  }

  ## formatting
  data(Universe_id_bta)
  IDtype_all = c('ens','entrez','symbol')
  IDprovd = which(IDtype_all == IDtype)
  sigGene_All = list()

  for (i in seq_along(GeneSetNames)){

    siglist = unlist(SigGene_list[[i]],use.names = F)
    totallist = unlist(TotalGene_list[[i]],use.names = F)

    Gather = data.frame(Gene = totallist,stringsAsFactors=FALSE) %>%
      dplyr::left_join(Universe_id_bta,by=c('Gene'=names(Universe_id_bta)[IDprovd])) %>%
      mutate(SYMBOL_Suggested = check_symbol(SYMBOL)) %>%
      mutate(Sig = ifelse(Gene %in% siglist,'1','0'))
    sigGene_All[[i]] = Gather
    names(sigGene_All)[i] = GeneSetNames[i]
  }
  return(sigGene_All)
}






