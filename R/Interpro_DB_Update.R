#' Update/download new database - Interpro.
#' New database will be gathered and stored locally.
#'
#' @param biomart BioMart database name you want to connect to. Default \code{biomart="ensembl"}
#' @param dataset Mart object created with the useMart function. Default \code{dataset="btaurus_gene_ensembl"}
#' @param attributes Attributes you want to retrieve. Refer to \code{listAttributes}.
#' @param keyword Keyword used for naming. Default \code{keyword = "Interpro_DB"}
#' @param DB_location Store path of the new database. Default is current working directory.
#'
#' @return NA. New dataset will be packed and stored in .RData format.
#' @export
#' @import biomaRt
#' @examples Interpro_DB_Update()
#'
Interpro_DB_Update = function(biomart="ensembl",
                              dataset="btaurus_gene_ensembl",
                              attributes = c("ensembl_gene_id","interpro","interpro_description"),
                              keyword = "Interpro_DB",
                              DB_location = '.'){
  library(ggplot2);library(biomaRt);library(gage);library(magrittr)# load pkg
  ## download DB
  ptm <- proc.time();message(paste("Accessing BiomaRt..."));message(paste("Database: ",keyword," download starts!"))
  database = useMart(biomart)
  genome = useDataset(dataset, mart = database)
  gene = getBM(attributes,mart = genome)
  message("Downloads finished! Time used: ")
  print(proc.time() - ptm)
  ## parse into terms/records (for easier extraction)
  DB_List = gene %>%
    dplyr::filter(nchar(gene[,2]) != 0) %>%
    mutate(TermDescription = paste(.[,2],.[,3],sep = "---" )) %>%
    dplyr::select(names(gene)[1],TermDescription) %>%
    split_tibble(column = 'TermDescription',keep = names(gene)[1])
  #
  Interpro_DB = DB_List
  save(Interpro_DB,file = paste0(DB_location,'/',keyword,'.rda'))
  #
  file_name = paste0(keyword,'.rda')
  message(paste("Totally ",length(DB_List),'records were updated on ',Sys.time()))
  if (DB_location != '.'){
    message(paste("Database was saved in ",DB_location," in the name of",file_name))
  } else {
    pwd = getwd()
    message(paste("Database was saved in ",pwd," in the name of",file_name))
  }
  load(paste0(DB_location,'/',keyword,'.rda'))
}
