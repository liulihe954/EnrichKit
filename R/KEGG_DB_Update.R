#' Update/download new database - KEGG.
#' New database will be gathered and stored locally.
#'
#' @param keyword Keyword used for naming. Default \code{keyword = "KEGG_DB"}
#' @param DB_location Store path of the new database. Default is current working directory.
#' @param species Default \code{species = "bta"}
#' @param id.type Default \code{id.type = "kegg"}
#'
#' @return NA. New dataset will be packed and stored in .RData format.
#' @export
#' @import gage dplyr magrittr
#' @examples Kegg_DB_Update()
#'
Kegg_DB_Update  = function(species = "bta",
                           id.type = "kegg",
                           keyword = "KEGG_DB",
                           DB_location = "."){
  #
  ptm <- proc.time();message(paste("Database: ",keyword," download starts!"))
  sdb = kegg.gsets(species = species, id.type = id.type, check.new = F)
  kegg.gs = sdb$kg.sets[sdb$sigmet.id]
  message("Downloads finished! Time used: ")
  print(proc.time() - ptm)
  #
  KEGG_DB =  kegg.gs
  save(KEGG_DB,file = paste0(DB_location,'/',keyword,'.rda'))
  #
  pwd = getwd();file_name = paste0(keyword,'.rda')
  message(paste("Totally ",length( kegg.gs),'records were updated on ',Sys.time()))
  message(paste("Database was saved in",pwd," in the name of",file_name))
  load(paste0(DB_location,'/',keyword,'.rda'))
}
