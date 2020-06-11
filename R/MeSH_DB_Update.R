#' Update/download new database - Mesh.
#' New database will be gathered and stored locally.
#'
#' @param keyword Keyword used for naming. Default \code{keyword = "MeSH_DB"}
#' @param DB_location Store path of the new database. Default is current working directory.
#'
#' @return NA. New dataset will be packed and stored in .RData format.
#' @export
#' @import MeSH.db MeSH.Bta.eg.db dplyr
#' @examples MeSH_DB_Update()
#'
MeSH_DB_Update  =function(keyword = "MeSH_DB",DB_location = '.'){
## download DB
ptm <- proc.time();message(paste("Accessing Database..."))
message(paste("Database: ",keyword," download starts!"))
key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
#key_Bta = key_Bta[1:5]
List = MeSHDbi::select(MeSH.db,keys = key_Bta,keytype = "MESHID",columns = c("MESHID","MESHTERM"))
list_Bta = MeSHDbi::select(MeSH.Bta.eg.db,
                           keys = key_Bta,
                           columns = columns(MeSH.Bta.eg.db)[1:3],
                           keytype = "MESHID") %>%
  dplyr::filter(MESHCATEGORY %in% c("D","G")) %>%
  dplyr::left_join(List,by= c("MESHID" = "MESHID")) %>%
  dplyr::select(-MESHCATEGORY)

message("Downloads finished! Time used: ")
print(proc.time() - ptm)
## parse into terms/records (for easier extraction)
DB_List = list_Bta %>%
  dplyr::filter(nchar(list_Bta[,2]) != 0) %>%
  mutate(TermDescription = paste(.[,2],.[,3],sep = "---" )) %>%
  dplyr::select(names(list_Bta)[1],TermDescription) %>%
  split_tibble(column = 'TermDescription',keep = names(list_Bta)[1])
#
MeSH_DB = DB_List
save(MeSH_DB,file = paste0(DB_location,'/',keyword,'.rda'))
#
file_name = paste0(keyword,'.rda')
message(paste("Totally ",length(DB_List),'records were updated on ',Sys.time()))
if (DB_location != '.'){
  message(paste("Database was saved in ",DB_location," in the name of",file_name))
} else {
  pwd = getwd();
  message(paste("Database was saved in ",pwd," in the name of",file_name))
}
load(paste0(DB_location,'/',keyword,'.rda'))
}

