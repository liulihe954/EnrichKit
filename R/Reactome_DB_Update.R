#' Update/download new database - Reactome.
#' New database will be gathered and stored locally.
#'
#' @param keyword Keyword used for naming. Default \code{keyword = "Reactome_DB"}
#' @param DB_location Store path of the new database. Default is current working directory.
#' @param websource \code{websource = "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt"}
#' @param Species Default \code{Species = "Bos taurus"}
#'
#' @return NA. New dataset will be packed and stored in .RData format.
#' @export
#' @import biomaRt dplyr data.table
#' @examples Reactome_DB_Update()
#'
Reactome_DB_Update  =function(websource = "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt",
                              Species = "Bos taurus",
                              keyword = "Reactome_DB",
                              DB_location = '.'){
  ## download DB
  ptm <- proc.time();message(paste("Accessing",websource,'...'))
  message(paste("Database: ",keyword," download starts!"))
  entrezReactome_DB_all_path_bt = fread(websource) %>%
    dplyr::filter(V6 == Species) %>%
    dplyr::select(V1,V2,V4) %>%
    dplyr::rename(EntrezID = V1,
                  ReactomeID = V2,
                  Reactome_Description = V4) %>%
    arrange(ReactomeID)
  message("Downloads finished! Time used: ")
  print(proc.time() - ptm)
  ## parse into terms/records (for easier extraction)
  DB_List = entrezReactome_DB_all_path_bt %>% data.frame() %>%
    dplyr::filter(nchar(entrezReactome_DB_all_path_bt[,2]) != 0) %>%
    mutate(TermDescription = paste(.[,2],.[,3],sep = "---" )) %>%
    dplyr::select(names(entrezReactome_DB_all_path_bt)[1],TermDescription) %>%
    split_tibble(column = "TermDescription",keep = names(entrezReactome_DB_all_path_bt)[1])
  #
  Reactome_DB = DB_List
  save(Reactome_DB,file = paste0(DB_location,'/',keyword,'.rda'))
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
