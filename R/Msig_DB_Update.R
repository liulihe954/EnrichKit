#' Update/download new database - Msig.
#' New database will be gathered and stored locally.
#' Note this function could take quite a while...
#'
#' @param keyword Keyword used for naming. Default \code{keyword = "Msig_DB"}
#' @param DB_location Store path of the new database. Default is current working directory.
#'
#' @return NA. New dataset will be packed and stored in .RData format.
#' @export
#' @import msigdbr dplyr
#' @examples Msig_DB_Update()
#'
Msig_DB_Update  =function(keyword = "Msig_DB",DB_location = '.'){

  # specify database location
  url_template = 'http://software.broadinstitute.org/gsea/msigdb/cards/'
  # get pathway name index
  m_df = msigdbr(species = "Bos taurus")

  # obtain name index and paste to all urls
  Msig_name_index = unique(m_df$gs_name)
  Msig_urls = paste(url_template,Msig_name_index,sep = "")

  # prepare R function to retrive description
  Get_Descrip = function(URL,
                         selector = "td",
                         pos = 6){
    library(rvest)
    base_url <- URL
    webpage <- read_html(base_url)
    # Get the artist name
    target_raw <- html_nodes(webpage,selector)
    target_raw <- as.character(html_text( target_raw))
    text = str_replace_all( target_raw,"[\r\n]","")
    final = as.character(text[pos])
    message('try on ', URL)
    return(final)
  }
  #
  totalterm <- length(Msig_name_index)
  pb <- txtProgressBar(min = 0, # create progress bar
                       max = totalterm,
                       style = 3)
  # loop to retrive
  All_Descrip = c()
  for (i in seq_along(Msig_urls)){
    setTxtProgressBar(pb, i)
    All_Descrip[i] = Get_Descrip(Msig_urls[i])
  }
  close(pb)
  # massage
  Universe_Descrip = data.frame(cbind(Msig_name_index,All_Descrip))
  colnames(Universe_Descrip) = c('gs_name','gs_description')
  #
  m_df_all = m_df %>%
    #dplyr::filter(gs_name %in% Msig_name_index) %>%
    dplyr::left_join(Universe_Descrip,
                     by = c("gs_name" = "gs_name")) %>%
    dplyr::select(entrez_gene,gs_id,gs_description) %>%
    data.frame() %>%
    mutate(TermDescription = paste(.[,2],.[,3],sep = "---" )) %>%
    dplyr::select(entrez_gene,TermDescription)

  DB_List = m_df_all %>%
    split_tibble(column = 'TermDescription',keep = names(.)[1])
  #
  Msig_DB = DB_List
  save(Msig_DB,file = paste0(DB_location,'/',keyword,'.rda'))
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
