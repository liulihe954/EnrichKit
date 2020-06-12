# set wd
setwd('./data-raw/')

# ============================== #
# Obtain Universe Bta Identifier #
# ============================= #
Universe_id_bta = AnnotationDbi::select(org.Bt.eg.db,
                                        as.character(AnnotationDbi::keys(org.Bt.eg.db,keytype = c("ENSEMBL"))),
                                        columns = c("ENTREZID","SYMBOL"),
                                        keytype = "ENSEMBL") %>%
  dplyr::distinct(ENSEMBL,.keep_all= TRUE)

#save(Universe_id_bta, file = "data/Universe_id_bta.rda")

save(Universe_id_bta, file = "data-raw/Universe_id_bta.rda")

# ======= #
# == GO ==#
# ======= #
Go_DB_Update = function(biomart="ensembl",
                        dataset="btaurus_gene_ensembl",
                        attributes = c("ensembl_gene_id","go_id","name_1006"),
                        keyword = "GO_DB",
                        DB_location = '.'){
  library(tidyverse);library(biomaRt);library(gage);library(magrittr)# load pkg
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
  GO_DB = DB_List

  #save(GO_DB,file = paste0(DB_location,'/',keyword,'.rda'))

  save(GO_DB,file = paste0('../data/',keyword,'.rda'))
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
# Go_DB_Update()


# ========== #
# == KEGG  ==#
# ========= #
Kegg_DB_Update  = function(species = "bta",
                           id.type = "kegg",
                           keyword = "KEGG_DB",
                           DB_location = "."){
  library(biomaRt);library(gage);library(magrittr);library(tidyverse)
  #
  ptm <- proc.time();message(paste("Database: ",keyword," download starts!"))
  sdb = kegg.gsets(species = species, id.type = id.type, check.new = F)
  kegg.gs = sdb$kg.sets[sdb$sigmet.id]
  message("Downloads finished! Time used: ")
  print(proc.time() - ptm)
  #
  KEGG_DB =  kegg.gs

  #save(KEGG_DB,file = paste0(DB_location,'/',keyword,'.rda'))

  save(GO_DB,file = paste0('../data/',keyword,'.rda'))
  #
  pwd = getwd();file_name = paste0(keyword,'.rda')
  message(paste("Totally ",length( kegg.gs),'records were updated on ',Sys.time()))
  message(paste("Database was saved in",pwd," in the name of",file_name))
  load(paste0(DB_location,'/',keyword,'.rda'))
}

# Kegg_DB_Update()


# ============== #
# == Interpro  =#
# ============= #
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

  #save(Interpro_DB,file = paste0(DB_location,'/',keyword,'.rda'))

  save(GO_DB,file = paste0('../data/',keyword,'.rda'))
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
# Interpro_DB_Update()


# ============== #
# == Reactome  =#
# ============= #
Reactome_DB_Update  =function(websource = "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt",
                              Species = "Bos taurus",
                              keyword = "Reactome_DB",
                              DB_location = '.'){
  library(data.table);library(biomaRt);library(tidyverse) #load pkg
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
  #
  #save(Reactome_DB,file = paste0(DB_location,'/',keyword,'.rda'))

  save(GO_DB,file = paste0('../data/',keyword,'.rda'))
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

# Reactome_DB_Update()

# ============ #
# ==  MeSH  =#
# =========== #
MeSH_DB_Update  =function(keyword = "MeSH_DB",DB_location = '.'){
  library(MeSH.db);library(MeSH.Bta.eg.db);library(tidyverse);library(gage);library(magrittr)
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

  #save(MeSH_DB,file = paste0(DB_location,'/',keyword,'.rda'))

  save(GO_DB,file = paste0('../data/',keyword,'.rda'))
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
# MeSH_DB_Update()
# ============ #
# ==   Msig   =#
# =========== #
Msig_DB_Update  =function(keyword = "Msig_DB",DB_location = '.'){

  # specify database location
  url_template = 'http://software.broadinstitute.org/gsea/msigdb/cards/'
  # get pathway name index
  library(tidyverse);library(msigdbr)
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

  # loop to retrive
  All_Descrip = c()
  for (i in seq_along(Msig_urls)){
    #message(i)
    All_Descrip[i] = Get_Descrip(Msig_urls[i])
  }
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

  #save(Msig_DB,file = paste0(DB_location,'/',keyword,'.rda'))

  save(GO_DB,file = paste0('../data/',keyword,'.rda'))
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
#Msig_DB_Update()

