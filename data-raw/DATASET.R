library(tidyverse)
library(org.Bt.eg.db)


# Obtain Universe Bta Identifier
Universe_id_bta = AnnotationDbi::select(org.Bt.eg.db,
                                        as.character(AnnotationDbi::keys(org.Bt.eg.db,keytype = c("ENSEMBL"))),
                                        columns = c("ENTREZID","SYMBOL"),
                                        keytype = "ENSEMBL") %>%
  dplyr::distinct(ENSEMBL,.keep_all= TRUE)
save(Universe_id_bta, file = "data/Universe_id_bta.rda")
save(Universe_id_bta, file = "data-raw/Universe_id_bta.rda")


# obtain GO

