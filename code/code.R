## ENRICHMENTS
## Pulling in Gen-sequences

## KFP 2022-06-13



# load packages -----------------------------------------------------------
library(tidyverse)
#library(ape)
#BiocManager::install("rtracklayer")
#library(rtracklayer)

#
# load data files ---------------------------------------------------------
## these are all the csv files in the target directory

# 1. pull all the csv files 
FILEPATH = "bracky_bsub_analysis"
filePaths <- list.files(path = FILEPATH, pattern = ".csv", full.names = TRUE)

data_brach = 
  do.call(rbind, lapply(filePaths, function(path) {
    df <- read.csv(path)
    # then add a new column `source` to denote the file name
    df[["source"]] <- rep(path, nrow(df))
    df}))


#
# load sequence files -----------------------------------------------------
## these are the xlsm files, can be imported using readxl::read_xlsx

filePaths_seq <- list.files(path = FILEPATH, pattern = ".xlsm", full.names = TRUE)

data_seq = 
  do.call(bind_rows, lapply(filePaths_seq, function(path) {
    df <- readxl::read_xlsx(path,
                            # set column names
                            col_names = c("sequence", "refseq", "region", "a", "b", "c", "d", "e", "ID"))
    # then add a new column `source` to denote the file name
    df[["source"]] <- rep(path, nrow(df))
    df}))

## clean the sequence data
## from the ID column, we want to extract only Genbank, product, and gbkey info

data_seq_cleaned = 
  data_seq %>% 
  mutate(Genbank = str_extract_all(ID, "Genbank:.+?;"),
         product = str_extract_all(ID, "product=.+?;"),
         gbkey = str_extract_all(ID, "gbkey=.+?;")) %>% 
  distinct(Genbank, product, gbkey) %>% 
  mutate_all(as.character) %>% 
  mutate(across(where(is.character), ~na_if(., "character(0)"))) %>% 
  mutate(Genbank = str_remove(Genbank, "Genbank:"),
         Genbank = str_remove(Genbank, ";"),
         product = str_remove(product, "product="),
         product = str_remove(product, ";"),
         gbkey = str_remove(gbkey, "gbkey="),
         gbkey = str_remove(gbkey, ";"))

#
# merge the data with the sequence data -----------------------------------

data_merged = 
  left_join(data_brach, data_seq_cleaned, by = c("names" = "Genbank"))

#
# export the merged file --------------------------------------------------
data_merged %>% write.csv("processed/data/combined_brachy_sequences.csv", row.names = FALSE, na = "")

