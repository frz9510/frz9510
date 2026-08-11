rm(list = ls());

library(tidyr) 

library(dplyr)
setwd('/lustre09/project/6033529/ARSACS/results/ht-Farnaz')
getwd()

## Your ht file with hgs
samples_file <-'patients_carrying_0-Copy1.8MBhaplotypes001.asc'
samples_df <- read.table(samples_file, header = T, sep ="")

df_new <- samples_df %>%
  mutate(
    # build the prefix from the first 6 chars
    prefix = substr(GeneID1_trimmed, 1, 6),
    # create full label
    new_value = paste0(prefix, "_", GeneID1_trimmed)
  ) %>%
  # duplicate rows, one with _1 and one with _2
  rowwise() %>%
  mutate(
    expanded = list(c(paste0(new_value, "_1"),
                      paste0(new_value, "_2")))
  ) %>%
  unnest(expanded) %>%
  select(expanded)

print(df_new)

myfile<-paste('/lustre09/project/6033529/ARSACS/results/ht-Farnaz/patients_to_generate_16mb_raxmltree.asc')
write.table(df_new,myfile,quote=FALSE,row.names = FALSE)




