# Joe Boktor
# Caltech - Mazmanian Lab
# 2020

source("src/_load_packages.R")
source("src/_misc_functions.R")
#_______________________________________________________________________________
#      Create Phyloseq Object
#_______________________________________________________________________________

ps_obj <- 
  qza_to_phyloseq(
    features = "data/16S/analysis-11-07-21/table.qza",
    tree = "data/16S/analysis-11-07-21/rooted-tree.qza",
    taxonomy = "data/16S/analysis-11-07-21/taxonomy.qza",
    metadata = "input_files/NHP_16S_metadata_qiime_input.tsv")
ps_obj

sample_data(ps_obj)$read_depth <- sample_sums(ps_obj)
sample_names(ps_obj) <- paste(meta(ps_obj)$monkey, meta(ps_obj)$tissue, sep = " ")

##### Rarefy data ---- 
ps_obj_rare <-  ps_obj %>% phyloseqRarefy()
ps_obj_stool_rare <- ps_obj %>% subset_samples(tissue =="Stool") %>% phyloseqRarefy()
ps_obj_jej_rare <- ps_obj %>% subset_samples(tissue =="Jejunum") %>% phyloseqRarefy()
ps_obj_ile_rare <- ps_obj %>% subset_samples(tissue =="Ileum") %>% phyloseqRarefy()

ps_obj_rare <- 
  list(
  "All" = ps_obj_rare,
  "Stool" = ps_obj_stool_rare,
  "Jejunum" = ps_obj_jej_rare,
  "Ileum" =  ps_obj_ile_rare
)

saveRDS(ps_obj, 
        file = "data/16S/phyloseq.rds")
saveRDS(ps_obj_rare, 
        file = "data/16S/phyloseq_rarified.rds")

