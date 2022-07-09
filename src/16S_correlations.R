# Joe Boktor
# Caltech - Mazmanian Lab
# 2020

source("src/_load_packages.R")
source("src/_misc_functions.R")
source("src/_16S_functions.R")
library(progress)

# Load 16S data
ps_obj <- readRDS(file = "data/16S/phyloseq.rds")
# Load metabolome data
df_all <- readRDS(file = "input_files/Metabolomics/scaled_&_imputed_abundance.rds")
# Load Behavior data
behaviors <- readRDS(file = "input_files/Behavior_Immune_data/behaviors.rds")
# Load GI cytokine data
gi_cytokines <- readRDS("input_files/Behavior_Immune_data/GI_cytokines_all.rds")


# Prep 16S abundance data
ps_df_stool <- ps_obj %>%
  subset_samples(tissue == "Stool") %>%
  prep16SAbundTable()
ps_df_stool %>% dim()
ps_df_jej <- ps_obj %>%
  subset_samples(tissue == "Jejunum") %>%
  prep16SAbundTable()
ps_df_jej %>% dim()
ps_df_ile <- ps_obj %>%
  subset_samples(tissue == "Ileum") %>%
  prep16SAbundTable()
ps_df_ile %>% dim()

# Prep metabolomics  data
stool_mets <- df_all %>%
  select(contains("__Feces"))
jejunum_mets <- df_all %>%
  select(contains("__Jejunum"))
ileum_mets <- df_all %>%
  select(contains("__Ileum"))

# Prep GI cytokine data
colon_cytokines <- list()
jejunum_cytokines <- list()
ileum_cytokines <- list()
for (condition in names(gi_cytokines)) {
  print(condition)
  colon_cytokines[[condition]] <- gi_cytokines[[condition]] %>%
    select(contains("__colon"))
  jejunum_cytokines[[condition]] <- gi_cytokines[[condition]] %>%
    select(contains("__Jejunum"))
  ileum_cytokines[[condition]] <- gi_cytokines[[condition]] %>%
    select(contains("__Ileum"))
}

# Select NHPs with paired 16S and metabolome samples and make row-order identical
stool_mets_input <- stool_mets[rownames(ps_df_stool), ]
cat("Checking that sample ids are aligned: ", all(rownames(ps_df_stool) == rownames(stool_mets_input)), "\n")
stool_mets_input %>% dim()
jejunum_mets_input <- jejunum_mets[rownames(ps_df_jej), ]
cat("Checking that sample ids are aligned: ", all(rownames(ps_df_jej) == rownames(jejunum_mets_input)), "\n")
jejunum_mets_input %>% dim()
ileum_mets_input <- ileum_mets[rownames(ps_df_ile), ]
cat("Checking that sample ids are aligned: ", all(rownames(ps_df_ile) == rownames(ileum_mets_input)), "\n")
ileum_mets_input %>% dim()

# Select NHPs with paired 16S and PBMC cytokine gut samples and make row-order identical
colon_cyt_input <- list()
jejunum_cyt_input <- list()
ileum_cyt_input <- list()

for (condition in names(gi_cytokines)) {
  colon_cyt_input[[condition]] <- colon_cytokines[[condition]][rownames(ps_df_stool), ]
  cat("Checking that sample ids are aligned: ", all(rownames(ps_df_stool) == rownames(colon_cyt_input[[condition]])), "\n")
  jejunum_cyt_input[[condition]] <- jejunum_cytokines[[condition]][rownames(ps_df_jej), ]
  cat("Checking that sample ids are aligned: ", all(rownames(ps_df_jej) == rownames(jejunum_cyt_input[[condition]])), "\n")
  ileum_cyt_input[[condition]] <- ileum_cytokines[[condition]][rownames(ps_df_ile), ]
  cat("Checking that sample ids are aligned: ", all(rownames(ps_df_ile) == rownames(ileum_cyt_input[[condition]])), "\n")
}


# _______________________________________________________________________________

# microbe metabolite correlations ----
microbe_mets_stool <- corr_loop_parallel(ps_df_stool, stool_mets_input, obj.name = "16S_Metabolome__Stool")
microbe_mets_jejunum <- corr_loop_parallel(ps_df_jej, jejunum_mets_input, obj.name = "16S_Metabolome__Jejunum")
microbe_mets_ileum <- corr_loop_parallel(ps_df_ile, ileum_mets_input, obj.name = "16S_Metabolome__Ileum")

microbe_metabolite_cors <-
  bind_rows(
    microbe_mets_stool %>% corr_FDR_16S(),
    microbe_mets_jejunum %>% corr_FDR_16S(),
    microbe_mets_ileum %>% corr_FDR_16S()
  )
saveRDS(microbe_metabolite_cors, file = "data/correlations/microbe-to-metabolite/Microbe_Metabolite_Correlations.rds")
write.xlsx(microbe_metabolite_cors,
  file = "data/correlations/microbe-to-metabolite/Microbe_Metabolite_Correlations.xlsx",
  overwrite = T
)

# microbe cytokine correlations ----
microbe_cyt_stool <- list()
microbe_cyt_jejunum <- list()
microbe_cyt_ileum <- list()
iterlength <- length(names(gi_cytokines))
pb <- progress_bar$new(
  format = "  Calculating Correlations [:bar] :current/:total :percent in :elapsed eta::eta",
  total = iterlength, clear = FALSE, width = 90
)
for (condition in names(gi_cytokines)) {
  pb$tick()
  microbe_cyt_stool[[condition]] <-
    corr_loop_parallel(ps_df_stool, colon_cyt_input[[condition]],
      obj.name = paste0("16S_Cytokines__", condition, "__StoolxColon")
    )
  microbe_cyt_jejunum[[condition]] <-
    corr_loop_parallel(ps_df_jej, jejunum_cyt_input[[condition]],
      obj.name = paste0("16S_Cytokines__", condition, "__Jejunum")
    )
  microbe_cyt_ileum[[condition]] <-
    corr_loop_parallel(ps_df_ile, ileum_cyt_input[[condition]],
      obj.name = paste0("16S_Cytokines__", condition, "__Ileum")
    )
}

# add q-value and collapse lists
microbe_cyt_stool_df <- map(microbe_cyt_stool, corr_FDR_16S_cyt) %>% map_df(as.data.frame)
microbe_cyt_jejunum_df <- map(microbe_cyt_jejunum, corr_FDR_16S_cyt) %>% map_df(as.data.frame)
microbe_cyt_ileum_df <- map(microbe_cyt_ileum, corr_FDR_16S_cyt) %>% map_df(as.data.frame)

microbe_cytokine_cors <-
  bind_rows(
    microbe_cyt_stool_df,
    microbe_cyt_jejunum_df,
    microbe_cyt_ileum_df
  )
saveRDS(microbe_cytokine_cors, file = "data/correlations/microbe-to-metabolite/Microbe_Cytokine_Correlations.rds")
write.xlsx(microbe_cytokine_cors,
  file = "data/correlations/microbe-to-metabolite/Microbe_Cytokine_Correlations.xlsx",
  overwrite = T
)
