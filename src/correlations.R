
source("src/_load_packages.R")
source("src/_misc_functions.R")

#' This script generates metabolite-to-metabolite correlations within and
#' across tissues

# load data ----
sample_meta <-
  read_xlsx("input_files/misc/sample_keys.xlsx", sheet = "keys")
meta_groups <- sample_meta %>% select(SampleID, group, treatment)
meta_animal_id <- sample_meta %>%
  select(subject_id, SampleID) %>%
  dplyr::rename("Animal" = "subject_id")

df_scaled_imp <-
  readRDS(file = "data/metabolomics/2022-07-05_scaled_imputed_abundance_annotated.rds")
df_all <-
  readRDS(file = "data/metabolomics/2022-07-05_scaled_imputed_abundance.rds")
pathway_map <-
  readRDS(file = "data/metabolomics/2022-07-05_pathway-map.rds")

# _______________________________________________________________________________
#                Metabolite x Metabolite Correlations --------
# _______________________________________________________________________________

df_cor_all <- rcorr(as.matrix(df_all), type = "spearman") %>%
  tidy_cor_output() %>%
  separate(feature_A, c("feature_A", "tissue_A"), sep = "__") %>%
  separate(feature_B, c("feature_B", "tissue_B"), sep = "__")

write.xlsx(pathway_map, "input_files/misc/pathway_map.xlsx")
write.xlsx(
  df_cor_all,
  "data/correlations/metabolite-to-metabolite/2022-07-05_all_metabolite_corrs.xlsx"
)

# _______________________________________________________________________________
#                Metabolite x Cytokine Correlations --------
# _______________________________________________________________________________

# Read trimmed and imputed cytokine profiles
all_cytokines_trimmed_imputed <-
  readRDS(file = "data/cytokines/2022-07-05_cytokine-profiles_imputed.rds")

cytokine_df_long <-
  all_cytokines_trimmed_imputed %>%
  pivot_wider(
    names_from = c("cytokine", "Tissue"),
    values_from = "value",
    names_sep = "__"
  ) %>%
  select(-Animal) %>%
  arrange(match(SampleID, rownames(df_all))) %>%
  column_to_rownames(var = "SampleID")

# Remove all NA columns created by pivot wider
cytokine_df_long <- cytokine_df_long[, colSums(is.na(cytokine_df_long)) < nrow(cytokine_df_long)] # Remove rows with NA only
cytokine_df_long %>% as_tibble()

# Parallelized calculation of correlations
cytokine_corrs <-
  corr_loop_parallel(cytokine_df_long, df_all, obj.name = "Cytokine__NoStimulant")

cytokine_corrs_FDR <- cytokine_corrs %>%
  corr_feature_tissue_FDR() %>%
  select(-object_name)

cyt_to_mets <- "data/correlations/cytokine-to-metabolite/"
write.xlsx(
  cytokine_corrs_FDR,
  glue(
    "{cyt_to_mets}/2022-07-06_cytokine-metabolite_correlations.xlsx"
  ),
  overwrite = TRUE
)

saveRDS(cytokine_corrs_FDR,
  file = glue("{cyt_to_mets}/2022-07-06_cytokine-metabolite_correlations.rds")
)

cytokine_corrs_excel_out <- list()
cytokine_corrs_excel_out[['Plasma_cytokine_corrs']] <-
  cytokine_corrs_FDR %>%
  filter(tissue_A %in% c("Plasma", "Plasma-PBMCs"))
cytokine_corrs_excel_out[['Brain_cytokine_corrs']] <-
  cytokine_corrs_FDR %>%
  filter(tissue_A %in% c("CSF", "Cortex", "HIPP", "AC"))
cytokine_corrs_excel_out[['GI_cytokine_corrs']] <-
  cytokine_corrs_FDR %>%
  filter(tissue_A %in% c("Colon-PBMCs", "Ileum-PBMCs", "Jejunum-PBMCs", "MLN-PBMCs"))

write.xlsx(
  cytokine_corrs_excel_out,
  glue(
    "{cyt_to_mets}/2022-07-06_cytokine-metabolite_correlations_by_tissue.xlsx"
  ),
  overwrite = TRUE
)

# _______________________________________________________________________________
#                Metabolite x Behavior Correlations --------
# _______________________________________________________________________________

behavior_data <- readRDS(file = "data/behaviors/2022-07-05_behaviors.rds")
df_behavior <- behavior_data$df_behavior

behavior_all_mets <- corr_loop_parallel(df_behavior, df_all, obj.name = "Behavior")
behavior_all_mets_fdr <- corr_heatmap_FDRv3_behavior(behavior_all_mets)
heatmap_behavior_all_mets <- corr_heatmap_facet_behav(behavior_all_mets_fdr)

ggsave(
  heatmap_behavior_all_mets,
  filename = "figures/correlations/heatmap_behavior_all_mets.svg",
  width = 8,
  height = 6.5
)

behavior_corrs <- list(
  "behaviors_allMetabolites" = behavior_all_mets_fdr
)
write.xlsx(
  behavior_corrs,
  "data/correlations/behavior_metabolite/Behavior_Metabolite_Correlations.xlsx"
)
