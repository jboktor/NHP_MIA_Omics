
source("src/_load_packages.R")
source("src/_misc_functions.R")

#' This script generates metabolite-to-metabolite correlations within and 
#' across tissues

sample_meta <-
  read_xlsx('input_files/misc/sample_keys.xlsx', sheet = 'keys')
meta_groups <- sample_meta %>% select(SampleID, group, treatment)
meta_animal_id <- sample_meta %>%
  select(subject_id, SampleID) %>%
  dplyr::rename("Animal" = "subject_id")

df_scaled_imp <-
  readRDS(file = "data/metabolomics/2022-07-05_scaled_imputed_abundance_annotated.rds")
df_all <-
  readRDS(file = "data/metabolomics/2022-07-05_scaled_imputed_abundance.rds")
pathway_map <- df_scaled_imp %>%
  select(BIOCHEMICAL, `SUPER PATHWAY`, `SUB PATHWAY`) %>%
  distinct()
#_______________________________________________________________________________
#                Metabolite x Metabolite Correlations --------
#_______________________________________________________________________________

df_cor_all <- rcorr(as.matrix(df_all), type = "spearman") %>% 
  tidy_cor_output() %>% 
  separate(feature_A, c("feature_A", "tissue_A"), sep = "__") %>% 
  separate(feature_B, c("feature_B", "tissue_B"), sep = "__")

write.xlsx(pathway_map, 'input_files/misc/pathway_map.xlsx')
write.xlsx(
  df_cor_all,
  'data/correlations/metabolite-to-metabolite/2022-07-05_all_metabolite_corrs.xlsx'
)

#_______________________________________________________________________________
#                Metabolite x Cytokine Correlations --------
#_______________________________________________________________________________

# Read trimmed and imputed cytokine profiles
all_cytokines_trimmed_imputed <-
  readRDS(file = 'data/cytokines/2022-07-05_cytokine-profiles_imputed.rds')

cytokine_df_long <-
  all_cytokines_trimmed_imputed %>%
  pivot_wider(
    names_from = c('cytokine', 'Tissue'),
    values_from = 'value',
    names_sep = '__' 
  ) %>% 
  select(-Animal) %>%
  arrange(match(SampleID, rownames(df_all))) %>%
  column_to_rownames(var = "SampleID")

# Remove all NA columns created by pivot wider
cytokine_df_long <- cytokine_df_long[ , colSums(is.na(cytokine_df_long) ) < nrow(cytokine_df_long)]  # Remove rows with NA only
cytokine_df_long %>% as_tibble()

# Parallelized calculation of correlations
cytokine_corrs <-
  corr_loop_parallel(cytokine_df_long, df_all, obj.name = "Cytokine__NoStimulant")

cytokine_corrs_FDR <- cytokine_corrs %>% 
  corr_feature_tissue_FDR() %>% 
  select(-object_name)

cyt_to_mets <- 'data/correlations/cytokine-to-metabolite/'
write.xlsx(
  cytokine_corrs_FDR,
  glue(
    '{cyt_to_mets}/2022-07-06_cytokine-metabolite_correlations.xlsx'
  ),
  overwrite = TRUE
)
saveRDS(cytokine_corrs_FDR, 
        file = glue('{cyt_to_mets}/2022-07-06_cytokine-metabolite_correlations.rds'))


# GI_bioplex_corrs <- list(
#   "gi_media_all_mets" = gi_media_all_mets_fdr,
# )
# # write.xlsx(GI_bioplex_corrs, 'data/correlations/cytokine_metabolite/GI_Bioplex_Correlations.xlsx')
# 
# 
# heatmap_gi_media_all_mets <- corr_heatmap_facet(gi_media_all_mets_fdr)
# ggsave(heatmap_gi_media_all_mets, filename = "figures/correlations/heatmap_gi_media_all_mets.svg",
#        width = 17, height = 6.5)
# 
# # Brain Cytokine Profiles ----
# 
# df_brain_media<-
#   read_xlsx('input_files/Behavior_Immune_data/AMA27 - brain cytokines.xlsx', sheet = 'Sheet1') %>% 
#   select(-c(Study, Study_group)) %>% 
#   right_join(meta_animal_id, by = "Animal") %>% 
#   select(-c(Animal)) %>% 
#   arrange(match(SampleID, rownames(df_all))) %>% 
#   column_to_rownames(var = "SampleID") %>% 
#   mutate_all(as.numeric)
# 
# # select cytokines with a minimum of 4 obtained values
# df_brain_media_trim <- df_brain_media[, colSums(!is.na(df_brain_media)) > 4]
# 
# brain_all_mets <- corr_loop_parallel(df_brain_media_trim, df_all, obj.name = "Cytokine__Brain__NoStimulant")
# 
# brain_all_mets_fdr <- brain_all_mets %>% 
#   group_by(feature_A) %>%
#   mutate(q = p.adjust(p, method = 'BH')) %>%
#   ungroup() %>% 
#   mutate(feature_A = gsub("_y4", "", feature_A)) %>% 
#   mutate(feature_A = toupper(feature_A)) %>% 
#   separate(feature_A, c("tissue_A", "feature_A_obj"), sep = "_", remove = F) %>% 
#   separate(feature_B, c("feature_B_obj", "tissue_B"), sep = "__", remove = F)
# 
# # brain_bioplex_corrs <- list(
# #   "Brain_tissue_all_mets" = brain_all_mets_fdr
# # )
# # write.xlsx(brain_bioplex_corrs,
# #            'data/correlations/cytokine_metabolite/Brain_tissue_Bioplex_Correlations.xlsx')
# 
# 
# heatmap_brain_media_all_mets <-
#   corr_heatmap_facet(brain_all_mets_fdr,
#                      facet_ord = c("PLASMA", "CSF", "CORTEX", "HIPP", "AC"))
# 
# ggsave(heatmap_brain_media_all_mets, filename = "figures/correlations/heatmap_brain_media_all_mets.svg",
#        width = 17, height = 6.5)
# 
# 
# # Blood Cytokine Profiles ----
# 
# plasma_dfprep_helper <- function(df){
#   df %>% 
#     select(-c(Gender, treatment)) %>% 
#     right_join(meta_animal_id, by = "Animal") %>% 
#     select(-c(Animal)) %>% 
#     arrange(match(SampleID, rownames(df_all))) %>% 
#     column_to_rownames(var = "SampleID")
# }
# 
# df_plasma_media<-
#   read_xlsx('input_files/Behavior_Immune_data/Year 4 Blood Cytokine data.xlsx', sheet = 'Media') %>% 
#   plasma_dfprep_helper()
# 
# 
# plasma_media_all_mets <- corr_loop_parallel(df_plasma_media, df_all, obj.name = "Cytokine__Plasma__NoStimulant")
# plasma_media_all_mets_fdr <- corr_heatmap_FDRv2_plasma(plasma_media_all_mets)
# 
# 
# # plasma_bioplex_corrs <- list(
# #   "plasma_Media_all_mets" = plasma_media_all_mets_fdr,
# # )
# # write.xlsx(plasma_bioplex_corrs,
# #            'data/correlations/cytokine_metabolite/Blood_Bioplex_Correlations.xlsx')
# 
# heatmap_plasma_media_all_mets <- corr_heatmap(plasma_media_all_mets_fdr)
# ggsave(heatmap_plasma_media_all_mets, filename = "figures/correlations/heatmap_plasma_media_all_mets.svg",
#        width = 8, height = 6.5)


#_______________________________________________________________________________
#                Metabolite x Behavior Correlations --------
#_______________________________________________________________________________

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
write.xlsx(behavior_corrs,
           'data/correlations/behavior_metabolite/Behavior_Metabolite_Correlations.xlsx')




