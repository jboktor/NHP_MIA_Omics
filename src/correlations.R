
source("src/_load_packages.R")
source("src/_misc_functions.R")

#' This script generates metabolite-to-metabolite correlations within and 
#' across tissues

# Create one dataframe with all scaled / imputed data
df_plasma_scaled_imp <-
  read_xlsx('input_files/Metabolomics/NHP Plasma HD4 Metabolon.xlsx', sheet = 'ScaledImpData') %>% 
  mutate(tissue = "Plasma")
df_csf_scaled_imp <-
  read_xlsx('input_files/Metabolomics/NHP CSF Metabolon.xlsx', sheet = 'ScaledImpData') %>% 
  mutate(tissue = "CSF")
df_jejunum_scaled_imp <-
  read_xlsx('input_files/Metabolomics/NHP Jejunum Metabolon.xlsx', sheet = 'ScaledImpData') %>% 
  mutate(tissue = "Jejunum")
df_ileum_scaled_imp <-
  read_xlsx('input_files/Metabolomics/NHP Ileum Metabolon.xlsx', sheet = 'ScaledImpData') %>% 
  mutate(tissue = "Ileum")
df_colon_scaled_imp <-
  read_xlsx('input_files/Metabolomics/NHP Colon Metabolon.xlsx', sheet = 'ScaledImpData') %>% 
  mutate(tissue = "Colon")
df_feces_scaled_imp <-
  read_xlsx('input_files/Metabolomics/NHP Feces Metabolon.xlsx', sheet = 'ScaledImpData') %>% 
  mutate(tissue = "Feces")
keys <- read_xlsx('input_files/misc/sample_keys.xlsx', sheet = 'keys') %>%
  select(SampleID, group, treatment)
keys2 <- read_xlsx('input_files/misc/sample_keys.xlsx', sheet = 'keys') %>%
  select(subject_id, SampleID) %>% 
  dplyr::rename("Animal" = "subject_id" )

df_scaled_imp <- bind_rows(df_plasma_scaled_imp,
                           df_csf_scaled_imp,
                           df_jejunum_scaled_imp,
                           df_ileum_scaled_imp,
                           df_colon_scaled_imp,
                           df_feces_scaled_imp)

df_all <- df_scaled_imp %>% 
  select(-c(`SUPER PATHWAY`, `SUB PATHWAY`)) %>% 
  pivot_longer(!c(BIOCHEMICAL, tissue), names_to = "SampleID") %>% 
  left_join(keys, by = "SampleID") %>%
  mutate(analyte__tissue = paste(BIOCHEMICAL, tissue, sep = "__")) %>% 
  select(SampleID, value, analyte__tissue) %>%
  pivot_wider(values_from = "value", names_from = "analyte__tissue") %>% 
  column_to_rownames(var = "SampleID")

# saveRDS(df_scaled_imp, file = "input_files/Metabolomics/scaled_&_imputed_abundance_annotated.rds")
# saveRDS(df_all, file = "input_files/Metabolomics/scaled_&_imputed_abundance.rds")

# Superpathway means ----
df_superpath <- super_path_avg(df_scaled_imp, keys = keys) %>% 
  mutate(path__tissue = paste(`SUPER PATHWAY`, tissue, sep = "__")) %>% 
  ungroup() %>% 
  select(SampleID, path_mean, path__tissue) %>%
  pivot_wider(values_from = "path_mean", names_from = "path__tissue") %>% 
  column_to_rownames(var = "SampleID")

# Subpathway means ----
df_subpath <- sub_path_avg(df_scaled_imp, keys = keys) %>% 
  mutate(path__tissue = paste(`SUB PATHWAY`, tissue, sep = "__")) %>% 
  ungroup() %>% 
  select(SampleID, path_mean, path__tissue) %>%
  pivot_wider(values_from = "path_mean", names_from = "path__tissue") %>% 
  column_to_rownames(var = "SampleID")


#_______________________________________________________________________________
#                Metabolite x Metabolite Correlations --------
#_______________________________________________________________________________

# Correlation calculations
cor_superpath <- rcorr(as.matrix(df_superpath), type = "spearman")
df_cor_superpath <- tidy_cor_output(cor_superpath) %>% 
  separate(feature_A, c("feature_A", "tissue_A"), sep = "__") %>% 
  separate(feature_B, c("feature_B", "tissue_B"), sep = "__")

cor_subpath <- rcorr(as.matrix(df_subpath), type = "spearman")
df_cor_subpath <- tidy_cor_output(cor_superpath) %>% 
  separate(feature_A, c("feature_A", "tissue_A"), sep = "__") %>% 
  separate(feature_B, c("feature_B", "tissue_B"), sep = "__")

# cor_all <- rcorr(as.matrix(df_all), type = "spearman")
df_cor_all <- rcorr(as.matrix(df_all), type = "spearman") %>% 
  tidy_cor_output() %>% 
  separate(feature_A, c("feature_A", "tissue_A"), sep = "__") %>% 
  separate(feature_B, c("feature_B", "tissue_B"), sep = "__")

pathway_map <- df_scaled_imp %>% 
  select(BIOCHEMICAL, `SUPER PATHWAY`, `SUB PATHWAY`) %>% 
  distinct()

# write.xlsx(pathway_map, 'input_files/misc/pathway_map.xlsx')
# write.xlsx(df_cor_superpath, 'data/correlations/metabolite-to-metabolite/superpathway_corrs.xlsx')
# write.xlsx(df_cor_subpath, 'data/correlations/metabolite-to-metabolite/subpathway_corrs.xlsx')
# write.xlsx(df_cor_all, 'data/correlations/metabolite-to-metabolite/all_metabolite_corrs.xlsx')

#_______________________________________________________________________________
#                Metabolite x Cytokine Correlations --------
#_______________________________________________________________________________

# GI Cytokine Profiles ----
df_GI_ConA <-
  read_xlsx('input_files/Behavior_Immune_data/AMA 27 GI tissues ConA 072513.xlsx', sheet = 'Sheet1') %>% 
  select(-c(Gender, Study_group, Description)) %>% 
  pivot_longer(!c(Animal, `GI Section`), names_to = "cytokine") %>% 
  separate(cytokine, c("cytokine", "trash"), sep = " ") %>% 
  mutate(cytokine__tissue = paste(cytokine, `GI Section`, sep = "__")) %>% 
  right_join(keys2, by = "Animal") %>% 
  select(SampleID, value, cytokine__tissue) %>% 
  pivot_wider(names_from = "cytokine__tissue", values_from = "value") %>% 
  arrange(match(SampleID, rownames(df_all))) %>% 
  column_to_rownames(var = "SampleID") %>% 
  select(-"NA")

df_GI_PolyIC <-
  read_xlsx('input_files/Behavior_Immune_data/AMA 27 GI tissues Poly IC 073013.xlsx', sheet = 'Sheet1') %>% 
  select(-c(Gender, Study_group, Description)) %>% 
  pivot_longer(!c(Animal, `GI Section`), names_to = "cytokine") %>% 
  separate(cytokine, c("cytokine", "trash"), sep = " ") %>% 
  mutate(cytokine__tissue = paste(cytokine, `GI Section`, sep = "__")) %>% 
  right_join(keys2, by = "Animal") %>% 
  select(SampleID, value, cytokine__tissue) %>% 
  pivot_wider(names_from = "cytokine__tissue", values_from = "value") %>% 
  arrange(match(SampleID, rownames(df_all))) %>% 
  column_to_rownames(var = "SampleID") %>% 
  select(-"NA")

df_GI_PMAiono<-
  read_xlsx('input_files/Behavior_Immune_data/AMA 27 GI tissues PMA iono 073013.xlsx', sheet = 'Sheet1') %>% 
  select(-c(Gender, Study_group, Description)) %>% 
  pivot_longer(!c(Animal, `GI Section`), names_to = "cytokine") %>% 
  separate(cytokine, c("cytokine", "trash"), sep = " ") %>% 
  mutate(cytokine__tissue = paste(cytokine, `GI Section`, sep = "__")) %>% 
  right_join(keys2, by = "Animal") %>% 
  select(SampleID, value, cytokine__tissue) %>% 
  pivot_wider(names_from = "cytokine__tissue", values_from = "value") %>% 
  arrange(match(SampleID, rownames(df_all))) %>% 
  column_to_rownames(var = "SampleID") %>% 
  select(-"NA")

df_GI_media<-
  read_xlsx('input_files/Behavior_Immune_data/AMA 27 GI tissues media 072513.xlsx', sheet = 'Sheet1') %>% 
  select(-c(Gender, Study_group, Description)) %>% 
  pivot_longer(!c(Animal, `GI Section`), names_to = "cytokine") %>% 
  separate(cytokine, c("cytokine", "trash"), sep = " ") %>% 
  mutate(cytokine__tissue = paste(cytokine, `GI Section`, sep = "__")) %>% 
  right_join(keys2, by = "Animal") %>% 
  select(SampleID, value, cytokine__tissue) %>% 
  pivot_wider(names_from = "cytokine__tissue", values_from = "value") %>% 
  arrange(match(SampleID, rownames(df_all))) %>% 
  column_to_rownames(var = "SampleID") %>% 
  select(-"NA")

GI_cytokines <- list("Media"=df_GI_media,
                     "PMAiono" = df_GI_PMAiono,
                     "PolyIC" = df_GI_PolyIC,
                     "ConA" = df_GI_ConA)
saveRDS(GI_cytokines, file = "input_files/Behavior_Immune_data/GI_cytokines_all.rds")
saveRDS(df_GI_ConA, file = "input_files/Behavior_Immune_data/GI_ConA.rds")
saveRDS(df_GI_PolyIC, file = "input_files/Behavior_Immune_data/GI_PolyIC.rds")
saveRDS(df_GI_PMAiono, file = "input_files/Behavior_Immune_data/GI_PMAiono.rds")
saveRDS(df_GI_media, file = "input_files/Behavior_Immune_data/GI_media.rds")



gi_cona_sup <- corr_loop_parallel(df_GI_ConA, df_superpath, obj.name = "Cytokine__ConA")
gi_cona_sub <- corr_loop_parallel(df_GI_ConA, df_subpath, obj.name = "Cytokine__ConA")

gi_polyic_sup <- corr_loop_parallel(df_GI_PolyIC, df_superpath, obj.name = "Cytokine_Poly(I:C)")
gi_polyic_sub <- corr_loop_parallel(df_GI_PolyIC, df_subpath, obj.name = "Cytokine_Poly(I:C)")

gi_pmaiono_sup <- corr_loop_parallel(df_GI_PMAiono, df_superpath, obj.name = "Cytokine__PMA_ionomyocin")
gi_pmaiono_sub <- corr_loop_parallel(df_GI_PMAiono, df_subpath, obj.name = "Cytokine__PMA_ionomyocin")

gi_media_allMets <- corr_loop_parallel(df_GI_media, df_all, obj.name = "Cytokine__NoStimulant")
gi_media_sup <- corr_loop_parallel(df_GI_media, df_superpath, obj.name = "Cytokine__NoStimulant")
gi_media_sub <- corr_loop_parallel(df_GI_media, df_subpath, obj.name = "Cytokine__NoStimulant")


gi_cona_sup_fdr <- corr_heatmap_FDR(gi_cona_sup)
gi_cona_sub_fdr <- corr_heatmap_FDR(gi_cona_sub)
gi_polyic_sup_fdr <- corr_heatmap_FDR(gi_polyic_sup)
gi_polyic_sub_fdr <- corr_heatmap_FDR(gi_polyic_sub)
gi_pmaiono_sup_fdr <- corr_heatmap_FDR(gi_pmaiono_sup)
gi_pmaiono_sub_fdr <- corr_heatmap_FDR(gi_pmaiono_sub)
gi_media_sup_fdr <- corr_heatmap_FDR(gi_media_sup)
gi_media_sub_fdr <- corr_heatmap_FDR(gi_media_sub)
gi_media_allMets_fdr <- corr_heatmap_FDR(gi_media_allMets)

GI_bioplex_corrs <- list(
  "GI_Media_allMets" = gi_media_allMets_fdr,
  "GI_Media_Superpathways" = gi_media_sup_fdr,
  "GI_Media_Subpathways" = gi_media_sub_fdr,
  "GI_ConA_Superpathways" = gi_cona_sup_fdr,
  "GI_ConA_Subpathways" = gi_cona_sub_fdr,
  "GI_PolyIC_Superpathways" = gi_polyic_sup_fdr,
  "GI_PolyIC_Subpathways" = gi_polyic_sub_fdr,
  "GI_PMAionomycin_Superpathways" = gi_pmaiono_sup_fdr,
  "GI_PMAionomycin_Subpathways" = gi_pmaiono_sub_fdr
)


# write.xlsx(GI_bioplex_corrs, 'data/correlations/cytokine_metabolite/GI_Bioplex_Correlations.xlsx')

heatmap_gi_cona_superpath <- corr_heatmap_facet(gi_cona_sup_fdr)
ggsave(heatmap_gi_cona_superpath, filename = "figures/correlations/heatmap_GI_ConA_superpaths.svg",
       width = 17, height = 6.5)
heatmap_gi_cona_subpath <- corr_heatmap_facet(gi_cona_sub_fdr)
ggsave(heatmap_gi_cona_subpath, filename = "figures/correlations/heatmap_GI_ConA_subpaths.svg",
       width = 17, height = 6.5)

heatmap_gi_polyic_superpath <- corr_heatmap_facet(gi_polyic_sup_fdr)
ggsave(heatmap_gi_polyic_superpath, filename = "figures/correlations/heatmap_GI_PolyIC_superpaths.svg",
       width = 17, height = 6.5)
heatmap_gi_polyic_subpath <- corr_heatmap_facet(gi_polyic_sub_fdr)
ggsave(heatmap_gi_polyic_subpath, filename = "figures/correlations/heatmap_GI_PolyIC_subpaths.svg",
       width = 17, height = 6.5)

heatmap_gi_pmaiono_superpath <- corr_heatmap_facet(gi_pmaiono_sup_fdr)
ggsave(heatmap_gi_pmaiono_superpath, filename = "figures/correlations/heatmap_GI_PMAiono_superpaths.svg",
       width = 17, height = 6.5)
heatmap_gi_pmaiono_subpath <- corr_heatmap_facet(gi_pmaiono_sub_fdr)
ggsave(heatmap_gi_pmaiono_subpath, filename = "figures/correlations/heatmap_GI_PMAiono_subpaths.svg",
       width = 17, height = 6.5)

heatmap_gi_media_superpath <- corr_heatmap_facet(gi_media_sup_fdr)
ggsave(heatmap_gi_media_superpath, filename = "figures/correlations/heatmap_GI_media_superpaths.svg",
       width = 17, height = 6.5)
heatmap_gi_media_subpath <- corr_heatmap_facet(gi_media_sub_fdr)
ggsave(heatmap_gi_media_subpath, filename = "figures/correlations/heatmap_GI_media_subpaths.svg",
       width = 17, height = 6.5)
heatmap_gi_media_allMets <- corr_heatmap_facet(gi_media_allMets_fdr)
ggsave(heatmap_gi_media_allMets, filename = "figures/correlations/heatmap_GI_media_allMets.svg",
       width = 17, height = 6.5)

# Brain Cytokine Profiles ----
df_brain_media<-
  read_xlsx('input_files/Behavior_Immune_data/AMA27 - brain cytokines.xlsx', sheet = 'Sheet1') %>% 
  select(-c(Study, Study_group)) %>% 
  right_join(keys2, by = "Animal") %>% 
  select(-c(Animal)) %>% 
  arrange(match(SampleID, rownames(df_all))) %>% 
  column_to_rownames(var = "SampleID") %>% 
  mutate_all(as.numeric)


# select cytokines with a minimum of 4 obtained values
df_brain_media_trim <- df_brain_media[, colSums(!is.na(df_brain_media)) > 4]

brain_sup <- corr_loop_parallel(df_brain_media_trim, df_superpath, obj.name = "Cytokine__Brain__NoStimulant")
brain_sub <- corr_loop_parallel(df_brain_media_trim, df_subpath, obj.name = "Cytokine__Brain__NoStimulant")
brain_allMets <- corr_loop_parallel(df_brain_media_trim, df_all, obj.name = "Cytokine__Brain__NoStimulant")

brain_sup_fdr <- brain_sup %>% 
  group_by(feature_A) %>%
  mutate(q = p.adjust(p, method = 'BH')) %>%
  ungroup() %>% 
  mutate(feature_A = gsub("_y4", "", feature_A)) %>% 
  mutate(feature_A = toupper(feature_A)) %>% 
  separate(feature_A, c("tissue_A", "feature_A_obj"), sep = "_", remove = F) %>% 
  separate(feature_B, c("feature_B_obj", "tissue_B"), sep = "__", remove = F)

brain_sub_fdr <- brain_sub %>% 
  group_by(feature_A) %>%
  mutate(q = p.adjust(p, method = 'BH')) %>%
  ungroup() %>% 
  mutate(feature_A = gsub("_y4", "", feature_A)) %>% 
  mutate(feature_A = toupper(feature_A)) %>% 
  separate(feature_A, c("tissue_A", "feature_A_obj"), sep = "_", remove = F) %>% 
  separate(feature_B, c("feature_B_obj", "tissue_B"), sep = "__", remove = F)

brain_allMets_fdr <- brain_allMets %>% 
  group_by(feature_A) %>%
  mutate(q = p.adjust(p, method = 'BH')) %>%
  ungroup() %>% 
  mutate(feature_A = gsub("_y4", "", feature_A)) %>% 
  mutate(feature_A = toupper(feature_A)) %>% 
  separate(feature_A, c("tissue_A", "feature_A_obj"), sep = "_", remove = F) %>% 
  separate(feature_B, c("feature_B_obj", "tissue_B"), sep = "__", remove = F)

brain_bioplex_corrs <- list(
  "Brain_tissue_allMets" = brain_allMets_fdr,
  "Brain_tissue_Superpathways" = brain_sup_fdr,
  "Brain_tissue_Subpathways" = brain_sub_fdr
)

# write.xlsx(brain_bioplex_corrs,
#            'data/correlations/cytokine_metabolite/Brain_tissue_Bioplex_Correlations.xlsx')


heatmap_brain_media_superpath <- 
  corr_heatmap_facet(brain_sup_fdr, 
                     facet_ord = c("PLASMA", "CSF", "CORTEX", "HIPP", "AC"))
ggsave(heatmap_brain_media_superpath, filename = "figures/correlations/heatmap_brain_media_superpaths.svg",
       width = 17, height = 6.5)

heatmap_brain_media_subpath <-
  corr_heatmap_facet(brain_sub_fdr,
                     facet_ord = c("PLASMA", "CSF", "CORTEX", "HIPP", "AC"))
ggsave(heatmap_brain_media_subpath, filename = "figures/correlations/heatmap_brain_media_subpaths.svg",
       width = 17, height = 6.5)

heatmap_brain_media_allMets <-
  corr_heatmap_facet(brain_allMets_fdr,
                     facet_ord = c("PLASMA", "CSF", "CORTEX", "HIPP", "AC"))
ggsave(heatmap_brain_media_allMets, filename = "figures/correlations/heatmap_brain_media_allMets.svg",
       width = 17, height = 6.5)


# Blood Cytokine Profiles ----

plasma_dfprep_helper <- function(df){
  df %>% 
    select(-c(Gender, treatment)) %>% 
    right_join(keys2, by = "Animal") %>% 
    select(-c(Animal)) %>% 
    arrange(match(SampleID, rownames(df_all))) %>% 
    column_to_rownames(var = "SampleID")
}

df_plasma_media<-
  read_xlsx('input_files/Behavior_Immune_data/Year 4 Blood Cytokine data.xlsx', sheet = 'Media') %>% 
  plasma_dfprep_helper()
df_plasma_lps <-
  read_xlsx('input_files/Behavior_Immune_data/Year 4 Blood Cytokine data.xlsx', sheet = 'LPS') %>% 
  plasma_dfprep_helper()
df_plasma_cona <-
  read_xlsx('input_files/Behavior_Immune_data/Year 4 Blood Cytokine data.xlsx', sheet = 'ConA') %>% 
  plasma_dfprep_helper()
df_plasma_polyic <-
  read_xlsx('input_files/Behavior_Immune_data/Year 4 Blood Cytokine data.xlsx', sheet = 'Poly IC') %>% 
  plasma_dfprep_helper()
df_plasma_pha <-
  read_xlsx('input_files/Behavior_Immune_data/Year 4 Blood Cytokine data.xlsx', sheet = 'PHA') %>% 
  plasma_dfprep_helper()
df_plasma_pma24 <-
  read_xlsx('input_files/Behavior_Immune_data/Year 4 Blood Cytokine data.xlsx', sheet = 'PMA-iono 24') %>% 
  plasma_dfprep_helper()
df_plasma_pma48 <-
  read_xlsx('input_files/Behavior_Immune_data/Year 4 Blood Cytokine data.xlsx', sheet = 'PMA-iono 48') %>% 
  plasma_dfprep_helper()

cytokine_profiles_noStim <- list(
  "GI_cytokines" = df_GI_media,
  "Brain_cytokines" = df_brain_media,
  "Plasma_cytokines" = df_plasma_media
)
saveRDS(cytokine_profiles_noStim, file = "input_files/Behavior_Immune_data/cytokine_profiles_all_tissues_no_stimulant.rds")

plasma_media_sup <- corr_loop_parallel(df_plasma_media, df_superpath, obj.name = "Cytokine__Plasma__NoStimulant")
plasma_media_sub <- corr_loop_parallel(df_plasma_media, df_subpath, obj.name = "Cytokine__Plasma__NoStimulant")
plasma_media_allMets <- corr_loop_parallel(df_plasma_media, df_all, obj.name = "Cytokine__Plasma__NoStimulant")

plasma_lps_sup <- corr_loop_parallel(df_plasma_cona, df_superpath, obj.name = "Cytokine__Plasma__LPS")
plasma_lps_sub <- corr_loop_parallel(df_plasma_cona, df_subpath, obj.name = "Cytokine__Plasma__LPS")

plasma_cona_sup <- corr_loop_parallel(df_plasma_cona, df_superpath, obj.name = "Cytokine__Plasma__ConA")
plasma_cona_sub <- corr_loop_parallel(df_plasma_cona, df_subpath, obj.name = "Cytokine__Plasma__ConA")
                             
plasma_polyic_sup <- corr_loop_parallel(df_plasma_polyic, df_superpath, obj.name = "Cytokine__Plasma__PolyIC")
plasma_polyic_sub <- corr_loop_parallel(df_plasma_polyic, df_subpath, obj.name = "Cytokine__Plasma__PolyIC")

plasma_pha_sup <- corr_loop_parallel(df_plasma_pha, df_superpath, obj.name = "Cytokine__Plasma__PHA")
plasma_pha_sub <- corr_loop_parallel(df_plasma_pha, df_subpath, obj.name = "Cytokine__Plasma__PHA")

plasma_pma24_sup <- corr_loop_parallel(df_plasma_pma24, df_superpath, obj.name = "Cytokine__Plasma__PMA_iono_24")
plasma_pma24_sub <- corr_loop_parallel(df_plasma_pma24, df_subpath, obj.name = "Cytokine__Plasma__PMA_iono_24")

plasma_pma48_sup <- corr_loop_parallel(df_plasma_pma48, df_superpath, obj.name = "Cytokine__Plasma__PMA_iono_48")
plasma_pma48_sub <- corr_loop_parallel(df_plasma_pma48, df_subpath, obj.name = "Cytokine__Plasma__PMA_iono_48")


plasma_media_sup_fdr <- corr_heatmap_FDRv2_plasma(plasma_media_sup)
plasma_media_sub_fdr <- corr_heatmap_FDRv2_plasma(plasma_media_sub)
plasma_media_allMets_fdr <- corr_heatmap_FDRv2_plasma(plasma_media_allMets)
plasma_lps_sup_fdr <- corr_heatmap_FDRv2_plasma(plasma_lps_sup)
plasma_lps_sub_fdr <- corr_heatmap_FDRv2_plasma(plasma_lps_sub)
plasma_cona_sup_fdr <- corr_heatmap_FDRv2_plasma(plasma_cona_sup)
plasma_cona_sub_fdr <- corr_heatmap_FDRv2_plasma(plasma_cona_sub)
plasma_polyic_sup_fdr <- corr_heatmap_FDRv2_plasma(plasma_polyic_sup)
plasma_polyic_sub_fdr <- corr_heatmap_FDRv2_plasma(plasma_polyic_sub)
plasma_pha_sup_fdr <- corr_heatmap_FDRv2_plasma(plasma_pha_sup)
plasma_pha_sub_fdr <- corr_heatmap_FDRv2_plasma(plasma_pha_sub)
plasma_pma24_sup_fdr <- corr_heatmap_FDRv2_plasma(plasma_pma24_sup)
plasma_pma24_sub_fdr <- corr_heatmap_FDRv2_plasma(plasma_pma24_sub)
plasma_pma48_sup_fdr <- corr_heatmap_FDRv2_plasma(plasma_pma48_sup)
plasma_pma48_sub_fdr <- corr_heatmap_FDRv2_plasma(plasma_pma48_sub)


plasma_bioplex_corrs <- list(
  "plasma_Media_allMets" = plasma_media_allMets_fdr,
  "plasma_Media_Superpathways" = plasma_media_sup_fdr,
  "plasma_Media_Subpathways" = plasma_media_sub_fdr,
  "plasma_LPS_Superpathways" = plasma_lps_sup_fdr,
  "plasma_LPS_Subpathways" = plasma_lps_sub_fdr,
  "plasma_ConA_Superpathways" = plasma_cona_sup_fdr,
  "plasma_ConA_Subpathways" = plasma_cona_sub_fdr,
  "plasma_PolyIC_Superpathways" = plasma_polyic_sup_fdr,
  "plasma_PolyIC_Subpathways" = plasma_polyic_sub_fdr,
  "plasma_PHA_Superpathways" = plasma_pha_sup_fdr,
  "plasma_PHA_Subpathways" = plasma_pha_sub_fdr,
  "plasma_PMA24hr_Superpathways" = plasma_pma24_sup_fdr,
  "plasma_PMA24hr_Subpathways" = plasma_pma24_sub_fdr,
  "plasma_PMA48hr_Superpathways" = plasma_pma48_sup_fdr,
  "plasma_PMA48hr_Subpathways" = plasma_pma48_sub_fdr
)

# write.xlsx(plasma_bioplex_corrs,
#            'data/correlations/cytokine_metabolite/Blood_Bioplex_Correlations.xlsx')

heatmap_plasma_media_allMets <- corr_heatmap(plasma_media_allMets_fdr)
ggsave(heatmap_plasma_media_allMets, filename = "figures/correlations/heatmap_plasma_media_allMets.svg",
       width = 8, height = 6.5)
heatmap_plasma_media_superpath <- corr_heatmap(plasma_media_sup_fdr)
ggsave(heatmap_plasma_media_superpath, filename = "figures/correlations/heatmap_plasma_media_superpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_media_subpath <- corr_heatmap(plasma_media_sub_fdr)
ggsave(heatmap_plasma_media_subpath, filename = "figures/correlations/heatmap_plasma_media_subpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_lps_superpath <- corr_heatmap(plasma_lps_sup_fdr)
ggsave(heatmap_plasma_lps_superpath, filename = "figures/correlations/heatmap_plasma_LPS_superpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_lps_subpath <- corr_heatmap(plasma_lps_sub_fdr)
ggsave(heatmap_plasma_lps_subpath, filename = "figures/correlations/heatmap_plasma_LPS_subpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_cona_superpath <- corr_heatmap(plasma_cona_sup_fdr)
ggsave(heatmap_plasma_cona_superpath, filename = "figures/correlations/heatmap_plasma_ConA_superpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_cona_subpath <- corr_heatmap(plasma_cona_sub_fdr)
ggsave(heatmap_plasma_cona_subpath, filename = "figures/correlations/heatmap_plasma_ConA_subpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_polyic_superpath <- corr_heatmap(plasma_polyic_sup_fdr)
ggsave(heatmap_plasma_polyic_superpath, filename = "figures/correlations/heatmap_plasma_PolyIC_superpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_polyic_subpath <- corr_heatmap(plasma_polyic_sub_fdr)
ggsave(heatmap_plasma_polyic_subpath, filename = "figures/correlations/heatmap_plasma_PolyIC_subpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_pha_superpath <- corr_heatmap(plasma_pha_sup_fdr)
ggsave(heatmap_plasma_pha_superpath, filename = "figures/correlations/heatmap_plasma_PHA_superpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_pha_subpath <- corr_heatmap(plasma_pha_sub_fdr)
ggsave(heatmap_plasma_pha_subpath, filename = "figures/correlations/heatmap_plasma_PHA_subpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_pma24_superpath <- corr_heatmap(plasma_pma24_sup_fdr)
ggsave(heatmap_plasma_pma24_superpath, filename = "figures/correlations/heatmap_plasma_pma24_superpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_pma24_subpath <- corr_heatmap(plasma_pma24_sub_fdr)
ggsave(heatmap_plasma_pma24_subpath, filename = "figures/correlations/heatmap_plasma_pma24_subpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_pma48_superpath <- corr_heatmap(plasma_pma48_sup_fdr)
ggsave(heatmap_plasma_pma48_superpath, filename = "figures/correlations/heatmap_plasma_pma48_superpaths.svg",
       width = 8, height = 6.5)
heatmap_plasma_pma48_subpath <- corr_heatmap(plasma_pma48_sub_fdr)
ggsave(heatmap_plasma_pma48_subpath, filename = "figures/correlations/heatmap_plasma_pma48_subpaths.svg",
       width = 8, height = 6.5)


# Behavioral data ----
df_behavior_postwean <-
  read_xlsx(
    'input_files/Behavior_Immune_data/LongStereo - Year four stereotyped behaviors.xlsx',
    sheet = 'PostWeanStereowithCodes') %>%
  select(-c(COHORT, FOCALDYE,GENDER, EXPCODE, EXPCODE2)) %>% 
  pivot_longer(!c(Animal, Observation), names_to = "behavior") %>% 
  mutate(behavior__timepoint = paste(behavior, Observation, sep = "__")) %>% 
  right_join(keys2, by = "Animal") %>% 
  right_join(keys, by = "SampleID")
df_behavior_juvenile <-
  read_xlsx(
    'input_files/Behavior_Immune_data/LongStereo - Year four stereotyped behaviors.xlsx',
    sheet = 'JuvStereowithCodes') %>%
  select(-c(COHORT, FOCALDYE,GENDER, EXPCODE, EXPCODE2)) %>% 
  pivot_longer(!c(Animal, Observation), names_to = "behavior") %>% 
  mutate(behavior__timepoint = paste(behavior, Observation, sep = "__")) %>% 
  right_join(keys2, by = "Animal") %>% 
  right_join(keys, by = "SampleID")
df_behavior_adult <-
  read_xlsx(
    'input_files/Behavior_Immune_data/LongStereo - Year four stereotyped behaviors.xlsx',
    sheet = 'SubAdultwithCodes') %>%
  select(-c(GENDER, EXPCODE, EXPCODE2)) %>% 
  pivot_longer(!c(Animal, Observation), names_to = "behavior") %>% 
  mutate(behavior__timepoint = paste(behavior, Observation, sep = "__")) %>% 
  right_join(keys2, by = "Animal") %>% 
  right_join(keys, by = "SampleID")

# Data frame for exploration of behavior data
df_behav_eda <-
  bind_rows(df_behavior_postwean,
            df_behavior_juvenile,
            df_behavior_adult)

df_behav_eda %>% 
  filter(behavior == "all_stereos") %>% 
  ggplot(aes(x=group, y= value )) +
  geom_boxplot() +
  geom_point()

df_behavior <- df_behav_eda %>% 
  select(SampleID, value, behavior__timepoint) %>%
  pivot_wider(names_from = "behavior__timepoint", values_from = "value") %>%
  arrange(match(SampleID, rownames(df_all))) %>%
  column_to_rownames(var = "SampleID")

behavior_data <- list("df_behavior" = df_behavior, 
                      "df_behav_eda" = df_behav_eda)

saveRDS(behavior_data, file = "input_files/Behavior_Immune_data/behaviors.rds")

behavior_sup <- corr_loop_parallel(df_behavior, df_superpath, obj.name = "Behavior")
behavior_sub <- corr_loop_parallel(df_behavior, df_subpath, obj.name = "Behavior")
behavior_all_mets <- corr_loop_parallel(df_behavior, df_all, obj.name = "Behavior")

behavior_sup_fdr <- corr_heatmap_FDRv3_behavior(behavior_sup)
behavior_sub_fdr <- corr_heatmap_FDRv3_behavior(behavior_sub)
behavior_all_mets_fdr <- corr_heatmap_FDRv3_behavior(behavior_all_mets)

heatmap_behavior_superath <- corr_heatmap_facet_behav(behavior_sup_fdr)
heatmap_behavior_subpath <- corr_heatmap_facet_behav(behavior_sub_fdr)
heatmap_behavior_allMets <- corr_heatmap_facet_behav(behavior_all_mets_fdr)

ggsave(
  heatmap_behavior_superath,
  filename = "figures/correlations/heatmap_behavior_superpaths.svg",
  width = 7,
  height = 6.5
)
ggsave( 
  heatmap_behavior_subpath,
  filename = "figures/correlations/heatmap_behavior_subpaths.svg",
  width = 7,
  height = 6.5
)
ggsave(
  heatmap_behavior_allMets,
  filename = "figures/correlations/heatmap_behavior_allMets.svg",
  width = 8,
  height = 6.5
)

behavior_corrs <- list(
  "behaviors__Superpathways" = behavior_sup_fdr,
  "behaviors__Subpathways" = behavior_sub_fdr,
  "behaviors_allMetabolites" = behavior_all_mets_fdr)

# write.xlsx(behavior_corrs,
#            'data/correlations/behavior_metabolite/Behavior_Metabolite_Correlations.xlsx')



#_______________________________________________________________________________
# Saving a clean cytokine correlation datatable for supplement ----


# GI Cytokine Profiles 
df_GI_media <-
  read_xlsx('input_files/Behavior_Immune_data/AMA 27 GI tissues media 072513.xlsx', sheet = 'Sheet1') %>% 
  select(-c(Gender, Study_group, Description)) %>% 
  pivot_longer(!c(Animal, `GI Section`), names_to = "cytokine") %>% 
  separate(cytokine, c("cytokine", "trash"), sep = " ") %>% 
  right_join(keys2, by = "Animal") %>% 
  select(-trash) %>% 
  pivot_wider(names_from = "cytokine", values_from = "value") %>%
  select(-"NA") %>% 
  dplyr::rename(Tissue = `GI Section`) %>% 
  drop_na(Tissue)

# Brain Cytokine Profiles 
df_brain_media <-
  read_xlsx('input_files/Behavior_Immune_data/AMA27 - brain cytokines.xlsx', sheet = 'Sheet1') %>% 
  select(-c(Study, Study_group)) %>% 
  pivot_longer(!Animal, names_to = "cytokine_longname") %>% 
  mutate(cytokine_longname = gsub("_y4", "", cytokine_longname)) %>% 
  separate(cytokine_longname, c("Tissue", "cytokine"), sep = "_") %>% 
  mutate(cytokine = toupper(cytokine)) %>% 
  mutate(cytokine = str_replace(cytokine, "IL-", "IL")) %>% 
  mutate(cytokine = str_replace(cytokine, "GCSF", "G-CSF")) %>% 
  mutate(cytokine = str_replace(cytokine, "IL12P40", "IL12p40")) %>% 
  mutate(cytokine = str_replace(cytokine, "MCP1", "MCP-1")) %>% 
  mutate(cytokine = str_replace(cytokine, "MIP-1A", "MIP1A")) %>% 
  mutate(cytokine = str_replace(cytokine, "MCP1", "MCP-1")) %>% 
  mutate(cytokine = str_replace(cytokine, "SCD40L", "sCD40L")) %>% 
  pivot_wider(names_from = "cytokine", values_from = "value") %>%
  right_join(keys2, by = "Animal") %>% 
  relocate(SampleID, .after = Tissue) %>% 
  arrange(Tissue)

# Plasma Profiles 
df_plasma_media <-
  read_xlsx('input_files/Behavior_Immune_data/Year 4 Blood Cytokine data.xlsx', sheet = 'Media') %>% 
  select(-c(Gender, treatment)) %>% 
  pivot_longer(!Animal, names_to = "cytokine") %>% 
  mutate(cytokine = gsub("m_y4_", "", cytokine)) %>% 
  mutate(cytokine = toupper(cytokine)) %>% 
  mutate(cytokine = str_replace(cytokine, "IL-", "IL")) %>% 
  mutate(cytokine = str_replace(cytokine, "GCSF", "G-CSF")) %>% 
  mutate(cytokine = str_replace(cytokine, "IL12P40", "IL12p40")) %>% 
  mutate(cytokine = str_replace(cytokine, "MCP1", "MCP-1")) %>% 
  mutate(cytokine = str_replace(cytokine, "MIP-1A", "MIP1A")) %>% 
  mutate(cytokine = str_replace(cytokine, "MCP1", "MCP-1")) %>% 
  mutate(cytokine = str_replace(cytokine, "SCD40L", "sCD40L")) %>% 
  pivot_wider(names_from = "cytokine", values_from = "value") %>%
  right_join(keys2, by = "Animal") %>% 
  mutate(Tissue = "Blood") %>% 
  relocate(Tissue, .after = Animal) %>% 
  relocate(SampleID, .after = Tissue) 

  
final_cytokine_profile_inputs <- list(
  "GI Tissue" = df_GI_media,
  "Brain Tissue" = df_brain_media,
  "Blood Tissue" = df_plasma_media
)

write.xlsx(final_cytokine_profile_inputs,
           'data/correlations/cytokine-to-metabolite/Cytokine_profiles.xlsx', 
           overwrite = T)


