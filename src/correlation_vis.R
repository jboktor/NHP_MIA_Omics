
source("src/_load_packages.R")
source("src/_misc_functions.R")

# Load data 
pathway_map <-
  read_excel("input_files/misc/pathway_map.xlsx", sheet = "Sheet 1")
pathway_map_join <-
  pathway_map %>% dplyr::rename("feature_A" = "BIOCHEMICAL")
pathway_map_sub2sup <- pathway_map %>%
  janitor::clean_names() %>%
  select(super_pathway, sub_pathway) %>%
  distinct()
#_______________________________________________________________________________
#                           Loading Correlation data                        ---- 
#_______________________________________________________________________________

# gi_media_allMets_fdr <-
#   read_excel("data/correlations_cytokine_metabolite/GI_Bioplex_Correlations.xlsx",
#              sheet = "GI_Media_allMets") %>% 
#   mutate(tissue_A_group = "GI")
# 
# brain_allMets_fdr <-
#   read_excel("data/correlations_cytokine_metabolite/Brain_tissue_Bioplex_Correlations.xlsx",
#              sheet = "Brain_tissue_allMets") %>% 
#   mutate(tissue_A = if_else(grepl("PLASMA", tissue_A), paste0("Brain Plasma"), tissue_A)) %>% 
#   mutate(tissue_A_group = "Brain")
# 
# plasma_media_allMets_fdr <-
#   read_excel("data/correlations_cytokine_metabolite/Blood_Bioplex_Correlations.xlsx",
#              sheet = "plasma_Media_allMets") %>% 
#   mutate(tissue_A_group = "Plasma")

# cytokine_cors <- bind_rows(gi_media_allMets_fdr, brain_allMets_fdr, plasma_media_allMets_fdr)

cyt_to_mets <- 'data/correlations/cytokine-to-metabolite/'

cytokine_cors <- readRDS(glue(
  '{cyt_to_mets}/2022-07-06_cytokine-metabolite_correlations.rds'
))

behavior_all_mets_fdr <-
  read_excel("data/correlations/behavior-to-metabolite/Behavior_Metabolite_Correlations.xlsx",
             sheet = "behaviors_allMetabolites")

cyt_tissues <- 
  c("Plasma-PBMCs", "Plasma",  "CSF", "Cortex", "HIPP", "AC",
  "Jejunum-PBMCs", "Ileum-PBMCs", "MLN-PBMCs", "Colon-PBMCs")

plasma_cyt_cors <- 
  cytokine_cors %>% 
  filter(tissue_B == "Plasma")

csf_cyt_cors <- 
  cytokine_cors %>% 
  filter(tissue_B == "CSF")

GI_cyt_cors <- 
  cytokine_cors %>% 
  filter(tissue_B %in%  c("Jejunum", "Ileum", "Colon"))

feces_cyt_cors <- 
  cytokine_cors %>% 
  filter(tissue_B == "Feces")

#_______________________________________________________________________________
#                    Tissue source specific Heatmaps  ----
#_______________________________________________________________________________
 
plasma_cyt_cors.plot <- corr_heatmap_facet_v2(plasma_cyt_cors, facet_ord =cyt_tissues)
ggsave(plasma_cyt_cors.plot, filename = "figures/correlations/2022-07-05_heatmap_Plasma_allmets_allcytokines.svg",
       width = 25, height =5)
plasma_cyt_cors.plot <- corr_heatmap_facet_v2(plasma_cyt_cors, facet_ord =cyt_tissues)
ggsave(plasma_cyt_cors.plot, filename = "figures/correlations/2022-07-05_heatmap_Plasma_allmets_allcytokines.svg",
       width = 25, height =5)
csf_cyt_cors.plot <- corr_heatmap_facet_v2(csf_cyt_cors, facet_ord =cyt_tissues)
ggsave(csf_cyt_cors.plot, filename = "figures/correlations/2022-07-05_heatmap_CSF_allmets_allcytokines.svg",
       width = 25, height =4)
GI_cyt_cors.plot <- corr_heatmap_facet_v2(GI_cyt_cors, facet_ord =cyt_tissues)
ggsave(GI_cyt_cors.plot, filename = "figures/correlations/2022-07-05_heatmap_GI_allmets_allcytokines.svg",
       width = 25, height =6)
feces_cyt_cors.plot <- corr_heatmap_facet_v2(feces_cyt_cors, facet_ord =cyt_tissues)
ggsave(feces_cyt_cors.plot, filename = "figures/correlations/2022-07-05_heatmap_Feces_allmets_allcytokines.svg",
       width = 25, height =6)


#_______________________________________________________________________________
#                          Correlation Scatterplots                        ---- 
#_______________________________________________________________________________
load("input_files/scaled_&_imputed_abundance.RData")
load("input_files/cytokine_profiles_no_stimulant.RData")
load("input_files/behaviors.RData")

keys2join <- read_xlsx('input_files/misc/sample_keys.xlsx', sheet = 'keys') %>%
  select(SampleID, group, treatment) %>% 
  column_to_rownames(var = "SampleID")

# Cytokine scatterplot loop ----
# Selects top 5 most significant associations per tissue source (of mets)
cytokine_hits <- 
  cytokine_cors %>% 
  filter(q < 0.25) %>%
  filter(n > 15) %>%
  group_by(tissue_B) %>%
  slice_min(order_by = q, n = 5)

cytokine_hits_sametissue <- 
  cytokine_cors %>% 
  filter(tissue_A == tissue_B) %>%
  filter(q < 0.25) %>%
  filter(n > 15) %>%
  group_by(tissue_B) %>%
  slice_min(order_by = q, n = 5)

cytokine_data <-
  cbind(
    cytokine_profiles_noStim$GI_cytokines,
    cytokine_profiles_noStim$Brain_cytokines %>% 
      clean.cols.y4() %>%
      clean.cols.upper(),
    cytokine_profiles_noStim$Plasma_cytokines %>% 
      clean.cols.y4() %>%
      clean.cols.upper())
merged_cyt_data <- cbind(keys2join, df_all, cytokine_data)


top_n_scatterplots(merged_data = merged_cyt_data, 
                   hitslist = cytokine_hits, 
                   data_type = "PBMC_profile_")

top_n_scatterplots(merged_data = merged_cyt_data, 
                   hitslist = cytokine_hits_sametissue, 
                   data_type = "PBMC_profile_",
                   folder_prefix = "figures/correlations/scatterplots_tissue_specific/")


# TO MANUALLY CHECK SPECIFIC CORRELATIONS:  
cytokine_hits_manual <- 
  cytokine_cors %>% 
  filter(tissue_B == "Ileum") %>%
  filter(feature_A == "HIPP_IL5") %>% 
  filter(q <= 0.2)

top_n_scatterplots(merged_data = merged_cyt_data, 
                   hitslist = cytokine_hits_manual, 
                   data_type = "PBMC_profile_",
                   folder_prefix = "figures/correlations/manual_scatterplots/")


# Behavior scatterplot loop ----
behavior_hits <- 
  behavior_all_mets_fdr %>% 
  # filter(q < 0.25) %>%
  group_by(tissue_B) %>%
  slice_min(order_by = q, n = 5)

merged_behav_data <- cbind(keys2join, df_all, behavior_data$df_behavior)
top_n_scatterplots(merged_data = merged_behav_data, 
                   hitslist = behavior_hits, 
                   data_type = "Behavior_")

#_______________________________________________________________________________
#                Pathway Summary Stacked Barcharts            ---- 
#_______________________________________________________________________________

cytokine_cors_mapped <- cytokine_cors %>%
  dplyr::rename("BIOCHEMICAL" = "feature_B_obj") %>%
  left_join(pathway_map, by = "BIOCHEMICAL") %>%
  janitor::clean_names() 

corr_count_sup <- 
  cytokine_cors_mapped %>% 
  filter(q < 0.25) %>%
  dplyr::group_by(feature_a_obj, super_pathway, tissue_b) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))
corr_count_sub <- 
  cytokine_cors_mapped %>% 
  filter(q < 0.25) %>%
  dplyr::group_by(feature_a_obj, sub_pathway, tissue_b) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))

level_order_sup <- 
  corr_count_sup %>% 
  group_by(super_pathway) %>% 
  dplyr::summarise(n = n()) %>%
  arrange(n) %>% 
  select(super_pathway) %>% unlist(use.names = F)

cytokine_colors <- 
  colorRampPalette(colormash)(length(base::unique(
    corr_count_sub$feature_a_obj)))
names(cytokine_colors) <- unique(
  corr_count_sub$feature_a_obj)

# Count figure - Superpathways
supcount <- 
  corr_count_sup %>% 
  ggplot(aes(x=n, y=fct_relevel(super_pathway, level_order), fill = feature_a_obj)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = cytokine_colors) +
  facet_wrap(~tissue_b, nrow = 1) +
  labs(x = "Count of Spearman's Correlation Associations (q<0.25)",
       y = NULL, fill = NULL) +
  my_clean_theme() 

ggsave(supcount, filename = glue('{path_enrich_loc}/2022-07-05_superpathway_cytokine_summary.png'),
       width = 9, height = 4)


# Count figure - Subpathways
subcount <- 
  corr_count_sub %>% 
  left_join(pathway_map_sub2sup, by = "sub_pathway") %>% 
  ggplot(aes(x=n, y=fct_relevel(sub_pathway, level_order), fill = feature_a_obj)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = cytokine_colors) +
  facet_grid(super_pathway~tissue_b, scales = "free_y", space = "free") +
  labs(x = "Count of Spearman's Correlation Associations (q<0.25)",
       y = NULL, fill = NULL) +
  my_clean_theme()

ggsave(subcount, filename = glue('{path_enrich_loc}/2022-07-05_subpathway_cytokine_summary.png'),
       width = 12, height = 12)


#_______________________________________________________________________________
#                    Pathway Enrichment Analysis            ---- 
#_______________________________________________________________________________

enrichment_df <- tibble()
for (metadata_var in unique(cytokine_cors_mapped$feature_a)) {
  for (tissue in unique(cytokine_cors_mapped$tissue_b)) {
    
    strat1 <- cytokine_cors_mapped %>%
      filter(tissue_b == tissue) %>%
      filter(feature_a == metadata_var) %>%
      distinct_at(vars(feature_a, feature_b, sub_pathway, tissue_b))
    strat1_sig <- cytokine_cors_mapped %>%
      filter(q < 0.25) %>%
      filter(tissue_b == tissue) %>%
      filter(feature_a == metadata_var) %>%
      distinct_at(vars(feature_a, feature_b, sub_pathway, tissue_b))
    
    for (path in unique(cytokine_cors_mapped$sub_pathway)) {
      total_lib <- strat1 %>% nrow()
      total_sig_lib <- strat1_sig %>% nrow()
      path_lib <- strat1 %>%
        filter(sub_pathway ==  path) %>%
        nrow()
      path_sig_lib <- strat1_sig %>%
        filter(sub_pathway ==  path) %>%
        nrow()
      
      enrichment_value <-
        enrichment_formula(N = total_lib,
                           n = total_sig_lib,
                           m = path_lib,
                           k = path_sig_lib)
      
      row2add <-
        cbind(
          "cytokine" = metadata_var,
          "tissue" = tissue,
          "path" = path,
          "N" = total_lib,
          "n" = total_sig_lib,
          "m" = path_lib,
          "k" = path_sig_lib,
          "enrichment" = enrichment_value
        )
      
      enrichment_df <- rbind(enrichment_df, row2add)
      cat(path, "enrichment_value: ", enrichment_value, "\n\n")
    }
  }
}

statvars <- c("N", "n", "m", "k", "enrichment")
enrichment_df <-
  enrichment_df %>%
  mutate(across(all_of(statvars), as.character),
         across(all_of(statvars), as.numeric))
enrichment_df <- enrichment_df %>% 
  dplyr::rename(sub_pathway = path) %>% 
  left_join(pathway_map_sub2sup, by = "sub_pathway")

# # Wrangling data format
cor_enrichment_df <- 
  enrichment_df %>% 
  dplyr::rename(metabolite_origin = tissue) %>% 
  separate(cytokine, into = c('cytokine', 'cytokine_origin'), sep = '__')

# # Wrangling data format
# cor_enrichment_df <- 
#   enrichment_df %>% 
#   mutate(cytokine = if_else(grepl("M_", cytokine), 
#                                   paste0(gsub("M_", "", cytokine), "__Plasma"),  
#                                   cytokine)) %>% 
#   mutate(cytokine = if_else(grepl("GCSF_", cytokine),
#                             paste0(gsub("GCSF_", "protect_gransf_", cytokine)),
#                             cytokine)) %>%
#   mutate(cytokine = if_else(grepl("G-CSF_", cytokine),
#                             paste0(gsub("G-CSF_", "protect_gransf_", cytokine)),
#                             cytokine)) %>% 
#   mutate(cytokine = if_else(grepl("GM-CSF_", cytokine),
#                             paste0(gsub("GM-CSF_", "protect_monf_", cytokine)),
#                             cytokine)) %>% 
#   mutate(cytokine = if_else(grepl("GMCSF_", cytokine),
#                             paste0(gsub("GMCSF_", "protect_monf_", cytokine)),
#                             cytokine)) %>%
#   mutate(cytokine = if_else(grepl("CSF_", cytokine),
#                             paste0(gsub("CSF_", "", cytokine), "__CSF"),
#                             cytokine)) %>%
#   mutate(cytokine = if_else(grepl("protect_gransf_", cytokine),
#                             paste0(gsub("protect_gransf_", "G-CSF_", cytokine)),
#                             cytokine)) %>% 
#   mutate(cytokine = if_else(grepl("protect_monf_", cytokine),
#                             paste0(gsub("protect_monf_", "GM-CSF_", cytokine)),
#                             cytokine)) %>% 
#   mutate(cytokine = if_else(grepl("AC_", cytokine), 
#                             paste0(gsub("AC_", "", cytokine), "__Anterior Cingulate"),  # DOUBLE CHECK
#                             cytokine)) %>% 
#   mutate(cytokine = if_else(grepl("HIPP_", cytokine), 
#                             paste0(gsub("HIPP_", "", cytokine), "__Hippocampus"),  
#                             cytokine)) %>% 
#   mutate(cytokine = if_else(grepl("PLASMA_", cytokine), 
#                             paste0(gsub("PLASMA_", "", cytokine), "__Brain Plasma"),  
#                             cytokine)) %>% 
#   mutate(cytokine_feature = if_else(grepl("CORTEX_", cytokine), 
#                                     paste0(gsub("CORTEX_", "", cytokine), "__Cortex"),  
#                                   cytokine)) %>% 
#   separate(cytokine_feature, c("cytokine", "cytokine_origin"), sep = "__", remove = F) %>% 
#   mutate(cytokine = toupper(cytokine)) %>% 
#   mutate(cytokine = str_replace(cytokine, "IL-", "IL")) %>% 
#   mutate(cytokine = str_replace(cytokine, "GCSF", "G-CSF")) %>% 
#   mutate(cytokine = str_replace(cytokine, "IL12P40", "IL12p40")) %>% 
#   mutate(cytokine = str_replace(cytokine, "MCP1", "MCP-1")) %>% 
#   mutate(cytokine = str_replace(cytokine, "MIP-1A", "MIP1A")) %>% 
#   mutate(cytokine = str_replace(cytokine, "MCP1", "MCP-1")) %>% 
#   mutate(cytokine = str_replace(cytokine, "SCD40L", "sCD40L")) %>% 
#   dplyr::rename(metabolite_origin = tissue)

# write.xlsx(cor_enrichment_df, 
#            file = 'data/correlations_cytokine_metabolite/2022-07-08_cytokine_correlation_pathway-enrichment-analysis.xlsx')



# Enrichment plots ----
enrichment_df_trim <- cor_enrichment_df %>% 
  filter(k > 0) %>%
  filter(n > 2) %>%
  filter(enrichment < 0.05)
enrichment_df_nonsig <-
  cor_enrichment_df %>% 
  filter(enrichment > 0.05)
enrichment_df_labels <-
  enrichment_df_trim %>% 
  group_by(metabolite_origin) %>% 
  slice_min(order_by = m/k, n = 5) %>% 
  slice_min(order_by = enrichment, n = 5)



# Function to display bubble plots
pathway_enrichment_bubbleplot <- function(cor_enrichment_df, tissue, ynudge = 0.5){
  
  # Tissue Specific bubble plot ----
  enrichment_df_trim <- cor_enrichment_df %>% 
    filter(metabolite_origin == tissue) %>% 
    filter(k > 0) %>%
    filter(n > 2) %>%
    filter(enrichment < 0.05)
  enrichment_df_nonsig <-
    cor_enrichment_df %>% 
    filter(metabolite_origin == tissue) %>% 
    filter(enrichment > 0.05)
  enrichment_df_labels <-
    enrichment_df_trim %>% 
    group_by(metabolite_origin) %>% 
    slice_min(order_by = m/k, n = 5) %>% 
    slice_min(order_by = enrichment, n = 5) %>% 
    mutate(lab = paste0(cytokine, " (", cytokine_origin, ")\n", sub_pathway))
  
  plot <- enrichment_df_nonsig %>% 
    ggplot() +
    geom_point(aes(y=k/m, x = -log10(enrichment), size =k), 
               fill = "#525252", alpha = 0.05 ) +
    geom_point(data = enrichment_df_trim, 
               aes(y=k/m, x = -log10(enrichment ), size = k, 
                   fill = super_pathway), alpha = 0.7, shape=21, color = "#393939", stroke = 0.25) +
    facet_wrap(~metabolite_origin, scales = "free") +
    scale_fill_manual(values = super_pathway_cols) +
    labs(y="Significant Pathway Elements (k) /\n Total Pathway Elements (m)",
         x = expression(paste(-log[10], "[ Enrichment Statistic ]")), 
         fill = "Super-pathway") +
    my_clean_theme() +
    scale_y_continuous(limits = c(0, 1)) +
    guides(size = guide_legend(order = 1),
           fill = guide_legend(order = 2)) + 
    geom_text_repel(data = enrichment_df_labels,
                    aes(y=k/m, x = -log10(enrichment ), label = lab),
                    segment.color	 =  "#929292",
                    segment.curvature = -0.5,
                    angle = 90,
                    nudge_y = ynudge,
                    direction = "x", 
                    vjust = 0.1,
                    force = 3, 
                    size = 2.75)
  print(plot);return(plot)
}

plasma_bubble <-
  pathway_enrichment_bubbleplot(cor_enrichment_df = cor_enrichment_df,
                                tissue = "Plasma",
                                ynudge = 0.7)
csf_bubble <-
  pathway_enrichment_bubbleplot(cor_enrichment_df = cor_enrichment_df, 
                                tissue = "CSF")
ileum_bubble <-
  pathway_enrichment_bubbleplot(cor_enrichment_df = cor_enrichment_df, 
                                tissue = "Ileum")
jejunum_bubble <-
  pathway_enrichment_bubbleplot(cor_enrichment_df = cor_enrichment_df, 
                                tissue = "Jejunum")
colon_bubble <-
  pathway_enrichment_bubbleplot(cor_enrichment_df = cor_enrichment_df, 
                                tissue = "Colon")
feces_bubble <-
  pathway_enrichment_bubbleplot(cor_enrichment_df = cor_enrichment_df, 
                                tissue = "Feces")


path_enrich_loc <- c('figures/correlations/pathway_enrichment')

ggsave(plasma_bubble, 
       filename = glue('{path_enrich_loc}/2022-07-05_bubble_plot_Plasma.svg'),
       width = 6.5, height = 5.75)
ggsave(csf_bubble, 
       filename = glue('{path_enrich_loc}/022-07-05_bubble_plot_CSF.svg'),
       width = 5.4, height = 5.75)
ggsave(ileum_bubble, 
       filename = glue('{path_enrich_loc}/022-07-05_bubble_plot_Ileum.svg'),
       width = 6.25, height = 5.75)
ggsave(jejunum_bubble, 
       filename = glue('{path_enrich_loc}/022-07-05_bubble_plot_Jejunum.svg'),
       width = 6.25, height = 5.75)
ggsave(colon_bubble, 
       filename = glue('{path_enrich_loc}/022-07-05_bubble_plot_Colon.svg'),
       width = 6, height = 5.75)
ggsave(feces_bubble, 
       filename = glue('{path_enrich_loc}/022-07-05_bubble_plot_Feces.svg'),
       width = 6, height = 5.75)





# figure - Superpathway enrichment barplots ----
sup_enrichment <- 
  cor_enrichment_df %>% 
  filter(k > 0) %>%
  filter(n > 2) %>% 
  ggplot(aes(x=-log10(enrichment), y=fct_relevel(super_pathway, level_order), fill = cytokine)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = cytokine_colors) +
  facet_wrap(~metabolite_origin, nrow = 1) +
  labs(x = expression(paste(-log[10], "[ Enrichment Statistic ]")),
       y = NULL, fill = NULL) +
  my_clean_theme() 

ggsave(sup_enrichment, 
       filename = 
         glue('{path_enrich_loc}/2022-07-05_superpathway_cytokine_enrichment_summary.png'),
       width = 9, height = 4)


# figure - Subpathways enrichment barplots ----
sub_enrichment <- 
  cor_enrichment_df %>% 
  filter(super_pathway !=  "Partially Characterized Molecules") %>% 
  filter(k > 0) %>%
  filter(n > 2) %>% 
  ggplot(aes(x=-log10(enrichment), y=fct_relevel(sub_pathway, level_order), fill = cytokine)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = cytokine_colors) +
  facet_grid(super_pathway~metabolite_origin, scales = "free_y", space = "free") +
  labs(x = expression(paste(-log[10], "[ enrichment statistic ]")),
       y = NULL, fill = NULL) +
  my_clean_theme() +
  theme(strip.text.y.right = element_text(angle = 0))

ggsave(sub_enrichment, 
filename = glue('{path_enrich_loc}/2022-07-05_subpathway_cytokine_enrichment_summary_new.svg'),
       width = 12, height = 12)



tissue_enrichment_subpaths <- function(df, tissue){
  plot <- df %>% 
    filter(metabolite_origin == tissue) %>% 
    filter(k > 0) %>%
    filter(n > 2) %>% 
    mutate(sub_pathway = forcats::fct_reorder(sub_pathway, -log10(enrichment), .fun = sum)) %>%
    ggplot(aes(x=-log10(enrichment), y=sub_pathway, fill = cytokine)) +
    geom_col(width = 0.5) +
    scale_fill_manual(values = cytokine_colors) +
    labs(x = expression(paste(-log[10], "[ enrichment statistic ]")),
         y = NULL, fill = NULL) +
    my_clean_theme()
  print(plot)
  return(plot)
}


plasma_subpath_enriched <- 
  tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "Plasma")

csf_subpath_enriched <- 
  tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "CSF")

jejunum_subpath_enriched <- 
  tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "Jejunum")

ileum_subpath_enriched <- 
  tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "Ileum")

colon_subpath_enriched <- 
  tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "Colon")

feces_subpath_enriched <- 
  tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "Feces")

ggsave(plasma_subpath_enriched, 
       filename = glue('{path_enrich_loc}/2022-07-05_subpath_barplot_Plasma.svg'),
       width = 7, height = 3)
ggsave(csf_subpath_enriched, 
       filename = glue('{path_enrich_loc}/2022-07-05_subpath_barplot_CSF.svg'),
       width = 7, height = 1.5)
ggsave(jejunum_subpath_enriched, 
       filename = glue('{path_enrich_loc}/2022-07-05_subpath_barplot_Jejunum.svg'),
       width = 7, height = 7)
ggsave(ileum_subpath_enriched, 
       filename = glue('{path_enrich_loc}/2022-07-05_subpath_barplot_Ileum.svg'),
       width = 7, height = 9)
ggsave(colon_subpath_enriched, 
       filename = glue('{path_enrich_loc}/2022-07-05_subpath_barplot_Colon.svg'),
       width = 7, height = 2.5)
ggsave(feces_subpath_enriched, 
       filename = glue('{path_enrich_loc}/2022-07-05_subpath_barplot_Feces.svg'),
       width = 7, height = 4)














# cytokine_cors_mapped %>% 
#   dplyr::group_by(feature_a_obj, super_pathway, tissue_b) %>%
#   dplyr::summarise(n = n()) %>%
#   arrange(desc(n))


# # Select for correlations to Plasma and identical metabolites across tissues
# df_cor_all2 <- df_cor_all %>% 
#   drop_na(P) %>% 
#   filter(tissue_A == "Plasma",
#          feature_A == feature_B) %>% 
#   left_join(pathway_map_join, by = "feature_A")
# 
# 
# df_cor_all2 %>% 
#   group_by(`SUPER PATHWAY`, tissue_B) %>% 
#   dplyr::summarise(super_path_mean = mean(rho)) %>% 
#   dplyr::rename("feature_B" = "SUPER PATHWAY") %>% 
#     ggplot(aes(x=fct_reorder(feature_B, super_path_mean, .desc = T), y=super_path_mean,
#                group = tissue_B)) +
#   geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=0.25, ymax=Inf), alpha = 0.01, fill = "grey") +
#   geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=-0.25, ymax=-Inf), alpha = 0.01, fill = "grey") +
#   geom_line(aes(color = tissue_B), size=1.5) +
#   labs(title = "Superpath correlations between Plasma levels", 
#        x = NULL, y = "Spearman's Rho", color = "Tissue") +
#   scale_color_npg() +
#   my_clean_theme() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# cor_avg_summary_subpath <-
#   df_cor_all2 %>% 
#   group_by(`SUPER PATHWAY`, `SUB PATHWAY`, tissue_B) %>% 
#   dplyr::summarise(super_path_mean = mean(rho)) %>% 
#   dplyr::rename("feature_B" = "SUB PATHWAY") %>% 
#   ggplot(aes(x=fct_reorder(feature_B, super_path_mean, .desc = T), y=super_path_mean)) +
#   geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=0.25, ymax=Inf), alpha = 0.01, fill = "grey") +
#   geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=-0.25, ymax=-Inf), alpha = 0.01, fill = "grey") +
#   geom_line(aes(color = tissue_B, group = tissue_B), size=1.5) +
#   geom_point() +
#   labs(title = "Subpath correlations between Plasma levels", 
#        x = NULL, y = "Spearman's Rho", color = "Tissue") +
#   # facet_grid(cols = vars(`SUPER PATHWAY`), scales = "free_x", space = "free_x") +
#   # facet_wrap(~`SUPER PATHWAY`, ncol = 1, scales = "free_x") +
#   scale_color_npg() +
#   my_clean_theme() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# cor_avg_summary_subpath
# 
# ggsave(cor_avg_summary_subpath, filename = 
#          "figures/correlations/2022-07-05_cor_avg_summary_subpath.png", 
#        width = 14, height = 8)
# 
# 
# plasma_analyte_corrs_v1 <- 
#   df_cor_all2 %>% 
#   janitor::clean_names() %>% 
#   mutate(xcol = paste(super_pathway, tissue_b, sep = "__")) %>% 
#   ggplot(aes(x=fct_reorder(xcol, rho, .desc = T), y=rho)) +
#   geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=0.25, ymax=Inf), fill = "#e9e9e9") +
#   geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=-0.25, ymax=-Inf), fill = "#e9e9e9") +
#   geom_point(aes(color = tissue_b), alpha = 0.2, position =
#                position_jitterdodge(jitter.height = 0, jitter.width = 0.2, seed = 42)) +
#   geom_boxplot(width = 0.1, outlier.alpha = 0, size = 0.3) +
#   geom_violin(aes(color = tissue_b), alpha = 0.1, size = 0.3) +
#   labs(x = NULL, y = "Spearman's Rho", color = "Tissue") +
#   facet_wrap(~super_pathway, scales = "free_x") +
#   scale_color_npg() +  scale_fill_npg() +
#   my_clean_theme() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x.bottom = element_blank())
# plasma_analyte_corrs_v1
# 
# ggsave(plasma_analyte_corrs_v1, 
#        filename = "figures/correlations/plasma_correlation_distributions_v1.svg",
#        width = 8, height = 6)
# 
# plasma_analyte_corrs_v2 <-
#   df_cor_all2 %>% 
#   janitor::clean_names() %>% 
#   mutate(xcol = paste(super_pathway, tissue_b, sep = "__")) %>% 
#   ggplot(aes(x=fct_reorder(xcol, rho, .desc = T), y=rho)) +
#   geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=0.25, ymax=Inf), fill = "#e9e9e9") +
#   geom_rect(aes(xmin=Inf, xmax=-Inf, ymin=-0.25, ymax=-Inf), fill = "#e9e9e9") +
#   geom_boxplot(width = 0.1, outlier.alpha = 0, size = 0.3) +
#   geom_violin(aes(color = tissue_b), alpha = 0.1, size = 0.3) +
#   labs(x = NULL, y = "Spearman's Rho", color = "Tissue") +
#   facet_wrap(~super_pathway, scales = "free_x") +
#   scale_color_npg() +  scale_fill_npg() +
#   my_clean_theme() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x.bottom = element_blank())
# plasma_analyte_corrs_v2
# 
# ggsave(plasma_analyte_corrs_v2, 
#        filename = "figures/correlations/plasma_correlation_distributions_v2.svg",
#        width = 8, height = 6)
