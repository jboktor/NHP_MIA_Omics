
source("src/_load_packages.R")
source("src/_misc_functions.R")

# Load data ----
pathway_map <-
  readRDS(file = "data/metabolomics/2022-07-05_pathway-map.rds")
pathway_map_join <-
  pathway_map %>% dplyr::rename("feature_A" = "BIOCHEMICAL")
pathway_map_sub2sup <- pathway_map %>%
  janitor::clean_names() %>%
  select(super_pathway, sub_pathway) %>%
  distinct()

cyt_to_mets <- "data/correlations/cytokine-to-metabolite/"
cytokine_cors <- readRDS(glue(
  "{cyt_to_mets}/2022-07-06_cytokine-metabolite_correlations.rds"
))

behavior_all_mets_fdr <-
  read_excel(
    "data/correlations/behavior-to-metabolite/Behavior_Metabolite_Correlations.xlsx",
    sheet = "behaviors_allMetabolites"
  )

cyt_tissues <-
  c(
    "Plasma-PBMCs", "Plasma", "CSF", "Cortex", "HIPP", "AC",
    "Jejunum-PBMCs", "Ileum-PBMCs", "MLN-PBMCs", "Colon-PBMCs"
  )

## partition metabolites by tissue ----
plasma_cyt_cors <-
  cytokine_cors %>%
  filter(tissue_B == "Plasma")

csf_cyt_cors <-
  cytokine_cors %>%
  filter(tissue_B == "CSF")

GI_cyt_cors <-
  cytokine_cors %>%
  filter(tissue_B %in% c("Jejunum", "Ileum", "Colon"))

feces_cyt_cors <-
  cytokine_cors %>%
  filter(tissue_B == "Feces")

# _______________________________________________________________________________
#                    Tissue source specific Heatmaps  ----
# _______________________________________________________________________________

plasma_cyt_cors.plot <- corr_heatmap_facet_v2(plasma_cyt_cors, facet_ord = cyt_tissues)
ggsave(plasma_cyt_cors.plot,
  filename = "figures/correlations/2022-07-05_heatmap_Plasma_allmets_allcytokines.pdf",
  width = 25, height = 5
)
plasma_cyt_cors.plot <- corr_heatmap_facet_v2(plasma_cyt_cors, facet_ord = cyt_tissues)
ggsave(plasma_cyt_cors.plot,
  filename = "figures/correlations/2022-07-05_heatmap_Plasma_allmets_allcytokines.pdf",
  width = 25, height = 5
)
csf_cyt_cors.plot <- corr_heatmap_facet_v2(csf_cyt_cors, facet_ord = cyt_tissues)
ggsave(csf_cyt_cors.plot,
  filename = "figures/correlations/2022-07-05_heatmap_CSF_allmets_allcytokines.pdf",
  width = 25, height = 4
)
GI_cyt_cors.plot <- corr_heatmap_facet_v2(GI_cyt_cors, facet_ord = cyt_tissues)
ggsave(GI_cyt_cors.plot,
  filename = "figures/correlations/2022-07-05_heatmap_GI_allmets_allcytokines.pdf",
  width = 25, height = 6
)
feces_cyt_cors.plot <- corr_heatmap_facet_v2(feces_cyt_cors, facet_ord = cyt_tissues)
ggsave(feces_cyt_cors.plot,
  filename = "figures/correlations/2022-07-05_heatmap_Feces_allmets_allcytokines.pdf",
  width = 25, height = 6
)

# _______________________________________________________________________________
#                          Correlation Scatterplots                        ----
# _______________________________________________________________________________

df_all <- readRDS("data/metabolomics/2022-07-05_scaled_imputed_abundance.rds")
behavior_data <- readRDS("data/behaviors/2022-07-05_behaviors.rds")
cytokine_data <- readRDS("data/cytokines/2022-07-05_cytokine-profiles_imputed.rds")
keys2join <- read_xlsx("input_files/misc/sample_keys.xlsx", sheet = "keys") %>%
  select(SampleID, group, treatment)

## Cytokine scatterplot loop ----
# Selects top 5 most significant associations per tissue source (of mets)
cytokine_hits <-
  cytokine_cors %>%
  filter(q < 0.25) %>%
  filter(n > 15) %>%
  group_by(tissue_B) %>%
  slice_min(order_by = q, n = 5)

cytokine_hits_sametissue <-
  cytokine_cors %>%
  mutate(temp_tissue_A = gsub('-PBMCs', '', tissue_A)) %>% 
  filter(temp_tissue_A == tissue_B) %>%
  select(-temp_tissue_A) %>% 
  filter(q < 0.25) %>%
  filter(n > 15) %>%
  group_by(tissue_B) %>%
  slice_min(order_by = q, n = 5)

cytokine_data_wide <-
  cytokine_data %>% 
  pivot_wider(names_from = c('cytokine', 'Tissue'), names_sep = '__')

merged_cyt_data <- 
  df_all %>% 
  rownames_to_column(var = 'SampleID') %>% 
  left_join(keys2join) %>% 
  left_join(cytokine_data_wide) %>% 
  column_to_rownames(var = 'SampleID')

top_n_scatterplots(
  merged_data = merged_cyt_data,
  hitslist = cytokine_hits,
  data_type = "PBMC_profile_", 
  folder_prefix = 'figures/correlations/scatterplots_2022-07-09/'
)

top_n_scatterplots(
  merged_data = merged_cyt_data,
  hitslist = cytokine_hits_sametissue,
  data_type = "PBMC_profile_",
  folder_prefix = "figures/correlations/scatterplots_tissue-specific_2022-07-09/"
)


# # TO MANUALLY CHECK SPECIFIC CORRELATIONS:
# cytokine_hits_manual <-
#   cytokine_cors %>%
#   filter(tissue_B == "Ileum") %>%
#   filter(feature_A == "HIPP_IL5") %>%
#   filter(q <= 0.2)
# 
# top_n_scatterplots(
#   merged_data = merged_cyt_data,
#   hitslist = cytokine_hits_manual,
#   data_type = "PBMC_profile_",
#   folder_prefix = "figures/correlations/manual_scatterplots/"
# )

## Behavior scatterplot loop ----
behavior_hits <-
  behavior_all_mets_fdr %>%
  filter(q < 0.25) %>%
  group_by(tissue_B) %>%
  slice_min(order_by = q, n = 5)

df_behavior <- 
  behavior_data$df_behavior %>% 
  rownames_to_column(var = 'SampleID')

merged_behav_data <- 
  df_all %>% 
  rownames_to_column(var = 'SampleID') %>% 
  left_join(keys2join) %>% 
  left_join(df_behavior) %>% 
  column_to_rownames(var = 'SampleID')

top_n_scatterplots(
  merged_data = merged_behav_data,
  hitslist = behavior_hits,
  data_type = "Behavior_"
)

# _______________________________________________________________________________
#                Pathway Summary Stacked Barcharts            ----
# _______________________________________________________________________________

cytokine_cors_mapped <- cytokine_cors %>%
  dplyr::rename("BIOCHEMICAL" = "feature_B_obj") %>%
  left_join(pathway_map, by = "BIOCHEMICAL") %>%
  janitor::clean_names()

# _______________________________________________________________________________
#                    Pathway Enrichment Analysis            ----
# _______________________________________________________________________________

# enrichment_df <- tibble()
# for (metadata_var in unique(cytokine_cors_mapped$feature_a)) {
#   for (tissue in unique(cytokine_cors_mapped$tissue_b)) {
#
#     strat1 <- cytokine_cors_mapped %>%
#       filter(tissue_b == tissue) %>%
#       filter(feature_a == metadata_var) %>%
#       distinct_at(vars(feature_a, feature_b, sub_pathway, tissue_b))
#     strat1_sig <- cytokine_cors_mapped %>%
#       filter(q < 0.25) %>%
#       filter(tissue_b == tissue) %>%
#       filter(feature_a == metadata_var) %>%
#       distinct_at(vars(feature_a, feature_b, sub_pathway, tissue_b))
#
#     for (path in unique(cytokine_cors_mapped$sub_pathway)) {
#       total_lib <- strat1 %>% nrow()
#       total_sig_lib <- strat1_sig %>% nrow()
#       path_lib <- strat1 %>%
#         filter(sub_pathway ==  path) %>%
#         nrow()
#       path_sig_lib <- strat1_sig %>%
#         filter(sub_pathway ==  path) %>%
#         nrow()
#
#       enrichment_value <-
#         enrichment_formula(N = total_lib,
#                            n = total_sig_lib,
#                            m = path_lib,
#                            k = path_sig_lib)
#
#       row2add <-
#         cbind(
#           "cytokine" = metadata_var,
#           "tissue" = tissue,
#           "path" = path,
#           "N" = total_lib,
#           "n" = total_sig_lib,
#           "m" = path_lib,
#           "k" = path_sig_lib,
#           "enrichment" = enrichment_value
#         )
#
#       enrichment_df <- rbind(enrichment_df, row2add)
#       cat(path, "enrichment_value: ", enrichment_value, "\n\n")
#     }
#   }
# }
#
# statvars <- c("N", "n", "m", "k", "enrichment")
# enrichment_df <-
#   enrichment_df %>%
#   mutate(across(all_of(statvars), as.character),
#          across(all_of(statvars), as.numeric))
# enrichment_df <- enrichment_df %>%
#   dplyr::rename(sub_pathway = path) %>%
#   left_join(pathway_map_sub2sup, by = "sub_pathway")
#
# # Wrangling data format
# cor_enrichment_df <-
#   enrichment_df %>%
#   dplyr::rename(metabolite_origin = tissue) %>%
#   separate(cytokine, into = c('cytokine', 'cytokine_origin'), sep = '__')

# saveRDS(cor_enrichment_df, file = 'data/cytokines/2022-07-08_cytokine-metabolite-corr-pathway-enrichment.rds')
# write.xlsx(cor_enrichment_df,
#            file = 'data/correlations/cytokine-to-metabolite/2022-07-08_cytokine_correlation_pathway-enrichment-analysis.xlsx')

cor_enrichment_df <- readRDS(
  file = "data/cytokines/2022-07-08_cytokine-metabolite-corr-pathway-enrichment.rds"
)



## Tissue Specific bubble plot ----
# Function to display bubble plots
pathway_enrichment_bubbleplot <- function(cor_enrichment_df, tissue, ynudge = 0.5) {
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
    slice_min(order_by = m / k, n = 5) %>%
    slice_min(order_by = enrichment, n = 5) %>%
    mutate(lab = paste0(cytokine, " (", cytokine_origin, ")\n", sub_pathway))

  plot <- enrichment_df_nonsig %>%
    ggplot() +
    geom_point(aes(y = k / m, x = -log10(enrichment), size = k),
      fill = "#525252", alpha = 0.05
    ) +
    geom_point(
      data = enrichment_df_trim,
      aes(
        y = k / m, x = -log10(enrichment), size = k,
        fill = super_pathway
      ), alpha = 0.7, shape = 21, color = "#393939", stroke = 0.25
    ) +
    facet_wrap(~metabolite_origin, scales = "free") +
    scale_fill_manual(values = super_pathway_cols) +
    labs(
      y = "Significant Pathway Elements (k) /\n Total Pathway Elements (m)",
      x = expression(paste(-log[10], "[ Enrichment Statistic ]")),
      fill = "Super-pathway"
    ) +
    my_clean_theme() +
    scale_y_continuous(limits = c(0, 1)) +
    guides(
      size = guide_legend(order = 1),
      fill = guide_legend(order = 2)
    ) +
    geom_text_repel(
      data = enrichment_df_labels,
      aes(y = k / m, x = -log10(enrichment), label = lab),
      segment.color = "#929292",
      segment.curvature = -0.5,
      angle = 90,
      nudge_y = ynudge,
      direction = "x",
      vjust = 0.1,
      force = 3,
      size = 2.75
    )
  print(plot)
  return(plot)
}

plasma_bubble <-
  pathway_enrichment_bubbleplot(
    cor_enrichment_df = cor_enrichment_df,
    tissue = "Plasma",
    ynudge = 0.7
  )
csf_bubble <-
  pathway_enrichment_bubbleplot(
    cor_enrichment_df = cor_enrichment_df,
    tissue = "CSF"
  )
ileum_bubble <-
  pathway_enrichment_bubbleplot(
    cor_enrichment_df = cor_enrichment_df,
    tissue = "Ileum"
  )
jejunum_bubble <-
  pathway_enrichment_bubbleplot(
    cor_enrichment_df = cor_enrichment_df,
    tissue = "Jejunum"
  )
colon_bubble <-
  pathway_enrichment_bubbleplot(
    cor_enrichment_df = cor_enrichment_df,
    tissue = "Colon"
  )
feces_bubble <-
  pathway_enrichment_bubbleplot(
    cor_enrichment_df = cor_enrichment_df,
    tissue = "Feces"
  )


path_enrich_loc <- c("figures/correlations/pathway_enrichment")

ggsave(plasma_bubble,
  filename = glue("{path_enrich_loc}/2022-07-05_bubble_plot_Plasma.pdf"),
  width = 6.5, height = 5.75
)
ggsave(csf_bubble,
  filename = glue("{path_enrich_loc}/022-07-05_bubble_plot_CSF.pdf"),
  width = 5.4, height = 5.75
)
ggsave(ileum_bubble,
  filename = glue("{path_enrich_loc}/022-07-05_bubble_plot_Ileum.pdf"),
  width = 6.25, height = 5.75
)
ggsave(jejunum_bubble,
  filename = glue("{path_enrich_loc}/022-07-05_bubble_plot_Jejunum.pdf"),
  width = 6.25, height = 5.75
)
ggsave(colon_bubble,
  filename = glue("{path_enrich_loc}/022-07-05_bubble_plot_Colon.pdf"),
  width = 6, height = 5.75
)
ggsave(feces_bubble,
  filename = glue("{path_enrich_loc}/022-07-05_bubble_plot_Feces.pdf"),
  width = 6, height = 5.75
)



## figure - Superpathway enrichment barplots ----
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

level_order_ <- 
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

## figure - Subpathways enrichment barplots ----
sub_enrichment <-
  cor_enrichment_df %>%
  filter(super_pathway != "Partially Characterized Molecules") %>%
  filter(k > 0) %>%
  filter(n > 2) %>%
  ggplot(aes(x = -log10(enrichment), 
             y = sub_pathway, 
             fill = cytokine)
         ) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = cytokine_colors) +
  facet_grid(super_pathway ~ metabolite_origin, scales = "free_y", space = "free") +
  labs(
    x = expression(paste(-log[10], "[ enrichment statistic ]")),
    y = NULL, fill = NULL
  ) +
  my_clean_theme() +
  theme(strip.text.y.right = element_text(angle = 0))

ggsave(sub_enrichment,
  filename = glue("{path_enrich_loc}/2022-07-05_subpathway_cytokine_enrichment_summary_new.pdf"),
  width = 12, height = 12
)


tissue_enrichment_subpaths <- function(df, tissue) {
  plot <- df %>%
    filter(metabolite_origin == tissue) %>%
    filter(k > 0) %>%
    filter(n > 2) %>%
    mutate(sub_pathway = forcats::fct_reorder(sub_pathway, -log10(enrichment), .fun = sum)) %>%
    ggplot(aes(x = -log10(enrichment), y = sub_pathway, fill = cytokine)) +
    geom_col(width = 0.5) +
    scale_fill_manual(values = cytokine_colors) +
    labs(
      x = expression(paste(-log[10], "[ enrichment statistic ]")),
      y = NULL, fill = NULL
    ) +
    my_clean_theme()
  print(plot)
  return(plot)
}


# plasma_subpath_enriched <-
#   tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "Plasma")
# 
# csf_subpath_enriched <-
#   tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "CSF")
# 
# jejunum_subpath_enriched <-
#   tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "Jejunum")
# 
# ileum_subpath_enriched <-
#   tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "Ileum")
# 
# colon_subpath_enriched <-
#   tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "Colon")
# 
# feces_subpath_enriched <-
#   tissue_enrichment_subpaths(df = cor_enrichment_df, tissue = "Feces")
# 
# ggsave(plasma_subpath_enriched,
#   filename = glue("{path_enrich_loc}/2022-07-05_subpath_barplot_Plasma.pdf"),
#   width = 7, height = 3
# )
# ggsave(csf_subpath_enriched,
#   filename = glue("{path_enrich_loc}/2022-07-05_subpath_barplot_CSF.pdf"),
#   width = 7, height = 1.5
# )
# ggsave(jejunum_subpath_enriched,
#   filename = glue("{path_enrich_loc}/2022-07-05_subpath_barplot_Jejunum.pdf"),
#   width = 7, height = 7
# )
# ggsave(ileum_subpath_enriched,
#   filename = glue("{path_enrich_loc}/2022-07-05_subpath_barplot_Ileum.pdf"),
#   width = 7, height = 9
# )
# ggsave(colon_subpath_enriched,
#   filename = glue("{path_enrich_loc}/2022-07-05_subpath_barplot_Colon.pdf"),
#   width = 7, height = 2.5
# )
# ggsave(feces_subpath_enriched,
#   filename = glue("{path_enrich_loc}/2022-07-05_subpath_barplot_Feces.pdf"),
#   width = 7, height = 4
# )
