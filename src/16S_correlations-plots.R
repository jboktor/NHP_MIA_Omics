# Joe Boktor
# Caltech - Mazmanian Lab
# 2020

source("src/_load_packages.R")
source("src/_misc_functions.R")
source("src/_16S_functions.R")


# Load microbe cytokine correlations
microbe_cytokine_cors <- readRDS(file = "data/correlations/microbe-to-cytokine/Microbe_Cytokine_Correlations.rds")
# Load microbe metabolite correlations
microbe_metabolite_cors <- readRDS(file = "data/correlations/microbe-to-metabolite/Microbe_Metabolite_Correlations.rds")
# Load metabolome data
df_scaled_imp <- readRDS(file = "input_files/Metabolomics/scaled_&_imputed_abundance_annotated.rds")

keys <- read_xlsx("input_files/misc/sample_keys.xlsx", sheet = "keys") %>%
  select(SampleID, group, treatment) %>%
  column_to_rownames(var = "SampleID")
pathway_map <- df_scaled_imp %>%
  select(BIOCHEMICAL, `SUPER PATHWAY`, `SUB PATHWAY`) %>%
  distinct() %>%
  janitor::clean_names()
asv_labels <- ps_obj %>%
  tax_table() %>%
  as.data.frame() %>%
  unite("asv_label", Family:Species) %>%
  rownames_to_column(var = "feature") %>%
  mutate(
    asv_label = gsub("_NA", "", asv_label),
    asv_label_unique = paste(feature %>% substr(0, 6), asv_label, sep = "_")
  )
asv_labels %>% glimpse()


# _______________________________________________________________________________
# Quality Control Vis ----

microbe_metabolite_cors_QC1 <-
  microbe_metabolite_cors %>%
  left_join(asv_labels %>% dplyr::rename(feature_A = feature)) %>%
  left_join(pathway_map %>% dplyr::rename(feature_B_obj = biochemical)) %>%
  ggplot(aes(x = fct_reorder(feature_A, q), y = q, color = super_pathway)) +
  geom_point() +
  geom_ysidedensity(aes(x = stat(density))) +
  facet_grid(~tissue_B, space = "free", scales = "free") +
  scale_color_manual(values = super_pathway_cols) +
  labs(x = "ASVs") +
  my_clean_theme() +
  theme(axis.text.x = element_blank())

ggsave(microbe_metabolite_cors_QC1,
  filename = "figures/correlations/QC/microbe_metabolite_cors_QC1.svg",
  width = 15, height = 7
)


microbe_cytokine_cors_QC1 <-
  microbe_cytokine_cors %>%
  mutate(condition = factor(condition, levels = c("Media", "ConA", "PolyIC", "PMAiono"))) %>%
  left_join(asv_labels %>% dplyr::rename(feature_A = feature)) %>%
  ggplot(aes(x = fct_reorder(feature_A, q), y = q, color = feature_B_obj)) +
  geom_point() +
  geom_ysidedensity(aes(x = stat(density))) +
  facet_grid(tissue ~ condition) +
  scale_color_manual(values = colormash) +
  labs(x = "ASVs") +
  my_clean_theme() +
  theme(axis.text.x = element_blank())
ggsave(microbe_cytokine_cors_QC1,
  filename = "figures/correlations/QC/microbe_cytokine_cors_QC1.svg",
  width = 15, height = 7
)
microbe_cytokine_cors_QC2 <-
  microbe_cytokine_cors %>%
  mutate(condition = factor(condition, levels = c("Media", "ConA", "PolyIC", "PMAiono"))) %>%
  left_join(asv_labels %>% dplyr::rename(feature_A = feature)) %>%
  ggplot(aes(x = rho, y = q, color = feature_B_obj)) +
  geom_point() +
  geom_ysidedensity(aes(x = stat(density))) +
  geom_xsidedensity(aes(y = stat(density))) +
  facet_grid(tissue ~ condition) +
  scale_color_manual(values = colormash) +
  my_clean_theme() +
  theme(axis.text.x = element_blank())
ggsave(microbe_cytokine_cors_QC2,
  filename = "figures/correlations/QC/microbe_cytokine_cors_QC2.svg",
  width = 15, height = 7
)


# _______________________________________________________________________________
# ASV x Metabolite Heatmaps ----

ASVxMetabolites_stool <- microbe_metabolite_cors %>%
  filter(tissue_B == "Feces") %>%
  left_join(asv_labels %>% dplyr::rename(feature_A = feature)) %>%
  corr_heatmap_16s_to_mets(featselection = "q") + theme(plot.margin = unit(c(1, 1, 1, 5), "cm"))
ASVxMetabolites_jej <- microbe_metabolite_cors %>%
  filter(tissue_B == "Jejunum") %>%
  left_join(asv_labels %>% dplyr::rename(feature_A = feature)) %>%
  corr_heatmap_16s_to_mets(featselection = "q") + theme(plot.margin = unit(c(1, 1, 1, 7), "cm"))
ASVxMetabolites_ile <- microbe_metabolite_cors %>%
  filter(tissue_B == "Ileum") %>%
  left_join(asv_labels %>% dplyr::rename(feature_A = feature)) %>%
  corr_heatmap_16s_to_mets(featselection = "q") + theme(plot.margin = unit(c(1, 1, 1, 7), "cm"))

ggsave(ASVxMetabolites_stool,
  filename = "figures/correlations/heatmaps_by_tissue/heatmap_16S_ASVsxMetabolites_stool.svg",
  width = 7, height = 7
)
ggsave(ASVxMetabolites_jej,
  filename = "figures/correlations/heatmaps_by_tissue/heatmap_16S_ASVsxMetabolites_Jejunum.svg",
  width = 9, height = 7.5
)
ggsave(ASVxMetabolites_ile,
  filename = "figures/correlations/heatmaps_by_tissue/heatmap_16S_ASVsxMetabolites_Ileum.svg",
  width = 10, height = 7
)


# _______________________________________________________________________________
# ASV x Cytokine Heatmaps ----


ASVxCytokines_stool <- microbe_cytokine_cors %>%
  filter(
    tissue == "StoolxColon",
    condition == "Media"
  ) %>%
  left_join(asv_labels %>% dplyr::rename(feature_A = feature)) %>%
  corr_heatmap_16s_to_cyts(trimfeats = F)
ASVxCytokines_jej <- microbe_cytokine_cors %>%
  filter(
    tissue == "Ileum",
    condition == "Media"
  ) %>%
  left_join(asv_labels %>% dplyr::rename(feature_A = feature)) %>%
  corr_heatmap_16s_to_cyts(trimfeats = F)
ASVxCytokines_ile <- microbe_cytokine_cors %>%
  filter(
    tissue == "Jejunum",
    condition == "Media"
  ) %>%
  left_join(asv_labels %>% dplyr::rename(feature_A = feature)) %>%
  corr_heatmap_16s_to_cyts(trimfeats = F)

ggsave(ASVxCytokines_stool,
  filename = "figures/correlations/heatmaps_by_tissue/heatmap_16S_ASVsxCytokines_stool.svg",
  width = 7, height = 4
)
ggsave(ASVxCytokines_jej,
  filename = "figures/correlations/heatmaps_by_tissue/heatmap_16S_ASVsxCytokines_Jejunum.svg",
  width = 4, height = 4
)
ggsave(ASVxCytokines_ile,
  filename = "figures/correlations/heatmaps_by_tissue/heatmap_16S_ASVsxCytokines_Ileum.svg",
  width = 4, height = 4
)



# _______________________________________________________________________________
# ASV x Cytokine Heatmaps ----

# Load 16S data
ps_obj <- readRDS(file = "data/16S/phyloseq.rds")
# Load metabolome data
df_all <- readRDS(file = "input_files/Metabolomics/scaled_&_imputed_abundance.rds")
# Load GI cytokine data
gi_cytokines <- readRDS("input_files/Behavior_Immune_data/GI_cytokines_all.rds")



# _______________________________________________________________________________
#  Load data for Scatterplots  ----

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
# Scatterplots ----

masterlist <-
  list(
    "Feces" = list(ps_df_stool, stool_mets_input, colon_cyt_input[["Media"]]),
    "Jejunum" = list(ps_df_jej, jejunum_mets_input, jejunum_cyt_input[["Media"]]),
    "Ileum" = list(ps_df_ile, ileum_mets_input, ileum_cyt_input[["Media"]])
  )


for (tissue_origin in names(masterlist)) {
  met_hits <-
    microbe_metabolite_cors %>%
    filter(tissue_B == tissue_origin) %>%
    filter(q < 0.1) %>%
    slice_min(order_by = q, n = 20, with_ties = F)
  met_hits %>% glimpse()

  if (tissue_origin == "Feces") {
    tissue_origin_cyt <- "Colon"
  } else {
    tissue_origin_cyt <- tissue_origin
  }

  cyt_hits <-
    microbe_cytokine_cors %>%
    filter(
      tissue_B == tissue_origin_cyt,
      condition == "Media"
    ) %>%
    filter(q < 0.1) %>%
    slice_min(order_by = q, n = 20, with_ties = F)
  cyt_hits %>% glimpse()


  keys_filtered <- keys[rownames(masterlist[[tissue_origin]][[1]]), ]
  merged_test_df <- cbind(
    keys_filtered,
    masterlist[[tissue_origin]][[1]],
    masterlist[[tissue_origin]][[2]],
    masterlist[[tissue_origin]][[3]]
  )

  top_n_scatterplots(
    merged_data = merged_test_df,
    hitslist = met_hits,
    data_type = "16S_",
    subdir = F,
    folder_prefix = paste0(
      "figures/correlations/scatterplots/",
      tissue_origin,
      "/16SxMetabolites"
    )
  )

  top_n_scatterplots(
    merged_data = merged_test_df,
    hitslist = cyt_hits,
    data_type = "16S_",
    subdir = F,
    folder_prefix = paste0(
      "figures/correlations/scatterplots/",
      tissue_origin,
      "/16SxCytokines"
    )
  )
}








# #________________________________________________________________________________
# # sankeyNetwork Plot
#
# library(networkD3)
#
# microbe_mets_stool_df <- microbe_metabolite_cors %>%
#   filter(tissue_B == "Jejunum",
#          q < 0.05, abs(rho) > 0.4)
# microbe_cyt_stool_media_df <- microbe_cytokine_cors %>%
#   filter(tissue_B == "Jejunum",
#          condition == "Media",
#          q < 0.05, abs(rho) > 0.4)
#
# edges <- data.frame(
#   source=c(microbe_mets_stool_df %>% pull(feature_A),
#            microbe_cyt_stool_media_df %>% pull(feature_A)),
#   target=c(microbe_mets_stool_df %>% pull(feature_B_obj),
#            microbe_cyt_stool_media_df %>% pull(feature_B_obj)),
#   values_charged=c(microbe_mets_stool_df %>% pull(rho),
#                    microbe_cyt_stool_media_df %>% pull(rho)),
#   assessment=c(rep("metabolites", nrow(microbe_mets_stool_df)),
#                rep("cytokines", nrow(microbe_cyt_stool_media_df)))) %>%
#   mutate(
#     value = abs(values_charged),
#     direction = if_else(values_charged > 0, "positive", "negative"))
#
#
# nodes <- data.frame(name = c(as.character(edges$source),
#                              as.character(edges$target)) %>% unique()) %>%
#   left_join(asv_labels %>% dplyr::rename(name = feature)) %>%
#   left_join(pathway_map %>% dplyr::rename(name = biochemical)) %>%
#   mutate(group = if_else(
#     name %in% microbe_mets_stool_df$feature_A |
#       name %in% microbe_cyt_stool_media_df$feature_A, Phylum,
#     if_else(name %in% microbe_mets_stool_df$feature_B_obj, super_pathway,
#             if_else(name %in% microbe_cyt_stool_media_df$feature_B_obj, "cytokine", "error"
#             ))))
#
#
# # With networkD3, connection must be provided using id, not using real name like in the edges dataframe.. So we need to reformat it.
# edges$IDsource <- match(edges$source, nodes$name)-1
# edges$IDtarget <- match(edges$target, nodes$name)-1
#
# edges %>% glimpse()
# nodes %>% glimpse()
# colourScale <-
#   c(
#     "cytokine" = "grey60",
#     "positive" = "#4483AA",
#     "negative" = "#C287C0"
#   )
#
# taxacols <- brewer.pal(length(na.omit(unique(nodes$Phylum))), "Accent")
# names(taxacols) <- na.omit(unique(nodes$Phylum))
#
# colourScale <- c(colourScale, super_pathway_cols, taxacols)
# colourScaleStr = networkD3::JS(paste0("d3.scaleOrdinal().domain(['", paste(names(colourScale), collapse="', '"), "']).range(['", paste(colourScale, collapse="', '"), "'])"))
#
#
# # Make the Network
# p <- sankeyNetwork(
#   Links = edges,
#   Nodes = nodes,
#   Source = "IDsource",
#   Target = "IDtarget",
#   Value = "value",
#   NodeID = "name",
#   LinkGroup="direction",
#   NodeGroup="group",
#   sinksRight =F,
#   colourScale=colourScaleStr,
#   fontSize = 20,
#   iterations=10000,
# )
# p
