# Joe Boktor
# Caltech - Mazmanian Lab
# 2020

source("src/_load_packages.R")
source("src/_misc_functions.R")
source("src/_16S_functions.R")
ps_obj <- readRDS(file = "data/16S/phyloseq.rds")
# aucs <- readRDS(file = "data/16S/differential-abundance/feature_AUROCs_treatment.rds")

#_______________________________________________________________________________ 
# MaAsLin2 ----

asv_labels <- ps_obj %>% tax_table() %>% as.data.frame() %>% 
  unite('asv_label', Family:Species) %>% 
  rownames_to_column(var = "feature") %>% 
  mutate(asv_label = gsub("_NA", "", asv_label))

maaslin2_stats <- readRDS(file = 'data/16S/differential-abundance/MaAsLin2/MaAsLin2_summary.rds')
maaslin2_stats %>% glimpse()

maaslin2_stats %>%
  filter(
    tissue == "Stool",
    metadata == "treatment",
    feature_level == "ASVs",
    qval <= 0.1) %>% 
  nrow()

maaslin2_stats_all <- 
  maaslin2_stats %>% 
  filter(metadata == "treatment",
         feature_level == "ASVs",
         qval <= 0.1) %>% 
  left_join(asv_labels) %>% 
  left_join(aucs) %>% 
  group_by(tissue) %>%
  # slice_min(n = 30, order_by = qval, with_ties = F) %>% 
  ggplot(aes(y = reorder(asv_label, coef),
             x = coef,
             xmin = coef - stderr, 
             xmax = coef + stderr,
             color = Phylum)) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.4) +
  geom_errorbar(width=0.6, size = 0.6) +
  my_clean_theme() +
  facet_grid(rows = vars(tissue), scales = "free", space = "free") +
  labs(x = "Regression coefficient [Poly(I:C) / Control]", y = NULL) +
  scale_color_d3() +
  geom_point(size = 0.6) +
  theme(axis.text.y = element_text(size = 8),
        axis.ticks = element_blank())
maaslin2_stats_all


maaslin2_stats_stool <- 
  maaslin2_stats %>% 
  filter(metadata == "treatment",
         tissue == "Stool",
         feature_level == "ASVs",
         qval <= 0.1) %>% 
  left_join(asv_labels) %>% 
  left_join(aucs) %>% 
  group_by(tissue) %>%
  ggplot(aes(y = reorder(asv_label, coef),
             x = coef,
             xmin = coef - stderr, 
             xmax = coef + stderr,
             color = Phylum)) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.4) +
  geom_errorbar(width=0.6, size = 0.6) +
  my_clean_theme() +
  # facet_grid(rows = vars(tissue), scales = "free", space = "free") +
  labs(x = "Regression coefficient [Poly(I:C) / Control]", y = NULL) +
  scale_color_d3() +
  geom_point(size = 0.6) +
  theme(axis.text.y = element_text(size = 8),
        axis.ticks = element_blank())
maaslin2_stats_stool




ggsave(maaslin2_stats_all, 
       filename = "figures/16S/differential-abundance/maaslin2_stats_summary.svg", 
       width = 8, height = 7)




# plot_df <- maaslin2_stats %>% 
#   filter(metadata == "treatment",
#          tissue == "Stool",
#          feature_level == "ASVs") %>% 
#   left_join(asv_labels) %>% 
#   left_join(aucs)
# plot_df_nonsig <- plot_df %>% filter(qval > 0.01 | abs(auroc-0.5) < .1)
# plot_df_sig <- plot_df %>% filter(qval <= 0.01 & abs(auroc-0.5) > .1)
# 
# 
# 
# plot_df_nonsig %>% 
#   ggplot() +
#   geom_vline(xintercept = 0.5) +
#   geom_hline(yintercept = -log10(0.05)) +
#   geom_point(aes(x = auroc, y = -log10(qval + 1e-100)), color = "grey60") +
#   geom_point(data = plot_df_sig,
#              aes( x = auroc, y = -log10(qval + 1e-100), color = Phylum)) +
#   my_clean_theme() +
#   facet_wrap(~tissue, scales = "free") +
#   labs(x = "AUROC", y = NULL) +
#   scale_color_futurama() +
#   theme(axis.text.y = element_text(size = 8),
#         axis.ticks = element_blank())
# 
# 


  

