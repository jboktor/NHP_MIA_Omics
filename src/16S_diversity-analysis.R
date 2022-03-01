# Joe Boktor
# Caltech - Mazmanian Lab
# 2020

source("src/_load_packages.R")
source("src/_misc_functions.R")
source("src/_16S_functions.R")

ps_obj_rare <- readRDS(file = "data/16S/phyloseq_rarified.rds")

#_______________________________________________________________________________
#                             Alpha-Diversity ----
#_______________________________________________________________________________


alpha_stool <- alpha_div_metrics(ps_obj_rare[["Stool"]])
alpha_jej <- alpha_div_metrics(ps_obj_rare[["Jejunum"]])
alpha_ile <- alpha_div_metrics(ps_obj_rare[["Ileum"]])

alpha_div_analysis_stool <- alpha_div_analysis(alpha_stool, info = "stool")
alpha_div_analysis_jej <- alpha_div_analysis(alpha_jej, info = "jejunum")
alpha_div_analysis_ile <- alpha_div_analysis(alpha_ile, info = "ileum")

alpha_diversity_stats <- c(
  list("metrics-stool" = alpha_stool),
  alpha_div_analysis_stool,
  list("metrics-jejunum" = alpha_jej),
  alpha_div_analysis_jej,
  list("metrics-ileum" = alpha_ile),
  alpha_div_analysis_ile
)

write.xlsx(x = alpha_diversity_stats, 
           file = 'data/16S/alpha_diversity_analysis.xlsx',
           overwrite = T)


#_______________________________________________________________________________
#                             Beta-Diversity ----
#_______________________________________________________________________________

# PCoA Visualization (Weighted and UnWeighted UniFrac)
pcoa_stool <- unifrac_plots(ps_obj_rare[["Stool"]])
pcoa_jej <- unifrac_plots(ps_obj_rare[["Jejunum"]])
pcoa_ile <- unifrac_plots(ps_obj_rare[["Ileum"]])

ggsave(pcoa_stool, filename = "figures/16S/beta-diversity/stool_unifrac_pcoa.svg", 
       width =  7, height = 3)
ggsave(pcoa_jej, filename = "figures/16S/beta-diversity/jejunum_unifrac_pcoa.svg", 
       width =  7, height = 3)
ggsave(pcoa_ile, filename = "figures/16S/beta-diversity/ileum_unifrac_pcoa.svg", 
       width =  7, height = 3)


#_______________________________________________________________________________
#                              PERMANOVA ----
#_______________________________________________________________________________


permanova_vars <-
  c("treatment",
    "gender",
    "tissue",
    "group")

permanova_stats <- list()
for (tissue_origin in c("Stool", "Jejunum", "Ileum")){
  for (d in c("Aitchisons", "uunifrac", "wunifrac")){
    permanova_stats[[tissue_origin]][[d]] <-
      phyloseq_permanova(
        ps_obj_rare[[tissue_origin]],
        metadata_list = permanova_vars,
        nperm = 1000, 
        dist = d) %>% 
      mutate(tissue = tissue_origin)
  }
}

permanova_stats_df <- do.call(bind_rows, unlist(permanova_stats, recursive=FALSE))
write.xlsx(x = permanova_stats_df, 
           file = 'data/16S/beta_diversity_analysis_PERMANOVA.xlsx', 
           overwrite = T)






