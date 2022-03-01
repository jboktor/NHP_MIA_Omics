# Joe Boktor
# Caltech - Mazmanian Lab
# 2020

source("src/_load_packages.R")
source("src/_misc_functions.R")
source("src/_16S_functions.R")
ps_obj <- readRDS(file = "data/16S/phyloseq.rds")

# load("input_files/16S_data/phyloseq_stool.RData") # dat.stool
# load("input_files/16S_data/phyloseq_ileum.RData")
# load("input_files/16S_data/phyloseq_jejunum.RData")

library(ANCOMBC)
library(corncob)
library(Maaslin2)
library(doParallel)


#_______________________________________________________________________________



maaslin2_stats <- list()

for (tissue_origin in c("Stool", "Jejunum", "Ileum")){
  ps_dat <- ps_obj %>% 
    subset_samples(tissue == tissue_origin) 
  
  ps_dat_list <- list(
    "ASVs"= ps_dat,
    "Phylum" = aggregate_taxa(ps_dat, "Phylum"),
    "Genus" = aggregate_taxa(ps_dat, "Genus"),
    "Species" = aggregate_taxa(ps_dat, "Species")
  )
  
  #setup parallel processing
  start_time <- Sys.time()
  cores = detectCores()
  cl <- makeCluster(cores[1] - 2) # to prevent computer overload
  registerDoParallel(cl)
  
  maaslin2_loop <- 
    foreach(
      level = names(ps_dat_list),
      .combine = 'rbind',
      .verbose = F,
      .packages = c("Maaslin2", "phyloseq", "microbiome",
                    "magrittr", "tibble", "dplyr")
    ) %dopar% {
      # pre-process phyloseq obj
      dat_test <- ps_dat_list[[level]] %>%
        microbiome::core(detection = 0, prevalence = 0.1)
      sample_data(dat_test)$treatment <-
        factor(sample_data(dat_test)$treatment,
               levels = c("Saline", "Poly(I:C)"))
      
      input_data_df <-
        data.frame(otu_table(dat_test)) %>% 
        rownames_to_column() %>% 
        make_rfriendly_rows("rowname")
      
      input_metadata_df <-
        data.frame(sample_data(dat_test))
      rownames(input_metadata_df) <-
        gsub(" ", ".", rownames(input_metadata_df))
      
      
      outdir <-
        paste0("data/16S/differential-abundance/MaAsLin2/", 
               tissue_origin, "_", level)
      
      
      maaslin2_mod <- Maaslin2(
        input_data = input_data_df,
        input_metadata = input_metadata_df,
        output = outdir,
        min_abundance = 0.0,
        min_prevalence = 0.0,
        normalization = "TMM",
        transform = "NONE",
        analysis_method = "NEGBIN",
        max_significance = 0.05,
        fixed_effects = c("treatment", "gender"),
        correction = "BH",
        standardize = FALSE,
        cores = 1
      )
      data.frame(tissue = tissue_origin,
                 feature_level = level,
                 maaslin2_mod$results)
      
    }
  
  stopCluster(cl)
  end_time <- Sys.time()
  cat(tissue_origin, "MaAsLin2 calculated in : ", end_time - start_time, 
      attr(end_time - start_time, "units"), "\n")
  
  maaslin2_stats %<>% bind_rows(maaslin2_loop)
}
maaslin2_stats %>% glimpse()

maaslin2_stats %<>%
  decode_rfriendly_rows("feature") %>%
  select(-feature) %>%
  dplyr::rename(feature = fullnames)

write.xlsx(x = maaslin2_stats, 
           file = 'data/16S/differential-abundance/MaAsLin2/MaAsLin2_summary.xlsx', 
           overwrite = T)
saveRDS(maaslin2_stats, 
        file = 'data/16S/differential-abundance/MaAsLin2/MaAsLin2_summary.rds')
  

maaslin2_stats %>% 
  filter(metadata == "treatment") %>%
  mutate(significant = if_else(qval <= 0.01, "Yes", "No")) %>% 
  ggplot(aes(x=coef, y=-log10(qval))) +
  geom_point(aes(color = significant)) +
  facet_grid(tissue~feature_level, scales = "free") +
  scale_x_log10() +
  my_clean_theme() +
  scale_color_aaas()



#_______________________________________________________________________________
# 
# 
# corncob_stats <- list()
# 
# for (tissue_origin in c("Stool", "Jejunum", "Ileum")) {
#   ps_dat <- ps_obj %>%
#     subset_samples(tissue == tissue_origin)
#   
#   ps_dat_list <- list(
#     "ASVs" = ps_dat,
#     "Phylum" = aggregate_taxa(ps_dat, "Phylum"),
#     "Genus" = aggregate_taxa(ps_dat, "Genus"),
#     "Species" = aggregate_taxa(ps_dat, "Species")
#   )
#   
#   for (level in names(ps_dat_list)) {
#     dat_test <- ps_dat_list[[level]] %>%
#       microbiome::core(detection = 0, prevalence = 0.1)
#     sample_data(dat_test)$treatment <-
#       factor(sample_data(dat_test)$treatment,
#              levels = c("Saline", "Poly(I:C)"))
#     
#     corn_da <-
#       differentialTest(
#         formula = ~ gender + treatment,
#         phi.formula = ~ gender + treatment,
#         formula_null = ~ 1,
#         phi.formula_null = ~ gender + treatment,
#         data = dat_test,
#         test = "Wald",
#         boot = TRUE,
#         B = 100,
#         fdr_cutoff = 0.05
#       )
#     
#     corncob_stats[[tissue_origin]][[level]] <- corncob_stats
#     
#   }
# }
# 
# saveRDS(corncob_stats, file = "data/16S/differential-abundance/corncob_stats.rds")
# 
# 
# 
# #_______________________________________________________________________________
# 
# ancombc_stats <- list()
# 
# for (tissue_origin in c("Stool", "Jejunum", "Ileum")){
#   ps_dat <- ps_obj %>% 
#     subset_samples(tissue == tissue_origin) 
#   
#   ps_dat_list <- list(
#     "ASVs"= ps_dat,
#     "Phylum" = aggregate_taxa(ps_dat, "Phylum"),
#     "Genus" = aggregate_taxa(ps_dat, "Genus"),
#     "Species" = aggregate_taxa(ps_dat, "Species")
#   )
#   
#   #setup parallel processing
#   start_time <- Sys.time()
#   cores = detectCores()
#   cl <- makeCluster(cores[1] - 2) # to prevent computer overload
#   registerDoParallel(cl)
#   
#   ancombc_loop <- 
#     foreach(level = names(ps_dat_list), 
#             .combine = 'rbind', .verbose = F,
#             .packages = c("ANCOMBC", "phyloseq", "microbiome", 
#                           "magrittr", "tibble")) %dopar% {
#                             
#                             # filter NA values from metadata and abundance df
#                             dat_test <- ps_dat_list[[level]] %>% 
#                               microbiome::core(detection = 0, prevalence = 0.1)
#                             
#                             ancom_da <-
#                               ancombc(
#                                 phyloseq = dat_test,
#                                 formula = "gender + treatment",
#                                 p_adj_method = "fdr",
#                                 zero_cut = 0.90,
#                                 lib_cut = 1000,
#                                 group = "group",
#                                 struc_zero = TRUE,
#                                 neg_lb = F,
#                                 tol = 1e-5,
#                                 max_iter = 100,
#                                 conserve = TRUE,
#                                 alpha = 0.05,
#                                 global = F
#                               )
#                             
#                             data.frame(
#                               tissue = tissue_origin,
#                               feature_level = level,
#                               features = row.names(ancom_da$res$beta),
#                               beta = unlist(ancom_da$res$beta),
#                               se = unlist(ancom_da$res$se),
#                               W = unlist(ancom_da$res$W),
#                               p_val = unlist(ancom_da$res$p_val),
#                               q_val = unlist(ancom_da$res$q_val),
#                               diff_abn = unlist(ancom_da$res$diff_abn)) %>%
#                               rownames_to_column(var = "covariate")
#                             
#                           }
#   
#   stopCluster(cl)
#   end_time <- Sys.time()
#   cat(tissue_origin, level, "ANCOMBC calculated in : ", end_time - start_time, 
#       attr(end_time - start_time, "units"), "\n")
#   
#   ancombc_stats %<>% bind_rows(ancombc_loop)
# }
# 
# ancombc_stats %>% glimpse()
# 
# write.xlsx(x = ancombc_stats, file = 'data/16S/differential-abundance/ancombc_stats.xlsx')



