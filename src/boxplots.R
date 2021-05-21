
source("src/_load_packages.R")
source("src/_misc_functions.R")


for (tissue in tissue_order) {
  
  print(print_line); print(tissue); print(print_line) 
  
  df_stats <-
    read_xlsx(paste0('input_files/NHP ', tissue, ' Metabolon.xlsx'), sheet = 'LogData') %>%
    janitor::clean_names()
  df_scaled_imp <-
    read_xlsx(paste0('input_files/NHP ', tissue, ' Metabolon.xlsx'), sheet = 'ScaledImpData_t')
  keys <-
    read_xlsx('input_files/sample_keys.xlsx', sheet = 'keys') %>%
    select(SampleID, group, treatment)
  
  # extract names of metabolites significant for treatment
  treatment_sig <- df_stats %>%
    dplyr::filter(all_poly_ic_all_saline_p_value <= 0.05) %>%
    select(biochemical_name) %>%
    unlist(use.names = F)
  
  # extract names of metabolites significant for gender specific treatment
  gender_treatment_sig <- df_stats %>%
    dplyr::filter((
      poly_ic_male_saline_male_p_value <= 0.05 |
        poly_ic_female_saline_female_p_value <= 0.05
    ) &
      all_poly_ic_all_saline_p_value > 0.05
    ) %>%
    select(biochemical_name) %>%
    unlist(use.names = F)
  
  # transform scaled/imputed data into long format and join with metadata
  df_scaled_imp <- df_scaled_imp %>%
    left_join(keys, by = "SampleID") %>%
    pivot_longer(!c(SampleID, group, treatment), names_to = "metabolite")
  
  boxplot_loop_treatment(
    significant_mets = treatment_sig,
    scaled_df = df_scaled_imp,
    stats_df = df_stats,
    tissue = tissue
  )
  
  boxplot_loop_treatment_sex(
    significant_mets = gender_treatment_sig,
    scaled_df = df_scaled_imp,
    stats_df = df_stats,
    tissue = tissue
  )
  
}

