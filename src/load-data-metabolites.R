source("src/_load_packages.R")
source("src/_misc_functions.R")

# Create one dataframe with all scaled / imputed data
df_plasma_scaled_imp <-
  read_xlsx("input_files/Metabolomics/NHP Plasma HD4 Metabolon.xlsx", sheet = "ScaledImpData") %>%
  mutate(tissue = "Plasma")
df_csf_scaled_imp <-
  read_xlsx("input_files/Metabolomics/NHP CSF Metabolon.xlsx", sheet = "ScaledImpData") %>%
  mutate(tissue = "CSF")
df_jejunum_scaled_imp <-
  read_xlsx("input_files/Metabolomics/NHP Jejunum Metabolon.xlsx", sheet = "ScaledImpData") %>%
  mutate(tissue = "Jejunum")
df_ileum_scaled_imp <-
  read_xlsx("input_files/Metabolomics/NHP Ileum Metabolon.xlsx", sheet = "ScaledImpData") %>%
  mutate(tissue = "Ileum")
df_colon_scaled_imp <-
  read_xlsx("input_files/Metabolomics/NHP Colon Metabolon.xlsx", sheet = "ScaledImpData") %>%
  mutate(tissue = "Colon")
df_feces_scaled_imp <-
  read_xlsx("input_files/Metabolomics/NHP Feces Metabolon.xlsx", sheet = "ScaledImpData") %>%
  mutate(tissue = "Feces")
keys <- read_xlsx("input_files/misc/sample_keys.xlsx", sheet = "keys") %>%
  select(SampleID, group, treatment)
# keys2 <- read_xlsx('input_files/misc/sample_keys.xlsx', sheet = 'keys') %>%
#   select(subject_id, SampleID) %>%
#   dplyr::rename("Animal" = "subject_id" )


df_scaled_imp <- bind_rows(
  df_plasma_scaled_imp,
  df_csf_scaled_imp,
  df_jejunum_scaled_imp,
  df_ileum_scaled_imp,
  df_colon_scaled_imp,
  df_feces_scaled_imp
)

df_all <- df_scaled_imp %>%
  select(-c(`SUPER PATHWAY`, `SUB PATHWAY`)) %>%
  pivot_longer(!c(BIOCHEMICAL, tissue), names_to = "SampleID") %>%
  left_join(keys, by = "SampleID") %>%
  mutate(analyte__tissue = paste(BIOCHEMICAL, tissue, sep = "__")) %>%
  select(SampleID, value, analyte__tissue) %>%
  pivot_wider(values_from = "value", names_from = "analyte__tissue") %>%
  column_to_rownames(var = "SampleID")

pathway_map <- df_scaled_imp %>%
  select(BIOCHEMICAL, `SUPER PATHWAY`, `SUB PATHWAY`) %>%
  distinct()

saveRDS(df_scaled_imp, file = "data/metabolomics/2022-07-05_scaled_imputed_abundance_annotated.rds")
saveRDS(df_all, file = "data/metabolomics/2022-07-05_scaled_imputed_abundance.rds")
saveRDS(pathway_map, file = "data/metabolomics/2022-07-05_pathway-map.rds")
