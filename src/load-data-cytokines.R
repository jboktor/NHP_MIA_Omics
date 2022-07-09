
source("src/_load_packages.R")
source("src/_misc_functions.R")

#' This script imputes cytokine profiles for downstream analyses

sample_meta <- read_xlsx("input_files/misc/sample_keys.xlsx", sheet = "keys")
meta_animal_id <- sample_meta %>%
  select(subject_id, SampleID) %>%
  dplyr::rename("Animal" = "subject_id")

# _______________________________________________________________________________
# Collate clean post-processed cytokine profiles for supplement ----

# GI Cytokine Profiles  --
df_GI_media <-
  read_xlsx("input_files/Behavior_Immune_data/AMA 27 GI tissues media 072513.xlsx",
    sheet = "Sheet1"
  ) %>%
  select(-c(Gender, Study_group, Description)) %>%
  pivot_longer(!c(Animal, `GI Section`), names_to = "cytokine") %>%
  separate(cytokine, c("cytokine", "trash"), sep = " ") %>%
  select(-trash) %>%
  right_join(meta_animal_id, by = "Animal") %>%
  dplyr::rename(Tissue = `GI Section`) %>%
  mutate(
    cytokine = toupper(cytokine),
    cytokine = gsub("\\-", "", cytokine),
    Tissue = glue("{Tissue}-PBMCs")
  ) %>%
  drop_na(cytokine) %>%
  pivot_wider(names_from = "cytokine", values_from = "value")

# Brain Cytokine Profiles
df_brain_media <-
  read_xlsx("input_files/Behavior_Immune_data/AMA27 - brain cytokines.xlsx",
    sheet = "Sheet1"
  ) %>%
  select(-c(Study, Study_group)) %>%
  pivot_longer(!Animal, names_to = "cytokine_longname") %>%
  mutate(cytokine_longname = gsub("_y4", "", cytokine_longname)) %>%
  separate(cytokine_longname, c("Tissue", "cytokine"), sep = "_") %>%
  mutate(Tissue = toupper(Tissue)) %>%
  mutate(cytokine = toupper(cytokine)) %>%
  pivot_wider(names_from = "cytokine", values_from = "value") %>%
  right_join(meta_animal_id, by = "Animal") %>%
  relocate(SampleID, .after = Tissue) %>%
  mutate(
    Tissue = gsub("PLASMA", "Plasma", Tissue),
    Tissue = gsub("CORTEX", "Cortex", Tissue)
  )

# Plasma Profiles
df_plasma_media <-
  read_xlsx("input_files/Behavior_Immune_data/Year 4 Blood Cytokine data.xlsx",
    sheet = "Media"
  ) %>%
  select(-c(Gender, treatment)) %>%
  pivot_longer(!Animal, names_to = "cytokine") %>%
  mutate(cytokine = gsub("m_y4_", "", cytokine)) %>%
  mutate(cytokine = toupper(cytokine)) %>%
  pivot_wider(names_from = "cytokine", values_from = "value") %>%
  right_join(meta_animal_id, by = "Animal") %>%
  mutate(Tissue = "Plasma-PBMCs") %>%
  relocate(Tissue, .after = Animal) %>%
  relocate(SampleID, .after = Tissue)

# Join all cytokine tables
all_cytokines <-
  bind_rows(
    df_GI_media,
    df_brain_media,
    df_plasma_media
  )

# Millipore PCYTMG-40k-PX23 Assay sensitivities
# (minimum detectable conc.) minDC (pg/mL) Average + 2SD

cytokine_lods <-
  tribble(
    ~cytokine, ~LOD_min,
    "IL4", 3.1,
    "IL10", 6.4,
    "IL18", 6.1,
    "SCD40L", 2.1,
    "TNFA", 1.6,
    "GMCSF", 1.8,
    "GCSF", 2.1,
    "IFNG", 1.6,
    "IL1RA", 2.4,
    "IL1B", 1.2,
    "IL2", 2.1,
    "IL5", 1.5,
    "IL6", 1.6,
    "IL8", 1.1,
    "IL12P40", 1.5,
    "IL13", 5.8,
    "IL15", 0.5,
    "IL17", 1.3,
    "MCP1", 3.1,
    "MIP1A", 4.9,
    "MIP1B", 1.6,
    "TGFA", 1.1,
    "VEGF", 13.6
  )

cytokine_lods %<>%
  mutate(impute_val = LOD_min / 2)

# Calculate the total number of measurements per tissue
# 1) in total and 2) requiring imputation

lod_stats <- tibble()

for (cyt in cytokine_lods$cytokine) {
  lod <- cytokine_lods %>%
    filter(cytokine == cyt) %>%
    pull(LOD_min)

  lod_stats2add <-
    all_cytokines %>%
    select(Tissue, all_of(cyt)) %>%
    drop_na(all_of(cyt)) %>%
    mutate(above_LOD = if_else(.data[[cyt]] >= lod, TRUE, FALSE)) %>%
    group_by(Tissue) %>%
    dplyr::summarize(
      total_samples = n(),
      samples_above_lod = sum(above_LOD),
      samples_below_lod = total_samples - samples_above_lod
    ) %>%
    mutate(cytokine = cyt)

  lod_stats %<>% bind_rows(lod_stats2add)
}

# Tissue cytokines with fewer than 5 measurements
blacklist <-
  lod_stats %>%
  filter(samples_above_lod < 5) %>%
  mutate(tissue_cyt = glue("{Tissue}_{cytokine}"))

all_cytokines_trimmed <-
  all_cytokines %>%
  pivot_longer(!c(Animal, Tissue, SampleID), names_to = "cytokine") %>%
  mutate(tissue_cyt = glue("{Tissue}_{cytokine}")) %>%
  filter(tissue_cyt %nin% blacklist$tissue_cyt) %>%
  select(-tissue_cyt)


# How many cytokines are measures in each tissue after removal of cytokines with low readouts?

cytokines_per_tissue <-
  all_cytokines_trimmed %>%
  select(Tissue, cytokine) %>%
  distinct() %>%
  group_by(Tissue) %>%
  dplyr::summarise(n = n())
cytokines_per_tissue

# A tibble: 10 × 2
# Tissue             n
# 1 AC                16
# 2 Colon PBMCs       23
# 3 Cortex            13
# 4 CSF                9
# 5 HIPP              19
# 6 Ileum PBMCs       23
# 7 Jejunum PBMCs     23
# 8 MLN PBMCs         12
# 9 Plasma            22
# 10 Plasma - PBMCs   23

# How many times is each cytokine measured in a tissue successfully?

cytokine_prevalence <-
  all_cytokines_trimmed %>%
  select(Tissue, cytokine) %>%
  distinct() %>%
  group_by(cytokine) %>%
  dplyr::summarise(n = n())

cytokine_prevalence %>% print(n = Inf)

# A tibble: 23 × 2
# cytokine     n
# 1 GCSF        10
# 2 GMCSF        7
# 3 IFNG         6
# 4 IL10         8
# 5 IL12P40      6
# 6 IL13         8
# 7 IL15         9
# 8 IL17         6
# 9 IL18         8
# 10 IL1B         8
# 11 IL1RA        9
# 12 IL2          8
# 13 IL4          5
# 14 IL5          8
# 15 IL6          9
# 16 IL8          9
# 17 MCP1         8
# 18 MIP1A       10
# 19 MIP1B        8
# 20 SCD40L       9
# 21 TGFA         9
# 22 TNFA         6
# 23 VEGF         9

# _______________________________________________________________________________
# Imputation ----
# replace cytokine measurements below the LOD with 1/2 LOD

all_cytokines_trimmed_imputed <- all_cytokines_trimmed

for (cyt in cytokine_lods$cytokine) {
  print(cyt)
  cyt_thres <- cytokine_lods %>% filter(cytokine == cyt)
  lod <- cyt_thres$LOD_min
  impute <- cyt_thres$impute_val

  all_cytokines_trimmed_imputed %<>%
    mutate(value = case_when(
      (cytokine == cyt & value < lod & !is.na(value)) ~ impute,
      TRUE ~ value
    ))
}


# Save imputation summary statistics ----

imputation_summary_stats <- list(
  "Measurement Quality" = lod_stats,
  "cytokines_per_tissue" = cytokines_per_tissue,
  "cytokine_prevalence" = cytokine_prevalence
)

write.xlsx(imputation_summary_stats,
  "data/cytokines/2022-07-05_cytokine-summary-qc.xlsx",
  overwrite = T
)

# Save filtered and imputed data ----
cytokine_profiles <- list(
  "imputed" = all_cytokines_trimmed_imputed,
  "raw" = all_cytokines %>%
    pivot_longer(!c(Animal, Tissue, SampleID), names_to = "cytokine")
)

write.xlsx(
  cytokine_profiles,
  "data/cytokines/2022-07-05_cytokine-profiles.xlsx",
  overwrite = T
)

saveRDS(all_cytokines_trimmed_imputed,
  file = "data/cytokines/2022-07-05_cytokine-profiles_imputed.rds"
)
