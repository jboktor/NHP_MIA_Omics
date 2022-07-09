
source("src/_load_packages.R")
source("src/_misc_functions.R")

# Behavioral data ----
df_behavior_postwean <-
  read_xlsx(
    'input_files/Behavior_Immune_data/LongStereo - Year four stereotyped behaviors.xlsx',
    sheet = 'PostWeanStereowithCodes'
  ) %>%
  select(-c(COHORT, FOCALDYE, GENDER, EXPCODE, EXPCODE2)) %>%
  pivot_longer(!c(Animal, Observation), names_to = "behavior") %>%
  mutate(behavior__timepoint = paste(behavior, Observation, sep = "__")) %>%
  right_join(meta_animal_id, by = "Animal") %>%
  right_join(meta_groups, by = "SampleID")

df_behavior_juvenile <-
  read_xlsx(
    'input_files/Behavior_Immune_data/LongStereo - Year four stereotyped behaviors.xlsx',
    sheet = 'JuvStereowithCodes'
  ) %>%
  select(-c(COHORT, FOCALDYE, GENDER, EXPCODE, EXPCODE2)) %>%
  pivot_longer(!c(Animal, Observation), names_to = "behavior") %>%
  mutate(behavior__timepoint = paste(behavior, Observation, sep = "__")) %>%
  right_join(meta_animal_id, by = "Animal") %>%
  right_join(meta_groups, by = "SampleID")

df_behavior_adult <-
  read_xlsx(
    'input_files/Behavior_Immune_data/LongStereo - Year four stereotyped behaviors.xlsx',
    sheet = 'SubAdultwithCodes'
  ) %>%
  select(-c(GENDER, EXPCODE, EXPCODE2)) %>%
  pivot_longer(!c(Animal, Observation), names_to = "behavior") %>%
  mutate(behavior__timepoint = paste(behavior, Observation, sep = "__")) %>%
  right_join(meta_animal_id, by = "Animal") %>%
  right_join(meta_groups, by = "SampleID")


# Data frame for exploration of behavior data
df_behav_eda <-
  bind_rows(df_behavior_postwean,
            df_behavior_juvenile,
            df_behavior_adult)

df_behavior <- df_behav_eda %>%
  select(SampleID, value, behavior__timepoint) %>%
  pivot_wider(names_from = "behavior__timepoint", values_from = "value") %>%
  arrange(match(SampleID, rownames(df_all))) %>%
  column_to_rownames(var = "SampleID")

behavior_data <- list("df_behavior" = df_behavior,
                      "df_behav_eda" = df_behav_eda)

saveRDS(behavior_data, file = "data/behaviors/2022-07-05_behaviors.rds")


