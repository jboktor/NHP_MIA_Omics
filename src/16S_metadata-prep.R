# Caltech - Mazmanian Lab
# Joe Boktor
# November 2021

source("src/_load_packages.R")
source("src/_misc_functions.R")

### Load metadata
sample_keys <- read_xlsx("input_files/sample_keys.xlsx", sheet = "keys")
seq_keys <- read_xlsx("input_files/sample_keys.xlsx", sheet = "16S-keys")


sample_keys %<>%
  select(subject_id, parent_sample_id, gender, sample_number, treatment, group, SampleID)

joined_meta <- seq_keys %>%
  left_join(sample_keys) %>%
  mutate_all(as.character) %>%
  dplyr::rename(monkey = SampleID)

id_row <- c("#q2:types", replicate(8, "categorical"))
names(id_row) <- joined_meta %>% colnames()

qiime_metadata_out <- id_row %>%
  bind_rows(joined_meta)

qiime_metadata_out %>% head()
qiime_metadata_out %>% glimpse()

# openxlsx::write.xlsx(qiime_metadata_out, file = 'input_files/NHP_16S_metadata_qiime_input.xlsx')
write.table(qiime_metadata_out,
  file = "input_files/NHP_16S_metadata_qiime_input.tsv",
  quote = FALSE, sep = "\t", col.names = T, row.names = F
)
