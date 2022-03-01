
source("src/_load_packages.R")
source("src/_misc_functions.R")

# Create one dataframe with all scaled / imputed data
df_plasma_stats <-
  read_xlsx('input_files/NHP Plasma HD4 Metabolon.xlsx', sheet = 'LogData') %>% 
  mutate(tissue = "Plasma")
df_csf_stats <-
  read_xlsx('input_files/NHP CSF Metabolon.xlsx', sheet = 'LogData') %>% 
  mutate(tissue = "CSF")
df_jejunum_stats <-
  read_xlsx('input_files/NHP Jejunum Metabolon.xlsx', sheet = 'LogData') %>% 
  mutate(tissue = "Jejunum")
df_ileum_stats <-
  read_xlsx('input_files/NHP Ileum Metabolon.xlsx', sheet = 'LogData') %>% 
  mutate(tissue = "Ileum")
df_colon_stats <-
  read_xlsx('input_files/NHP Colon Metabolon.xlsx', sheet = 'LogData') %>% 
  mutate(tissue = "Colon")
df_feces_stats <-
  read_xlsx('input_files/NHP Feces Metabolon.xlsx', sheet = 'LogData') %>% 
  mutate(tissue = "Feces")

df_stats <- 
  bind_rows(df_plasma_stats, 
            df_csf_stats,
            df_jejunum_stats,
            df_ileum_stats,
            df_colon_stats,
            df_feces_stats) %>% 
  janitor::clean_names()

library(gtsummary)
# make dataset with a few variables to summarize
df_stats %>% 
  filter(poly_ic_male_saline_male_p_value <= 0.05) %>% 
  select(tissue) %>% 
  tbl_summary()
df_stats %>% 
  filter(poly_ic_female_saline_female_p_value <= 0.05) %>% 
  select(tissue) %>% 
  tbl_summary()
df_stats %>% 
  filter(all_poly_ic_all_saline_p_value <= 0.05) %>% 
  select(tissue) %>% 
  tbl_summary()


df_heatmap_detected <-
  df_stats %>%
  select(biochemical_name, tissue) %>%
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1") %>% 
  pivot_longer(!biochemical_name, names_to = "tissue", values_to = "detected") %>% 
  mutate(variable = "detected")

df_heatmap_treatment <-
  df_stats %>%
  filter(all_poly_ic_all_saline_p_value < 0.05) %>% 
  select(biochemical_name, tissue) %>%
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1") %>% 
  pivot_longer(!biochemical_name, names_to = "tissue", values_to = "detected") %>% 
  mutate(variable = "Treatment Effect")

df_heatmap_male_treatment <-
  df_stats %>%
  filter(poly_ic_male_saline_male_p_value < 0.05) %>% 
  select(biochemical_name, tissue) %>%
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1") %>% 
  pivot_longer(!biochemical_name, names_to = "tissue", values_to = "detected") %>% 
  mutate(variable = "Male Treatment Effect")

df_heatmap_female_treatment <-
  df_stats %>%
  filter(poly_ic_female_saline_female_p_value < 0.05) %>% 
  select(biochemical_name, tissue) %>%
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1") %>% 
  pivot_longer(!biochemical_name, names_to = "tissue", values_to = "detected") %>% 
  mutate(variable = "Female Treatment Effect")

df_heatmap_sex <-
  df_stats %>%
  filter(poly_ic_female_saline_female < 0.05 | 
           poly_ic_male_saline_male < 0.05) %>% 
  select(biochemical_name, tissue) %>%
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1") %>% 
  pivot_longer(!biochemical_name, names_to = "tissue", values_to = "detected") %>% 
  mutate(variable = "Sex Effect")

df.heatmap <- 
  bind_rows(df_heatmap_detected,
            df_heatmap_treatment, 
            df_heatmap_male_treatment,
            df_heatmap_female_treatment,
            df_heatmap_sex) %>% 
  mutate(tissue = factor(tissue, levels = tissue_order))

#_______________________________________________________________________________
#                           Association heatmap ---- 
#_______________________________________________________________________________


shared.sig.heatmap <- 
  df.heatmap %>% 
  mutate(variable = factor(variable, levels = 
                             c("Treatment Effect", "Male Treatment Effect", 
                               "Female Treatment Effect", "Sex Effect"))) %>%
  ggplot(aes(x=tissue, y=biochemical_name, fill=detected)) +
  geom_tile() +
  theme_classic() +
  labs(x = NULL, y = "Metabolites (P < 0.05)") +
  scale_fill_manual(values = c("0" = "#f1f1f1", "1" = "#434343")) +
  facet_wrap(~variable, nrow = 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  )
shared.sig.heatmap
# ggsave(shared.sig.heatmap, filename = "figures/summary/Shared_Significance_Heatmap_facet.svg",
#        width = 7, height = 9)



#_______________________________________________________________________________
#                                 Upset plots ---- 
#_______________________________________________________________________________


nhp_meta <- data.frame(sets = colnames(df.upset.female.treatment),
                       tissue_group = c("Plasma", "CSF"," GI",
                                        "GI", "GI","GI"))


df.upset.treatment <-
  df.heatmap %>%
  filter(variable == "Treatment Effect") %>%
  pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
  select(-c(biochemical_name, variable)) %>% 
  mutate_if(is.character, as.numeric) %>% 
  as.data.frame()

upset.treatment <-
  upset(
    df.upset.treatment,
    nintersects = 30,
    nsets = 30,
    order.by = "freq",
    decreasing = T,
    mb.ratio = c(0.6, 0.4),
    number.angles = 0,
    text.scale = 1.1,
    point.size = 2.8,
    line.size = 1
  )

svg(file="figures/summary/treatment_upset.svg", width = 7, height = 4)
upset.treatment
dev.off()


df.upset.male.treatment <-
  df.heatmap %>%
  filter(variable == "Male Treatment Effect") %>%
  pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
  select(-c(biochemical_name, variable)) %>% 
  mutate_if(is.character, as.numeric) %>% 
  as.data.frame()

upset.male.treatment <-
  upset(
    df.upset.male.treatment,
    nintersects = 30,
    nsets = 30,
    order.by = "freq",
    decreasing = T,
    mb.ratio = c(0.6, 0.4),
    number.angles = 0,
    text.scale = 1.1,
    point.size = 2.8,
    line.size = 1
  )

svg(file="figures/summary/treatment_male_upset.svg", width = 7, height = 4)
upset.male.treatment
dev.off()


df.upset.female.treatment <-
  df.heatmap %>%
  filter(variable == "Female Treatment Effect") %>%
  pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
  select(-c(biochemical_name, variable)) %>% 
  mutate_if(is.character, as.numeric) %>% 
  as.data.frame()


upset.female.treatment <-
  upset(
    df.upset.female.treatment,
    nintersects = 30,
    nsets = 30,
    order.by = "freq",
    decreasing = T,
    mb.ratio = c(0.6, 0.4),
    number.angles = 0,
    text.scale = 1.1,
    point.size = 2.8,
    line.size = 1 #, 
    # set.metadata = 
    #   list(data = nhp_meta, plots = 
    #          list(
    #            list(
    #              type = "matrix_rows",
    #              column = "tissue_group",
    #              colors = c(Plasma = "green",
    #                         CSF = "navy", GI = "purple"),
    #              alpha = 0.5
    #            )))
    )

svg(file="figures/summary/treatment_female_upset.svg", width = 7, height = 4)
upset.female.treatment
dev.off()




df.upset.detected <-
  df.heatmap %>%
  filter(variable == "detected") %>%
  pivot_wider(names_from = 'tissue', values_from = 'detected') %>%
  select(-c(biochemical_name, variable)) %>% 
  mutate_if(is.character, as.numeric) %>% 
  as.data.frame()


upset.detected <-
  upset(
    df.upset.detected,
    nintersects = 30,
    nsets = 30,
    order.by = "freq",
    decreasing = T,
    mb.ratio = c(0.6, 0.4),
    number.angles = 0,
    text.scale = 1.1,
    point.size = 2.8,
    line.size = 1
  )

svg(file="figures/summary/all_detected_upset.svg", width = 7, height = 4)
upset.detected
dev.off()



df_heatmap_detected_output <-
  df_stats %>%
  select(biochemical_name, tissue) %>%
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1") 

df_heatmap_treatment_output <-
  df_stats %>%
  filter(all_poly_ic_all_saline_p_value < 0.05) %>% 
  select(biochemical_name, tissue) %>%
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1")

df_heatmap_male_treatment_output <-
  df_stats %>%
  filter(poly_ic_male_saline_male_p_value < 0.05) %>% 
  select(biochemical_name, tissue) %>%
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1")

df_heatmap_female_treatment_output <-
  df_stats %>%
  filter(poly_ic_female_saline_female_p_value < 0.05) %>% 
  select(biochemical_name, tissue) %>%
  pivot_wider(names_from = "tissue", 
              values_from = "tissue", 
              values_fill = "0",
              values_fn = function(x) "1")

# write.xlsx(df_heatmap_detected_output,
#            'data/summary_info/tissue_shared_metabolites_all_detected.xlsx')
# write.xlsx(df_heatmap_treatment_output,
#            'data/summary_info/tissue_shared_metabolites_polyic.vs.saline_significant.xlsx')
# write.xlsx(df_heatmap_male_treatment_output,
#            'data/summary_info/tissue_shared_metabolites_MALE_polyic.vs.saline_significant.xlsx')
# write.xlsx(df_heatmap_female_treatment_output,
#            'data/summary_info/tissue_shared_metabolites_FEMALE_polyic.vs.saline_significant.xlsx')





# upset(
#   movies,
#   set.metadata = list(data = metadata, plots = list(
#     list(
#       type= "hist",
#       column = "avgRottenTomatoesScore",
#       assign = 20
#     ),
#     list(
#       type = "bool",
#       column = "accepted",
#       assign = 5,
#       colors = c("#FF3333", "#006400")
#     ),
#     list(
#       type = "text",
#       column = "Cities",
#       assign = 5,
#       colors = c(Boston = "green", NYC = "navy", LA = "purple")
#     ),
#     list(
#       type = "matrix_rows",
#       column = "Cities",
#       colors = c(Boston = "green",
#                  NYC = "navy", LA = "purple"),
#       alpha = 0.5
#     )
#   )),
#   queries = list(
#     list(
#       query = intersects,
#       params = list("Drama"),
#       color = "red",
#       active = F
#     ),
#     list(
#       query = intersects,
#       params = list("Action", "Drama"),
#       active = T
#     ),
#     list(
#       query = intersects,
#       params = list("Drama", "Comedy", "Action"),
#       color = "orange",
#       active = T
#     )
#   ),
#   attribute.plots = list(
#     gridrows = 45,
#     plots = list(
#       list(
#         plot = scatter_plot,
#         x = "ReleaseDate",
#         y = "AvgRating",
#         queries = T
#       ),
#       list(
#         plot = scatter_plot,
#         x = "AvgRating",
#         y = "Watches",
#         queries = F
#       )
#     ),
#     ncols = 2
#   ),
#   query.legend = "bottom"
# )




