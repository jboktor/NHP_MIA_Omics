source("src/_load_packages.R")
source("src/_misc_functions.R")


for (tissue in tissue_order) {
  print(print_line)
  print(tissue)

  df_stats <-
    read_xlsx(paste0("input_files/NHP ", tissue, " Metabolon.xlsx"), sheet = "LogData")

  # annotate points of interest
  df_stats <-
    mutate(df_stats,
      male_feats = if_else(`Negative Log(PM/SM p-value)` >= 1.3 & `Log2 PM/SM` > 0, "A",
        if_else(`Negative Log(PM/SM p-value)` >= 1.3 &
          `Log2 PM/SM` < 0, "B", "C")
      )
    )
  df_stats <-
    mutate(df_stats,
      female_feats = if_else(`Negative Log(PF/SF p-value)` >= 1.3 & `Log2 PF/SF` > 0, "A",
        if_else(`Negative Log(PF/SF p-value)` >= 1.3 &
          `Log2 PF/SF` < 0, "B", "C")
      )
    )
  df_stats <-
    mutate(df_stats,
      all_sig = if_else(`Negative Log(P/S p-value)` >= 1.3 & `Log2 P/S` > 0, "A",
        if_else(`Negative Log(P/S p-value)` >= 1.3 &
          `Log2 P/S` < 0, "B", "C")
      )
    )

  male_features <- filter(df_stats, male_feats == "A" | male_feats == "B") %>%
    slice_max(order_by = `Negative Log(PM/SM p-value)`, n = 10)
  female_features <- filter(df_stats, female_feats == "A" | female_feats == "B") %>%
    slice_max(order_by = `Negative Log(PF/SF p-value)`, n = 10)
  all_features <- filter(df_stats, all_sig == "A" | all_sig == "B") %>%
    slice_max(order_by = `Negative Log(P/S p-value)`, n = 10)

  v1 <- volcano_plot(
    df = df_stats,
    x = df_stats$`Log2 PM/SM`,
    y = df_stats$`Negative Log(PM/SM p-value)`,
    xlabel = "[ Male Poly(I:C)/Control ]", colorvar = df_stats$male_feats
  ) +
    geom_label(
      df_labels = male_features,
      x = male_features$`Log2 PM/SM`,
      y = male_features$`Negative Log(PM/SM p-value)`,
      colorvar = male_features$male_feats
    )
  v1
  ggsave(v1,
    filename = paste0("figures/volcano_plots/", tissue, "_males.svg"),
    height = 5.5, width = 5
  )
  v2 <- volcano_plot(
    df = df_stats,
    x = df_stats$`Log2 PF/SF`,
    y = df_stats$`Negative Log(PF/SF p-value)`,
    xlabel = "[ Female Poly(I:C)/Control ]", colorvar = df_stats$female_feats
  ) +
    geom_label(
      df_labels = female_features,
      x = female_features$`Log2 PF/SF`,
      y = female_features$`Negative Log(PF/SF p-value)`,
      colorvar = female_features$female_feats
    )
  v2
  ggsave(v2,
    filename = paste0("figures/volcano_plots/", tissue, "_females.svg"),
    height = 5.5, width = 5
  )
  v3 <- volcano_plot(
    df = df_stats,
    x = df_stats$`Log2 P/S`,
    y = df_stats$`Negative Log(P/S p-value)`,
    xlabel = "[ Poly(I:C)/Control ]", colorvar = df_stats$all_sig
  ) +
    geom_label(
      df_labels = all_features,
      x = all_features$`Log2 P/S`,
      y = all_features$`Negative Log(P/S p-value)`,
      colorvar = all_features$all_sig
    )
  v3
  ggsave(v3,
    filename = paste0("figures/volcano_plots/", tissue, "_all.svg"),
    height = 5.5, width = 5
  )
}
