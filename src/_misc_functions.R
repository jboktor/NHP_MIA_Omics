# Miscellaneous functions

#_______________________________________________________________________________
#                          Aesthetic variables ----
#_______________________________________________________________________________

tissue_colors = c(
  "Colon" = "#E41A1C",
  "CSF" = "#377EB8",
  "Feces" = "#4DAF4A",
  "Ileum" = "#984EA3",
  "Jejunum"  = "#FF7F00",
  "Plasma" = "#A65628"
)

tissue_order <-
  c("CSF",
    "Plasma",
    "Jejunum",
    "Ileum",
    "Colon",
    "Feces")

SuperpathwayOrder <-
  c("Amino Acid",
    "Lipid",
    "Carbohydrate",
    "Nucleotide",
    "Cofactors and Vitamins",
    "Peptide",
    "Xenobiotics",
    "Energy",
    "Partially Characterized Molecules"
  )

super_pathways_alphetical <-
  c("Amino Acid",
    "Carbohydrate",
    "Cofactors and Vitamins",
    "Energy",
    "Lipid",
    "Nucleotide",
    "Partially Characterized Molecules",
    "Peptide",
    "Xenobiotics"
  )

super_pathway_cols<- brewer.pal(9, "RdYlBu")
names(super_pathway_cols) <- super_pathways_alphetical
super_pathway_cols2 <- pal_futurama()(9)
names(super_pathway_cols2) <- super_pathways_alphetical

# Plotting Labels & Colors
group_ord <- c("Poly(I:C) Male",
               "Saline Male",
               "Poly(I:C) Female",
               "Saline Female")
  
group_cols <-
  c(
    "Poly(I:C) Male" = "#1f78b4",
    "Saline Male" = "#a6cee3",
    "Poly(I:C) Female" = "#e31a1c",
    "Saline Female" = "#fb9a99"
  )
cols_gt <-
  c(
    "Poly(I:C) Male" = "#011c40",
    "Saline Male" = "#836d44",
    "Poly(I:C) Female" = "#3c6db4",
    "Saline Female" = "#c3ac84"
  )
cols_gtf <-
  c(
    "Poly(I:C) Male" = "#3c6db4",
    "Saline Male" = "#c3ac84",
    "Poly(I:C) Female" = "#3c6db4",
    "Saline Female" = "#c3ac84"
  )
cols_t <- c("Saline" = "#836d44", "Poly(I:C)" = "#011c40")
cols_tf <- c("Saline" = "#c3ac84", "Poly(I:C)" = "#3c6db4")
comparisons.Treatment = list(c("Saline", "Poly(I:C)"))
comparisons.Treatment.Male = list(c("Poly(I:C) Male", "Saline Male"))
comparisons.Treatment.Female = list(c("Poly(I:C) Female", "Saline Female"))

colormash <- c(
  # Set 1
  '#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
  '#ff7f00', '#ffff33', '#a65628', '#f781bf','#999999',
  # Accent
  '#7fc97f','#beaed4','#fdc086','#ffff99',
  '#386cb0','#f0027f','#bf5b17','#666666',
  # Dark
  '#1b9e77','#d95f02','#7570b3','#e7298a',
  '#66a61e','#e6ab02','#a6761d','#666666')

# _______________________________________________________________________________

my_clean_theme <- function() {
  th <- ggplot2::theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5)
    )
  return(th)
}
# _______________________________________________________________________________

sig_mapper <- function(pval, shh = F, porq = "p", symbols = T) {
  ###' Traditional mapping of p-value to symbol 
  ###' prints p-values if below significance
  
  if (symbols == T){
    if (is.na(pval)){
      sigvalue = ""
    } else if (pval <= .001) {
      sigvalue = "***"
    } else if (pval <= .01) {
      sigvalue = "**"
    } else if (pval <= .05) {
      sigvalue = "*"
    } else if (pval > .05 & shh == F) {
      sigvalue = paste0(porq, "=", format.pval(pval, digits=2)) 
    } else if (pval > .05 & shh == T) {
      sigvalue = ""
    }
  } else if (symbols == F){
    sigvalue = paste0(porq, "=", format.pval(pval, digits=2)) 
  }
  return(sigvalue)
}

# _______________________________________________________________________________

# Treatment specific boxplots

boxplot_loop_treatment <- function(significant_mets,
                                   scaled_df,
                                   stats_df,
                                   tissue){
  
  for (analyte in significant_mets) {
    
    sigvalue <- stats_df %>% 
      filter(biochemical_name == analyte) %>% 
      select(all_poly_ic_all_saline_p_value) %>% 
      sig_mapper() 
    
    plot <-
      scaled_df %>% 
      filter(metabolite == analyte) %>% 
      ggplot(aes(x=treatment, y=value)) +
      geom_boxplot(aes(colour=treatment), outlier.alpha = 0, width = 0.4) +
      geom_point(aes(fill=treatment), position = position_jitterdodge(),shape=21, size=2) +
      labs(y = "Scaled Intensity", title = analyte) +
      theme_classic() +
      scale_color_manual(values = cols_t) +
      scale_fill_manual(values = cols_tf) +
      geom_signif(comparisons = comparisons.Treatment, 
                  annotations=sigvalue, tip_length = 0) +
      theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),
            axis.title.x=element_blank(),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5))
    print(plot)
    analyte_label <- gsub('/', '_', analyte)
    output_label <- paste0("figures/boxplots/", tissue, "/", tissue, 
                           "_treatment_significant_", analyte_label,".svg")
    cat(output_label, "\n")
    ggsave(plot, filename = output_label, height = 5, width =3.5)
  }
  
}

# _______________________________________________________________________________

# Treatment gender boxplots


boxplot_loop_treatment_sex <- function(significant_mets,
                                   scaled_df,
                                   stats_df,
                                   tissue){

  
  for (analyte in significant_mets) {
    
    Msigvalue <- stats_df %>% 
      filter(biochemical_name == analyte) %>% 
      select(poly_ic_male_saline_male_p_value) %>% 
      sig_mapper()
    Fsigvalue <- stats_df %>% 
      filter(biochemical_name == analyte) %>% 
      select(poly_ic_female_saline_female_p_value) %>% 
      sig_mapper()
    
    plot <-
      scaled_df %>% 
      filter(metabolite == analyte) %>% 
      ggplot(aes(x=fct_relevel(group, group_ord), y=value)) +
      geom_boxplot(aes(colour=group), outlier.alpha = 0, width = 0.4) +
      geom_point(aes(fill=group), position = position_jitterdodge(),shape=21, size=2) +
      labs(y = "Scaled Intensity", title = analyte) +
      theme_classic() +
      scale_fill_manual(values = cols_gtf) +
      scale_color_manual(values = cols_gt) +
      geom_signif(comparisons = comparisons.Treatment.Male, 
                  annotations=Msigvalue, tip_length = 0) +
      geom_signif(comparisons = comparisons.Treatment.Female, 
                  annotations=Fsigvalue, tip_length = 0) +
      theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),
            axis.title.x=element_blank(),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5))
    print(plot)
    analyte_label <- gsub('/', '_', analyte)
    output_label <- paste0("figures/boxplots/", tissue, "/", tissue, 
                           "_treatment_gender_significant_", analyte_label,".svg")
    cat(output_label, "\n")
    ggsave(plot, filename = output_label, height = 5, width =4)
  }
  
}

# _______________________________________________________________________________

print_line <- paste(replicate(20, "---"), collapse = "")


# _______________________________________________________________________________

super_path_avg <- function(df, keys) {
  
  df_pathway <- df %>% 
    select(-c(BIOCHEMICAL, `SUB PATHWAY`)) %>% 
    pivot_longer(!c(`SUPER PATHWAY`, tissue), names_to = "SampleID") %>% 
    group_by(`SUPER PATHWAY`, SampleID, tissue) %>% 
    dplyr::summarize(path_mean = mean(value, na.rm = TRUE)) %>% 
    left_join(keys, by = "SampleID")
  
  return(df_pathway)
}

sub_path_avg <- function(df, keys) {
  
  df_pathway <- df %>% 
    select(-c(BIOCHEMICAL, `SUPER PATHWAY`)) %>% 
    pivot_longer(!c(`SUB PATHWAY`, tissue), names_to = "SampleID") %>% 
    group_by(`SUB PATHWAY`, SampleID, tissue) %>% 
    dplyr::summarize(path_mean = mean(value, na.rm = TRUE)) %>% 
    left_join(keys, by = "SampleID")
  
  return(df_pathway)
}


#_______________________________________________________________________________
#                          Correlation Functions ----
#_______________________________________________________________________________

tidy_cor_output <- function(rcorr_out){
  df_rho <-
    rcorr_out$r %>% as.data.frame() %>%
    rownames_to_column(var = "feature_A") %>%
    pivot_longer(!feature_A, names_to = "feature_B", values_to = "rho")
  df_n <-
    rcorr_out$n %>% as.data.frame() %>%
    rownames_to_column(var = "feature_A") %>%
    pivot_longer(!feature_A, names_to = "feature_B", values_to = "n")
  df_p <-
    rcorr_out$P %>% as.data.frame() %>%
    rownames_to_column(var = "feature_A") %>%
    pivot_longer(!feature_A, names_to = "feature_B", values_to = "P")
  
  df_out <- df_rho %>% 
    left_join(df_n, by = c("feature_A", "feature_B")) %>% 
    left_join(df_p, by = c("feature_A", "feature_B"))
  
  return(df_out)
  
}

# _______________________________________________________________________________


corr_loop <- function(dfA, dfB, obj.name) {

  # Set dfA as behavior or cytokine data and dfB as metabolic profile

  corr_output <- tibble()
  for (metavar in colnames(dfA)) {
    cat("Calculating correlations for: ", metavar[[1]], "\n")
    # for (feature in colnames(dfB))
    for(feature in colnames(dfB)) {
      # Calculate Spearman's Correlation
      spearman <-
        cor.test(
          x = dfA[[metavar]],
          y = dfB[[feature]],
          method = "spearman",
          na.action = na.exclude,
          alternative = "two.sided"
        )
      row2add <-
        cbind(
          "feature_A" = metavar,
          "feature_B" = feature,
          "object_name" = obj.name,
          "rho" = spearman$estimate[[1]],
          "S" = spearman$statistic[[1]],
          "n" = length(na.omit(dfA[[metavar]])),
          "p" = spearman$p.value[[1]]
        )
      corr_output <- rbind(corr_output, row2add)
    }
  }
  statvars <- c("rho", "S", "n", "p")
  corr_output <-
    corr_output %>%
    na.omit() %>%
    mutate(across(all_of(statvars), as.character),
           across(all_of(statvars), as.numeric))
  return(corr_output)
}


corr_loop_parallel <- function(dfA, dfB, obj.name) {

  # Set dfA as behavior or cytokine data and dfB as metabolic profile
  require(foreach)
  require(doParallel)

  #setup parallel backend to use many processors
  cores = detectCores()
  cl <- makeCluster(cores[1] - 1) 
  registerDoParallel(cl)

  corr_output <- tibble()
  loop <-
    foreach(metavar = colnames(dfA), .combine='rbind') %:%
    foreach(feature = colnames(dfB), .combine='rbind') %dopar% {

    # Calculate Spearman's Correlation
    spearman <-
      cor.test(
        x = dfA[[metavar]],
        y = dfB[[feature]],
        method = "spearman",
        na.action = na.exclude,
        alternative = "two.sided"
      )
    data.frame(
      "feature_A" = metavar,
      "feature_B" = feature,
      "object_name" = obj.name,
      "rho" = spearman$estimate[[1]],
      "S" = spearman$statistic[[1]],
      "n" = length(na.omit(dfA[[metavar]])),
      "p" = spearman$p.value[[1]]
    )
  }
  stopCluster(cl)
  return(loop)
}

#_______________________________________________________________________________


corr_heatmap <- function(corr.df){
  
  # Calculate number of significant (P-VALUE < 0.05) associations
  #  and select top 30 features
  corr.df.top <- corr.df %>% 
    dplyr::filter(p < 0.05) %>%
    dplyr::group_by(feature_B) %>%
    dplyr::summarise(n = n()) %>%
    slice_max(n = 30, order_by = n, with_ties = F)
  
  # Filter correlation df for top 30 features
  corr.df.trim <- corr.df %>% 
    mutate(feature_A = as.character(feature_A)) %>% 
    filter(feature_B %in% corr.df.top$feature_B)
  
  # create dataframe matrix of Rho correlation values for distance functions
  rho.df <- 
    corr.df.trim %>% 
    dplyr::select(feature_B, feature_A, rho) %>% 
    pivot_wider(names_from = feature_B, values_from = rho, values_fill = NA) %>% 
    column_to_rownames(var = "feature_A")
  
  # Hierarchical clustering of Rho values for features & feature_A
  feature_A.dendro <- 
    as.dendrogram(hclust(d = dist(x = rho.df), method = "complete"))
  feature_A.dendro.plot <- ggdendrogram(data = feature_A.dendro, rotate = TRUE)
  feature_B.dendro <- 
    as.dendrogram(hclust(d = dist(x = as.data.frame(t(rho.df))), method = "complete"))
  feature_B.dendro.plot <- ggdendrogram(data = feature_B.dendro, rotate = TRUE)
  
  ### Reorder Heatmap axis using order of dendrograms
  feature_B.order <- order.dendrogram(feature_B.dendro)
  feature_A.order <- order.dendrogram(feature_A.dendro)
  corr.df.trim.ordered <- corr.df.trim %>% 
    dplyr::mutate(feature = factor(feature_B, 
                                   ordered = TRUE,
                                   levels = unique(corr.df.trim$feature_B)[feature_B.order] )) %>% 
    dplyr::mutate(feature_A = factor(feature_A, 
                                    ordered = TRUE,
                                    levels = unique(corr.df.trim$feature_A)[feature_A.order] )) 
  
  ### Plot Heatmap
  h1 <- 
    corr.df.trim.ordered %>% 
    mutate(siglabel = if_else(q < 0.20, "*", "")) %>%
    ggplot(aes(x = feature_A_obj, 
               y = fct_relevel(feature_B, unique(corr.df.trim$feature_B)[feature_B.order]),
               fill = rho)) +
    geom_tile() + 
    geom_text(aes(label=siglabel), size=5, vjust = 0.77, color = "white") +
    labs(fill = "Spearman's correlation") +
    scale_y_discrete(position = "right") +
    scale_fill_distiller(palette = "RdBu") +
    theme(axis.text.x = element_text(angle=45, hjust =1),
          legend.position = "top",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(1, 1, 1, 2), "cm"))
  
  print(h1)
  return(h1)
}


#_______________________________________________________________________________

corr_heatmap_facet <- function(corr.df, 
                               facet_ord = c("Jejunum", "Ileum", "MLN", "Colon")){
  
  # Calculate number of significant (P-VALUE < 0.05) associations
  #  and select top 30 features
  corr.df.top <- corr.df %>% 
    dplyr::filter(p < 0.05) %>%
    dplyr::group_by(feature_B) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n)) %>%
    top_n(n = 30, wt = n) %>% 
    slice_head(n = 30)
  
  # Filter correlation df for top 30 features
  corr.df.trim <- corr.df %>% 
    mutate(feature_A = as.character(feature_A)) %>% 
    filter(feature_B %in% corr.df.top$feature_B)
  
  # create dataframe matrix of Rho correlation values for distance functions
  rho.df <- 
    corr.df.trim %>% 
    dplyr::select(feature_B, feature_A, rho) %>% 
    pivot_wider(names_from = feature_B, values_from = rho, values_fill = NA) %>% 
    column_to_rownames(var = "feature_A")
  
  # Hierarchical clustering of Rho values for features & feature_A
  feature_A.dendro <- 
    as.dendrogram(hclust(d = dist(x = rho.df), method = "complete"))
  feature_A.dendro.plot <- ggdendrogram(data = feature_A.dendro, rotate = TRUE)
  feature_B.dendro <- 
    as.dendrogram(hclust(d = dist(x = as.data.frame(t(rho.df))), method = "complete"))
  feature_B.dendro.plot <- ggdendrogram(data = feature_B.dendro, rotate = TRUE)
  
  ### Reorder Heatmap axis using order of dendrograms
  feature_B.order <- order.dendrogram(feature_B.dendro)
  feature_A.order <- order.dendrogram(feature_A.dendro)
  corr.df.trim.ordered <- corr.df.trim %>% 
    dplyr::mutate(feature = factor(feature_B, 
                                   ordered = TRUE,
                                   levels = unique(corr.df.trim$feature_B)[feature_B.order] )) %>% 
    dplyr::mutate(feature_A = factor(feature_A, 
                                     ordered = TRUE,
                                     levels = unique(corr.df.trim$feature_A)[feature_A.order] )) 
  
  ### Plot Heatmap
  h1 <- 
    corr.df.trim.ordered %>% 
    mutate(siglabel = if_else(q < 0.20, "*", "")) %>%
    mutate(tissue_A = factor(tissue_A, levels = facet_ord)) %>%
    ggplot(aes(x = feature_A_obj, 
               y = fct_relevel(feature_B, unique(corr.df.trim$feature_B)[feature_B.order]), 
               fill = rho)) +
    geom_tile() + 
    geom_text(aes(label=siglabel), size=5 ,vjust = 0.77, color = "white") +
    labs(fill = "Spearman's correlation") +
    scale_y_discrete(position = "right") +
    scale_fill_distiller(palette = "RdBu") +
    # facet_wrap(~tissue_A, scales = "free_x", nrow = 1) +
    facet_grid(cols = vars(tissue_A), scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=60, hjust =1),
      legend.position = "top",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.margin = unit(c(1, 1, 1, 2), "cm")
    )
  
  print(h1)
  return(h1)
}

corr_heatmap_facet_v2 <- function(corr.df, 
                                  facet_ord = c("Jejunum", "Ileum", "MLN", "Colon")){
  
  # Calculate number of significant (P-VALUE < 0.05) associations
  #  and select top 30 features
  corr.df.top <- corr.df %>% 
    dplyr::filter(q < 0.2) %>%
    dplyr::group_by(feature_B) %>%
    dplyr::summarise(n = n()) %>%
    slice_max(n = 30, order_by = n, with_ties = F)
  
  
  # Filter correlation df for top 30 features
  corr.df.trim <- corr.df %>% 
    mutate(feature_A = as.character(feature_A)) %>% 
    filter(feature_B %in% corr.df.top$feature_B)
  
  # create dataframe matrix of Rho correlation values for distance functions
  rho.df <- 
    corr.df.trim %>% 
    dplyr::select(feature_B, feature_A, rho) %>% 
    pivot_wider(names_from = feature_B, values_from = rho, values_fill = NA) %>% 
    column_to_rownames(var = "feature_A")
  
  # Hierarchical clustering of Rho values for features & feature_A
  feature_A.dendro <- 
    as.dendrogram(hclust(d = dist(x = rho.df), method = "complete"))
  feature_A.dendro.plot <- ggdendrogram(data = feature_A.dendro, rotate = TRUE)
  feature_B.dendro <- 
    as.dendrogram(hclust(d = dist(x = as.data.frame(t(rho.df))), method = "complete"))
  feature_B.dendro.plot <- ggdendrogram(data = feature_B.dendro, rotate = TRUE)
  
  ### Reorder Heatmap axis using order of dendrograms
  feature_B.order <- order.dendrogram(feature_B.dendro)
  feature_A.order <- order.dendrogram(feature_A.dendro)
  corr.df.trim.ordered <- corr.df.trim %>% 
    dplyr::mutate(feature = factor(feature_B, 
                                   ordered = TRUE,
                                   levels = unique(corr.df.trim$feature_B)[feature_B.order] )) %>% 
    dplyr::mutate(feature_A = factor(feature_A, 
                                     ordered = TRUE,
                                     levels = unique(corr.df.trim$feature_A)[feature_A.order] )) 
  
  ### Plot Heatmap
  h1 <- 
    corr.df.trim.ordered %>% 
    mutate(siglabel = if_else(q < 0.20, "*", "")) %>%
    mutate(tissue_A = factor(tissue_A, levels = facet_ord)) %>%
    ggplot(aes(x = feature_A_obj, 
               y = fct_relevel(feature_B, unique(corr.df.trim$feature_B)[feature_B.order]), 
               fill = rho)) +
    geom_tile() + 
    geom_text(aes(label=siglabel), size=5 ,vjust = 0.77, color = "white") +
    labs(fill = "Spearman's correlation") +
    scale_y_discrete(position = "right") +
    scale_fill_distiller(palette = "RdBu") +
    # scale_fill_gradient2(low="#538cb8", high="#e21d1d",
    #                      midpoint = 0, limit = c(-1,1), space = "lab") +
    facet_grid(cols = vars(tissue_A), scales = "free_x", space = "free_x") +
    # facet_grid(tissue_A_group~tissue_A, scales = "free", space = "free", drop = TRUE) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=90, hjust =1),
      legend.position = "top",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.margin = unit(c(1, 1, 1, 2), "cm")
    )
  
  print(h1)
  return(h1)
}


corr_heatmap_facet_behav <- function(corr.df, 
                               facet_ord = c("PostWean", "Juvenile", "Adult")){
  
  # Calculate number of significant (P-VALUE < 0.05) associations
  #  and select top 30 features
  corr.df.top <- corr.df %>% 
    dplyr::filter(p < 0.05) %>%
    dplyr::group_by(feature_B) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n)) %>%
    top_n(n = 30, wt = n) %>% 
    slice_head(n = 30)
  
  # Filter correlation df for top 30 features
  corr.df.trim <- corr.df %>% 
    mutate(feature_A = as.character(feature_A)) %>% 
    filter(feature_B %in% corr.df.top$feature_B)
  
  # create dataframe matrix of Rho correlation values for distance functions
  rho.df <- 
    corr.df.trim %>% 
    dplyr::select(feature_B, feature_A, rho) %>% 
    pivot_wider(names_from = feature_B, values_from = rho, values_fill = NA) %>% 
    column_to_rownames(var = "feature_A")
  
  # Hierarchical clustering of Rho values for features & feature_A
  feature_A.dendro <- 
    as.dendrogram(hclust(d = dist(x = rho.df), method = "complete"))
  feature_A.dendro.plot <- ggdendrogram(data = feature_A.dendro, rotate = TRUE)
  feature_B.dendro <- 
    as.dendrogram(hclust(d = dist(x = as.data.frame(t(rho.df))), method = "complete"))
  feature_B.dendro.plot <- ggdendrogram(data = feature_B.dendro, rotate = TRUE)
  
  ### Reorder Heatmap axis using order of dendrograms
  feature_B.order <- order.dendrogram(feature_B.dendro)
  feature_A.order <- order.dendrogram(feature_A.dendro)
  corr.df.trim.ordered <- corr.df.trim %>% 
    dplyr::mutate(feature = factor(feature_B, 
                                   ordered = TRUE,
                                   levels = unique(corr.df.trim$feature_B)[feature_B.order] )) %>% 
    dplyr::mutate(feature_A = factor(feature_A, 
                                     ordered = TRUE,
                                     levels = unique(corr.df.trim$feature_A)[feature_A.order] )) 
  
  ### Plot Heatmap
  h1 <- 
    corr.df.trim.ordered %>% 
    mutate(siglabel = if_else(q < 0.20, "*", "")) %>%
    mutate(timepoint = factor(timepoint, levels = facet_ord)) %>%
    ggplot(aes(x = feature_A_obj, 
               y = fct_relevel(feature_B, unique(corr.df.trim$feature_B)[feature_B.order]), 
               fill = rho)) +
    geom_tile() + 
    geom_text(aes(label=siglabel), size=5 ,vjust = 0.77, color = "white") +
    labs(fill = "Spearman's correlation") +
    scale_y_discrete(position = "right") +
    scale_fill_distiller(palette = "RdBu") +
    facet_grid(cols = vars(timepoint), scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=60, hjust =1),
      legend.position = "top",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.margin = unit(c(1, 1, 1, 2), "cm")
    )
  
  print(h1)
  return(h1)
}



# _______________________________________________________________________________

corr_heatmap_FDR <- function(df){
  df %>% 
    group_by(feature_A) %>%
    mutate(q = p.adjust(p, method = 'BH')) %>%
    ungroup() %>% 
    separate(feature_A, c("feature_A_obj", "tissue_A"), sep = "__", remove = F) %>% 
    separate(feature_B, c("feature_B_obj", "tissue_B"), sep = "__", remove = F)
}

corr_heatmap_FDRv2_plasma <- function(df){
  df %>% 
    group_by(feature_A) %>%
    mutate(q = p.adjust(p, method = 'BH')) %>%
    ungroup() %>% 
    mutate(feature_A = gsub("_y4", "", feature_A)) %>% 
    mutate(feature_A = toupper(feature_A)) %>% 
    separate(feature_A, c("stimulant", "feature_A_obj"), sep = "_", remove = F) %>% 
    separate(feature_B, c("feature_B_obj", "tissue_B"), sep = "__", remove = F) %>% 
    mutate(tissue_A = "Plasma") 
}

corr_heatmap_FDRv3_behavior <- function(df){
  df %>% 
    separate(feature_A, c("feature_A_obj", "timepoint"), sep = "__", remove = F) %>% 
    separate(feature_B, c("feature_B_obj", "tissue_B"), sep = "__", remove = F) %>% 
    group_by(feature_A, tissue_B) %>%
    mutate(q = p.adjust(p, method = 'BH')) %>%
    ungroup()
}

corr_FDR_16S <- function(df){
  df %>% 
    group_by(feature_A) %>%
    mutate(q = p.adjust(p, method = 'BH')) %>%
    ungroup() %>% 
    separate(feature_B, c("feature_B_obj", "tissue_B"), sep = "__", remove = F)
}

corr_FDR_16S_cyt <- function(df){
  df %>% 
    group_by(feature_A) %>%
    mutate(q = p.adjust(p, method = 'BH')) %>%
    ungroup() %>% 
    separate(feature_B, c("feature_B_obj", "tissue_B"), sep = "__", remove = F) %>% 
    separate(object_name, c("correlation_objects", "condition", "tissue"), sep = "__", remove = F)
}


# _______________________________________________________________________________

top_n_scatterplots <- function(hitslist,
                               merged_data,
                               data_type,
                               folder_prefix = "figures/correlations/scatterplots/",
                               subdir = T) {
  dir.create(file.path(folder_prefix), showWarnings = FALSE)
  
  for (corr in 1:nrow(hitslist)) {
    cor_row <- hitslist[corr, ]
    p <- corr_scatter_plot(
      df.plot = merged_data,
      corr_obj = hitslist,
      feature_var = cor_row$feature_A,
      metadata_var = cor_row$feature_B,
      data_type,
      folder_prefix,
      subdir
    )
    
  }
}

# _______________________________________________________________________________

corr_scatter_plot <- function(df.plot, corr_obj, feature_var, metadata_var, 
                              data_type, folder_prefix, subdir = T){
  
  #' Function creates a scatter plot of a given feature and a metadata column
  
  stat_col <-
    corr_obj %>%
    dplyr::filter(feature_A == sym(feature_var)) %>%
    dplyr::filter(feature_B == sym(metadata_var))
  stat_title <-
    paste0(
      "Spearman's Rho: ", round(stat_col$rho, digits = 3), "\n",
      "p-value: ", format(stat_col$p, scientific = T, digits = 3),
      "  FDR: ", format(stat_col$q, scientific = T, digits = 3)
    )
  cat("Rho: ", stat_col$rho, "\n")
  cat("P-value: ", format(stat_col$p, scientific = T, digits = 3), "\n")
  cat("Q-value: ", format(stat_col$q, scientific = T, digits = 3), "\n")
  
  plot <- 
    df.plot %>%
    mutate(group = factor(group, levels = group_ord)) %>% 
    drop_na(metadata_var) %>%
    ggplot(aes(x = .data[[feature_var]], y = .data[[metadata_var]])) +
    geom_point(aes(color = group), size = 1.5, alpha = 1) +
    geom_smooth(method = lm, color="darkgrey", linetype="dotted", se = F) +
    theme_bw() +
    labs(x = feature_var, 
         y = metadata_var, 
         title = stat_title,
         color = NULL) +
    scale_color_manual(values = group_cols) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 12))
  print(plot)
  figurename <- paste0(stat_col$feature_A[[1]], " vs ", stat_col$feature_B[[1]]) %>% 
    str_replace_all("[[:punct:]]", " ")
  if (subdir){
    plot.name = paste0(folder_prefix, stat_col$tissue_B, "/", data_type, figurename, ".svg")
  } else {
    plot.name = paste0(folder_prefix, "/", data_type, figurename, ".svg")
  }
  print(plot.name)
  ggsave(plot, filename = plot.name, height = 4, width = 5, )
}

# _______________________________________________________________________________
clean.cols.y4 <- function(x) {
  colnames(x) <- gsub("_y4", "", colnames(x)); x }
clean.cols.upper <- function(x) {
  colnames(x) <- toupper(colnames(x)); x }
# _______________________________________________________________________________

volcano_plot <- function(df, x, y, xlabel, colorvar){
  
  volc <- 
    df %>%
    ggplot(aes(x=x, y=y, color=colorvar)) + 
    geom_point(alpha = 0.9) + 
    xlab(bquote(log[2] ~ .(xlabel))) +
    ylab(expression(paste(-log[10], "[ P-value ]"))) +
    scale_colour_manual(values = c("A"= "#e31a1c", "B"= "#1f78b4", "C"= "black")) +
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
    geom_vline(xintercept = 1.00, linetype = 2, alpha = 0.5) +
    geom_vline(xintercept = -1.00, linetype = 2, alpha = 0.5) +
    theme_bw() +
    theme(legend.position="none",
          panel.grid = element_blank())
  return(volc)
}

geom_label <- function(df_labels, x, y, colorvar) {
  geom_label_repel(
    data = df_labels,
    aes(x = x, y = y, color = colorvar, label = `Biochemical Name`),
    segment.alpha = 0.75,
    segment.size = 0.2,
    segment.color = "grey50",
    # force = 1,
    max.time = 1,
    max.iter = Inf,
    max.overlaps = Inf
  )
}

# _______________________________________________________________________________
enrichment_formula <- function(N, n, m, k){
  
  #' N : the number of total metabolites in both the study and pathway library
  #' n : the number of significant metabolites in both the study and pathway library
  #' m : the number of detected metabolites in the pathway
  #' k : number of significant metabolites in the pathway
  
  cat("number of total metabolites: ", N, "\n")
  cat("number of total significant metabolites: ", n, "\n")
  cat("number of detected metabolites in the pathway: ", m, "\n")
  cat("number of significant metabolites in the pathway: ", k, "\n")
  
  i <- seq(from = 0, to = k, by = 1)
  binomal_coefs <- sum((choose(m, i)*choose((N-m), (n-i)) ) / choose(N, n))
  output <- 1 - binomal_coefs
  return(output)
}
# _______________________________________________________________________________
# _______________________________________________________________________________
# _______________________________________________________________________________

