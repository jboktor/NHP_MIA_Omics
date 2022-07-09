

# 16S Analysis functions ----
#_______________________________________________________________________________

phyloseq_rarefy <- function(ps) {
  # Returns Rarefied Phyloseq Object
  obj <- ps %>%
    # remove taxa not detected in this tissue
    core(detection = 0, prevalence = 1 / 1000)
  df <- obj %>% abundances() %>% as.data.frame() %>% t()
  print(paste0(
    "A difference of ",
    max(rowSums(df)) - min(rowSums(df)),
    " reads is detected in sample sequencing depth"
  ))
  df_rarefied <-
    vegan::rrarefy(df, sample = min(rowSums(df))) %>% t() %>% as.data.frame()
  print(paste0("Reads rarefied to: ", max(colSums(df_rarefied))))
  ps_rare <- phyloseq(
    otu_table(df_rarefied, taxa_are_rows = TRUE),
    sample_data(meta(obj)),
    tax_table(obj),
    phy_tree(obj)
  )
  return(ps_rare)
}

#_______________________________________________________________________________


phyloseq_permanova <-
  function(ps_object,
           metadata_list,
           nperm = 10,
           dist = "Aitchisons") {
  #' Function to analyze phyloseq object with euclidean distance PERMANOVA
  #' Input: a phyloseq object and list of metadata columns to test
  #' Output: data frame with analysis variables
  require(phyloseq)
  require(microbiome)
  require(vegan)
  require(dplyr)
  require(foreach)
  require(doParallel)
  permanova_df <- tibble()
  start_time <- Sys.time()
  # Transform count data
  if (dist == "Aitchisons") {
    count_data <- ps_object %>% 
      microbiome::transform("clr") %>% 
      microbiome::abundances() %>% t()
  } else {
    count_data <- ps_object %>% 
      microbiome::abundances() %>% t()
  } 
  # pull metadata of interest from ps object
  metadata_vars <- microbiome::meta(ps_object) %>% 
    dplyr::select(all_of(metadata_list))
  # remove metadata with less than 1 unique value
  columns2keep <- sapply(metadata_vars, function(x) length(unique(na.omit(x)))) > 1
  metadata_vars <- metadata_vars[, columns2keep]
  #setup parallel processing
  start_time <- Sys.time()
  cores = detectCores()
  cl <- makeCluster(cores[1] - 2) # to prevent computer overload
  registerDoParallel(cl)
  loop <- 
    foreach(i = 1:length(metadata_vars), 
            .combine = 'rbind', .verbose = F,
            .packages = c("vegan", "phyloseq")) %dopar% {
              
              # filter NA values from metadata and abundance df
              a <- metadata_vars[,i]
              a.narm <- na.omit(a)
              if (any(is.na(a))) {
                count_data.narm <- count_data[-attr(a.narm, "na.action"), ]
              } else {
                count_data.narm <- count_data
              }
              # Calculate PERMANOVA
              if (dist == "Aitchisons"){
                meta_ano <- vegan::adonis(vegan::vegdist(
                  count_data.narm, method = "euclidean") ~ a.narm, permutations = nperm)
              } else if (dist == "uunifrac") {
                meta_ano <- vegan::adonis(phyloseq::distance(
                  ps_object, method="uunifrac") ~ a.narm, permutations = nperm)
              } else if (dist == "wunifrac") {
                meta_ano <- vegan::adonis(phyloseq::distance(
                  ps_object, method="wunifrac") ~ a.narm, permutations = nperm)
              } 
              
              # update stats df
              data.frame(
                "metadata" = colnames(metadata_vars[i]),
                "Df" = meta_ano$aov.tab[1, ]$Df,
                "SumsOfSqs" = meta_ano$aov.tab[1, ]$SumsOfSqs,
                "MeanSqs" = meta_ano$aov.tab[1, ]$MeanSqs,
                "F.Model" = meta_ano$aov.tab[1, ]$F.Model,
                "R2" = meta_ano$aov.tab[1, ]$R2,
                "p_value" = meta_ano$aov.tab[1, ]$`Pr(>F)`,
                "distance" = dist,
                "n" = length(a.narm),
                "permutations" = nperm
              )
            }
  stopCluster(cl)
  end_time <- Sys.time()
  cat("PERMANOVA calculated in : ",
      end_time - start_time, attr(end_time - start_time, "units"), "\n")
  
  return(loop)
}


#_______________________________________________________________________________


unifrac_plots <- function(obj_dist){
  
  # Ordination
  iDist_wu <- phyloseq::distance(obj_dist, method="wunifrac")
  iMDS_wu  <- phyloseq::ordinate(obj_dist, "MDS", distance=iDist_wu)
  p_wu <- plot_ordination(obj_dist, iMDS_wu, 
                          color="group", shape = "tissue", axes = c(1, 2))
  iDist_uu <- phyloseq::distance(obj_dist, method="uunifrac")
  iMDS_uu  <- phyloseq::ordinate(obj_dist, "MDS", distance=iDist_uu)
  p_uu <- plot_ordination(obj_dist, iMDS_uu, 
                          color="group", shape = "tissue", axes = c(1, 2))
  
  wu_plot <- 
    p_wu$data %>% 
    ggplot(aes(Axis.1, Axis.2, fill=group)) +
    geom_point(shape=21, size=2, alpha=0.7, stroke = 0.4) +
    xlab(paste0("PCoA 1 (", round((iMDS_wu$values$Relative_eig[1])*100, digits = 2), "%)")) +
    ylab(paste0("PCoA 2 (", round((iMDS_wu$values$Relative_eig[2])*100, digits = 2), "%)")) +
    labs(fill="", title = "Weighted UniFrac") +
    geom_xsideboxplot(aes(y = group), orientation = "y", width = 0.6) +
    geom_ysideboxplot(aes(x =group), orientation = "x", width = 0.6) +
    scale_xsidey_discrete(guide = "none") + scale_ysidex_discrete(guide = "none") +
    theme_bw() + 
    scale_fill_manual(values = group_cols) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          panel.grid = element_blank())
  
  uu_plot <- 
    p_uu$data %>% 
    ggplot(aes(Axis.1, Axis.2, fill=group)) +
    geom_point(shape=21, size=2, alpha=0.7, stroke = 0.4) +
    xlab(paste0("PCoA 1 (", round((iMDS_uu$values$Relative_eig[1])*100, digits = 2), "%)")) +
    ylab(paste0("PCoA 2 (", round((iMDS_uu$values$Relative_eig[2])*100, digits = 2), "%)")) +
    labs(fill="", title = "Unweighted UniFrac") +
    geom_xsideboxplot(aes(y = group), orientation = "y", width = 0.6) +
    geom_ysideboxplot(aes(x =group), orientation = "x", width = 0.6) +
    scale_xsidey_discrete(guide = "none") + scale_ysidex_discrete(guide = "none") +
    theme_bw() + 
    scale_fill_manual(values = group_cols) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())
  
  wu_plot + uu_plot
  
}


#_______________________________________________________________________________

alpha_div_metrics <- function(ps){
  core_meta <- meta(ps) %>% rownames_to_column(var = "id")
  ps %>% 
    microbiome::abundances() %>% 
    microbiome::alpha() %>% 
    rownames_to_column(var = "id") %>% 
    left_join(core_meta)  
}

#_______________________________________________________________________________

lm_basic <- function(metadf, metric){
  
  formula <-
    as.formula(paste(metric, "~", paste(
      c("gender", "treatment", "gender*treatment"), collapse = "+"
    )))
  model <- lm(formula, data = metadf, na.action = na.exclude)
  
  return(model)
}

#_______________________________________________________________________________



alpha_div_analysis <- function(df, info){
  
  alpha_metrics <- c(
    "observed",
    "diversity_shannon",
    "evenness_simpson",
    "dominance_gini",
    "rarity_low_abundance"
  )
  
  alpha_lm_stats <- tibble()
  alpha_anova_stats <- tibble()
  alpha_tukeyHSD_stats <- tibble()
  
  for (metric in alpha_metrics){
    print(metric)
    
    df %<>%
      mutate(group = factor(group, labels = names(group_cols))) %>% 
      mutate(treatment = factor(treatment, labels = c("Saline", "Poly(I:C)")))
    
    plot <- df %>%
      mutate(group = factor(group, labels = names(group_cols))) %>% 
      mutate(treatment = factor(treatment, labels = c("Saline", "Poly(I:C)"))) %>%
      ggplot(aes(x=group, y= .data[[metric]] )) +
      stat_summary(aes(fill=group),
                   geom="bar", fun.data=mean_se, color="black", width = 0.6) +
      stat_summary(fun.data=mean_se, fun.args = list(mult=1),
                   geom="errorbar", color="black", width=0.1) +
      geom_jitter(width=0.1, height=0, size = 2, alpha = 0.7) +
      labs(x = NULL, y = metric) +
      theme_bw() +
      scale_fill_manual(values = group_cols) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position="none")
    plot %>% print()
    ggsave(plot, filename = 
             paste0("figures/16S/alpha-diversity/", info, "-", metric, ".svg"), 
           width = 2.5, height = 4)
    
    lm_mod <- lm_basic(df, metric) 
    metric_lm <- lm_mod %>% tidy() %>% mutate(alpha_metric = metric, tissue = info)
    alpha_lm_stats %<>% bind_rows(metric_lm)
    
    aov_mod <- aov(lm_mod)
    summary(aov_mod) %>% print()
    metric_aov <- aov_mod %>% tidy() %>% mutate(alpha_metric = metric, tissue = info)
    alpha_anova_stats %<>% bind_rows(metric_aov)
    
    tukey_test <- TukeyHSD(aov_mod) %>% tidy() %>% mutate(alpha_metric = metric, tissue = info)
    alpha_tukeyHSD_stats %<>% bind_rows(tukey_test)
    
  }
  
  alpha_diversity_statistics <-
    list(alpha_lm_stats,
         alpha_anova_stats,
         alpha_tukeyHSD_stats)
  
  names(alpha_diversity_statistics) <- 
    c(paste0("linear_models-", info),
      paste0("anova-", info),
      paste0("tukeyHSD-", info))
  
  return(alpha_diversity_statistics)
  
}


#_______________________________________________________________________________
#                   Functions to adjust feature names              ----
#_______________________________________________________________________________


make_rfriendly_rows <- function(df, passed_column) {
  features <- rlang::sym(passed_column)
  df.out <-   df %>% 
    mutate(simplenames = gsub(":", ".gc.", !!features)) %>% 
    mutate(simplenames = gsub("\\|", ".gp.", simplenames)) %>% 
    mutate(simplenames = gsub(" ", ".gs.", simplenames)) %>% 
    mutate(simplenames = gsub("-", ".gh.", simplenames)) %>% 
    mutate(simplenames = gsub("/", ".gd.", simplenames)) %>% 
    mutate(simplenames = gsub("\\]", ".gsqrr.", simplenames)) %>% 
    mutate(simplenames = gsub("\\[", ".gsqrl.", simplenames)) %>% 
    mutate(simplenames = gsub("\\)", ".gpr.", simplenames)) %>% 
    mutate(simplenames = gsub("\\(", ".gpl.", simplenames)) %>% 
    mutate(simplenames = gsub(",", ".gm.", simplenames)) %>% 
    mutate(simplenames = gsub("\\+", ".gplus.", simplenames)) %>% 
    mutate(simplenames = gsub("\\'", ".gpar.", simplenames)) %>% 
    mutate(simplenames = paste0("feat_", simplenames)) %>% 
    tibble::column_to_rownames(var = "simplenames") %>% 
    dplyr::select(-!!features)
  return(df.out)
}


decode_rfriendly_rows <- function(df, passed_column) {
  features <- rlang::sym(passed_column)
  df.out <- df %>% 
    mutate(fullnames = gsub("\\.gc.", ":", !!features)) %>% 
    mutate(fullnames = gsub("\\.gsqrr.", "\\]", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gsqrl.", "\\[", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gplus.", "\\+", fullnames)) %>%
    mutate(fullnames = gsub("\\.gpar.", "\\'", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gpr.", ")", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gpl.", "(", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gm.", ",", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gp.", "\\|", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gs.", " ", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gh.", "-", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gd.", "/", fullnames)) %>% 
    mutate(fullnames = gsub("feat_", "", fullnames))
  return(df.out)
}

decoded_abundance <- function(datObj) {
  df <- datObj %>%
    microbiome::abundances() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    decode_rfriendly_rows(passed_column = "rowname") %>%
    column_to_rownames(var = "fullnames") %>%
    dplyr::select(-rowname)
  return(df)
}


#_______________________________________________________________________________


corr_heatmap_16s_to_mets <- function(corr.df, trimfeats = TRUE, featselection = "q"){

  
  if (trimfeats){
    if (featselection == "q"){ # most sig feats
      corr.df.top.mets <- corr.df %>% 
        dplyr::filter(q < 0.2) %>%
        slice_min(n = 30, order_by = q, with_ties = F)
      corr.df.top.bugs <- corr.df %>% 
        dplyr::filter(q < 0.2) %>%
        slice_min(n = 50, order_by = q, with_ties = F)
    } else if (featselection == "nq"){ # number of sig feats
      # Calculate number of significant (P-VALUE < 0.05) associations
      # and select top 30 features
      corr.df.top.mets <- corr.df %>%
        dplyr::filter(q < 0.2) %>%
        dplyr::group_by(feature_B_obj) %>%
        dplyr::summarise(n = n()) %>%
        slice_max(n = 30, order_by = n, with_ties = F)
      corr.df.top.bugs <- corr.df %>%
        dplyr::filter(q < 0.2) %>%
        dplyr::group_by(feature_A) %>%
        dplyr::summarise(n = n()) %>%
        filter(n > 0) %>%
        slice_max(n = 50, order_by = n, with_ties = F)
    }
    
    # Filter correlation df for top 30 features
    corr.df.trim <- corr.df %>% 
      mutate(feature_A = as.character(feature_A)) %>% 
      filter(feature_B_obj %in% corr.df.top.mets$feature_B_obj) %>%
      filter(feature_A %in% corr.df.top.bugs$feature_A)
  } else {
    corr.df.trim <- corr.df %>% mutate(feature_A = as.character(feature_A))
  }
  
  # create dataframe matrix of Rho correlation values for distance functions
  rho.df <- 
    corr.df.trim %>% 
    dplyr::select(feature_B_obj, asv_label_unique, rho) %>% 
    pivot_wider(names_from = feature_B_obj, values_from = rho, values_fill = NA) %>% 
    column_to_rownames(var = "asv_label_unique") %>% 
    drop_na()
  
  # Hierarchical clustering of Rho values for features & feature_A
  asv_label_unique.dendro <- 
    as.dendrogram(hclust(d = dist(x = rho.df), method = "complete"))
  feature_B_obj.dendro <- 
    as.dendrogram(hclust(d = dist(x = as.data.frame(t(rho.df))), method = "complete"))

  ### Reorder Heatmap axis using order of dendrograms
  feature_B_obj.order <- order.dendrogram(feature_B_obj.dendro)
  asv_label.order <- order.dendrogram(asv_label_unique.dendro)
  corr.df.trim.ordered <- corr.df.trim %>% 
    dplyr::mutate(feature_B_obj = factor(feature_B_obj, 
                                   ordered = TRUE,
                                   levels = unique(corr.df.trim$feature_B_obj)[feature_B_obj.order] )) %>% 
    dplyr::mutate(asv_label_unique = factor(asv_label_unique, 
                                     ordered = TRUE,
                                     levels = unique(corr.df.trim$asv_label_unique)[asv_label.order] )) 
  
  ### Plot Heatmap
  h1 <-
    corr.df.trim.ordered %>% 
    mutate(siglabel = if_else(q <= 0.20, "*", "")) %>%
    ggplot(aes(x=asv_label_unique, y=feature_B_obj, fill = rho)) +
    geom_tile() + 
    geom_text(aes(label=siglabel), size=5 ,vjust = 0.77, color = "white") +
    labs(fill = "Spearman's correlation") +
    scale_y_discrete(position = "right") +
    scale_fill_distiller(palette = "RdBu") +
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

#_______________________________________________________________________________

corr_heatmap_16s_to_cyts <-
  function(corr.df,
           trimfeats = T,
           featselection = "q",
           remove.nonsig = FALSE) {
  
  if (trimfeats){
    if (featselection == "q"){ # most sig feats
      corr.df.top.mets <- corr.df %>% 
        dplyr::filter(q < 0.2) %>%
        slice_min(n = 30, order_by = q, with_ties = F)
      corr.df.top.bugs <- corr.df %>% 
        dplyr::filter(q < 0.2) %>%
        slice_min(n = 50, order_by = q, with_ties = F)
    } else if (featselection == "nq"){ 
      # Calculate number of significant (P-VALUE < 0.05) associations
      # and select top 30 features
      corr.df.top.mets <- corr.df %>%
        dplyr::filter(q < 0.2) %>%
        dplyr::group_by(feature_B_obj) %>%
        dplyr::summarise(n = n()) %>%
        slice_max(n = 30, order_by = n, with_ties = F)
      corr.df.top.bugs <- corr.df %>%
        dplyr::filter(q < 0.2) %>%
        dplyr::group_by(feature_A) %>%
        dplyr::summarise(n = n()) %>%
        filter(n > 0) %>%
        slice_max(n = 50, order_by = n, with_ties = F)
    }
    
    # Filter correlation df for top 30 features
    corr.df.trim <- corr.df %>% 
      mutate(feature_A = as.character(feature_A)) %>% 
      filter(feature_B_obj %in% corr.df.top.mets$feature_B_obj) %>%
      filter(feature_A %in% corr.df.top.bugs$feature_A)
  } else {
    corr.df.trim <- corr.df %>% mutate(feature_A = as.character(feature_A)) %>% 
      drop_na()
  }
  
  # create dataframe matrix of Rho correlation values for distance functions
  rho.df <- 
    corr.df.trim %>% 
    dplyr::select(feature_B_obj, asv_label_unique, rho) %>% 
    pivot_wider(names_from = feature_B_obj, values_from = rho, values_fill = NA) %>% 
    column_to_rownames(var = "asv_label_unique") %>% 
    drop_na()
  
  # Hierarchical clustering of Rho values for features & feature_A
  asv_label_unique.dendro <- 
    as.dendrogram(hclust(d = dist(x = rho.df), method = "complete"))
  feature_B_obj.dendro <- 
    as.dendrogram(hclust(d = dist(x = as.data.frame(t(rho.df))), method = "complete"))
  
  ### Reorder Heatmap axis using order of dendrograms
  feature_B_obj.order <- order.dendrogram(feature_B_obj.dendro)
  asv_label.order <- order.dendrogram(asv_label_unique.dendro)
  corr.df.trim.ordered <- corr.df.trim %>% 
    dplyr::mutate(feature_B_obj = factor(feature_B_obj, 
                                         ordered = TRUE,
                                         levels = unique(corr.df.trim$feature_B_obj)[feature_B_obj.order] )) %>% 
    dplyr::mutate(asv_label_unique = factor(asv_label_unique, 
                                            ordered = TRUE,
                                            levels = unique(corr.df.trim$asv_label_unique)[asv_label.order] )) 
  
  if(remove.nonsig){
   corr.df.trim.ordered %<>% 
      mutate(rho = replace(rho, q>0.2, NA))
  }
  
  ### Plot Heatmap
  h1 <-
    corr.df.trim.ordered %>% 
    ggplot(aes(x=asv_label_unique, y=feature_B_obj,
               fill = rho)) +
    geom_tile() + 
    labs(fill = "Spearman's correlation") +
    scale_y_discrete(position = "right") +
    scale_fill_distiller(palette = "RdBu", na.value = "white") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
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
#_______________________________________________________________________________

prep16SAbundTable <- function(ps_obj){
  ps_obj_out <- ps_obj %>% core(detection = 0, prevalence = 0.1)
  sample_names(ps_obj_out) <- sample_data(ps_obj_out)$monkey
  ps_obj_out_df <- ps_obj_out %>% abundances() %>% t() %>% as.data.frame() 
  return(ps_obj_out_df)
}

#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
