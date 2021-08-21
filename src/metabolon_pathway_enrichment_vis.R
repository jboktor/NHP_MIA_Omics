
source("src/_load_packages.R")
source("src/_misc_functions.R")

pthways <- read_xlsx('input_files/pathway_map.xlsx', sheet = 'Sheet 1') %>% 
  dplyr::rename("Superpathway" = "SUPER PATHWAY", "Subpathway" = "SUB PATHWAY") %>% 
  select(Superpathway, Subpathway) %>% distinct()
tissues <- c("CSF", "Plasma", "Jejunum", "Ileum", "Colon", "Feces")

#_______________________________________________________________________________
# Assemble Enrichment Analysis data

reportdir <- paste0("input_files/Metabolon_data/Enrichment_analysis/")
filepaths <- list.files(path =reportdir)
enrichment_treatment <- tibble()
enrichment_M_treatment <- tibble()
enrichment_F_treatment <- tibble()


for (file in filepaths){
  tissue_name <- str_split(file, ", ")[[1]][2] %>% str_split("_")
  df_report <- read_xlsx(paste0(reportdir, file)) %>% 
    mutate(tissue_source =  tissue_name[[1]][1])
  
  if (grepl(" male ", file)){
    print("Male Treatment File")
    enrichment_M_treatment <- enrichment_M_treatment %>% 
      bind_rows(df_report)
  } else if (grepl(" female ", file)){
    print("Female Treatment File")
    enrichment_F_treatment <- enrichment_F_treatment %>% 
      bind_rows(df_report)
  } else {
    print("Treatment File")
    enrichment_treatment <- enrichment_treatment %>% 
      bind_rows(df_report)
  }
}

prep_enrichment_stats <- function(df){
  df %>%
    dplyr::rename("Subpathway" = "Category") %>% 
    left_join(pthways, by = "Subpathway") %>% 
    mutate_at(c("k", "m", "n", "N", "Enrichment"), as.numeric) %>% 
    filter(m >= 5, Enrichment > 0)
}

enrichment_treatment <- prep_enrichment_stats(enrichment_treatment)
enrichment_F_treatment <- prep_enrichment_stats(enrichment_F_treatment)
enrichment_M_treatment <- prep_enrichment_stats(enrichment_M_treatment)

#_______________________________________________________________________________
##  Plotting functions

barplot_sums <- function(df, xlab, plot_title){
  df %>% 
    group_by(Superpathway, tissue_source) %>% 
    dplyr::summarise(score = sum(Enrichment)) %>% 
    ggplot() +
    geom_col(aes(y=factor(tissue_source, level = tissues), 
                 x=score, fill=Superpathway), 
             size = 10^-1, color = "black", width = 0.7) + 
    labs(x=xlab, title = plot_title) +
    scale_fill_manual(values = super_pathway_cols) +
    theme_classic() + 
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.95),
          axis.title.y =element_blank(),
          title = element_text(size = 8),
          axis.line.x.bottom = element_line(size=0.2),
          axis.line.y.left = element_line(size=0.2),
          legend.title = element_text(size = 10), 
          legend.text  = element_text(size = 8),
          legend.key.size = unit(0.7, "lines"))
}

barplot_sums_by_subpath <- function(df, xlab, plot_title){
  df %>% 
    group_by(Superpathway, tissue_source) %>% 
    dplyr::summarise(score = sum(Enrichment)) %>% 
    ggplot() +
    geom_col(aes(y=factor(tissue_source, level = tissues), 
                 x=score, fill=Superpathway), 
             size = 10^-1, color = "black", width = 0.7) + 
    labs(x=xlab, title = plot_title) +
    scale_fill_manual(values = super_pathway_cols) +
    theme_classic() + 
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.95),
          axis.title.y =element_blank(),
          title = element_text(size = 8),
          axis.line.x.bottom = element_line(size=0.2),
          axis.line.y.left = element_line(size=0.2),
          legend.title = element_text(size = 10), 
          legend.text  = element_text(size = 8),
          legend.key.size = unit(0.7, "lines"))
}

barplot_subpath_summary <- function(df, xlab, plot_title){
  
  subpathN <- df %>% pull(Subpathway) %>% unique() %>% length()
  print(paste0("  unique subpathways: ", subpathN))
  if(subpathN > 11){
    getPalette = colorRampPalette(brewer.pal(11, "RdYlBu"))
    fill_cols <- scale_fill_manual(values = getPalette(subpathN)) 
  } else {
    fill_cols <- scale_fill_brewer(palette = "RdYlBu") 
  }
  
  df %>% 
    ggplot() + 
    geom_col(aes(y=factor(tissue_source, level = tissues), 
                 x=Enrichment, fill=Subpathway), 
             size = 10^-1, color = "black", width = 0.7) + 
    labs(x=xlab, title = plot_title) +
    fill_cols +
    theme_classic() + 
    # guides(fill = guide_legend(override.aes = list(size = 0.2))) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.95),
          axis.title.y =element_blank(),
          title = element_text(size = 8),
          axis.line.x.bottom = element_line(size=0.2),
          axis.line.y.left = element_line(size=0.2),
          legend.title = element_text(size = 8), 
          legend.text  = element_text(size = 8),
          legend.key.size = unit(0.7, "lines"))
}

barplot_subpath <- function(df, xlab, plot_title){
  
  df %>% 
    ggplot() + 
    geom_col(aes(y=reorder(Subpathway, Enrichment), x=Enrichment, fill=Superpathway),
             size = 10^-1, color = "black", width = 0.7) +
    labs(x=xlab, title = plot_title) +
    # scale_fill_brewer(palette = "RdYlBu")  +
    scale_fill_manual(values = super_pathway_cols) +
    theme_classic() + 
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.95),
          axis.title.y =element_blank(),
          title = element_text(size = 8),
          axis.line.x.bottom = element_line(size=0.2),
          axis.line.y.left = element_line(size=0.2),
          legend.title = element_text(size = 10), 
          legend.text  = element_text(size = 10),
          legend.key.size = unit(0.7, "lines"))
}

barplot_subpath_all_tissues <- function(df, xlab, plot_title){
  
  df %>% 
    group_by(Subpathway) %>% 
    dplyr::summarise(score = sum(Enrichment)) %>% 
    dplyr::filter(score >= 5) %>% 
    left_join(pthways, by = "Subpathway") %>% 
    ungroup() %>% 
    ggplot() +
    geom_col(aes(y=reorder(Subpathway, score), 
                 x=score, fill=Superpathway), 
             size = 10^-1, color = "black", width = 0.7) + 
    labs(x=xlab, title = plot_title) +
    scale_fill_brewer(palette = "RdYlBu") +
    theme_classic() + 
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.95),
          axis.title.y =element_blank(),
          title = element_text(size = 8),
          axis.line.x.bottom = element_line(size=0.2),
          axis.line.y.left = element_line(size=0.2),
          legend.title = element_text(size = 10), 
          legend.text  = element_text(size = 10),
          legend.key.size = unit(0.7, "lines"))
}


# Save legend and plot as separate files to fix width:
save_me_pls <- function(ggobj, filename, plot_w, plot_h, leg_w, leg_h){
  
  require(cowplot)
  require(ggpubr)
  
  my_legend <- get_legend(ggobj) %>% as_ggplot()
  ggobj_out <- ggobj + theme(legend.position = "none")
  plot_filename <- paste0(filename, "_PLOT.svg")
  legend_filename <- paste0(filename, "_LEGEND.svg")
  
  ggsave(ggobj_out, filename = plot_filename, width = plot_w, height = plot_h)
  ggsave(my_legend, filename = legend_filename, width = leg_w, height = leg_h)
  
}


#_______________________________________________________________________________


f1d <- barplot_sums(df = enrichment_treatment, 
                    xlab = "Summated Subpathway Enrichment Poly(I:C)/Control", 
                    plot_title = " ")
fS1c <- barplot_sums(df = enrichment_M_treatment, 
                     xlab = "Summated Subpathway Enrichment Poly(I:C)/Control", 
                     plot_title = "Superpathway Enrichment in Male Treatment Group")
fS1d <- barplot_sums(df = enrichment_F_treatment,
                     xlab = "Summated Subpathway Enrichment Poly(I:C)/Control",
                     plot_title = "Superpathway Enrichment in Female Treatment Group")

save_me_pls(ggobj = f1d, 
            filename = "figures/pathway_enrichment/Fig1D_Combined_SuperPathwayByTissue", 
            plot_w = 4.5, plot_h = 2.5, leg_w = 5, leg_h = 5)
save_me_pls(ggobj = fS1c, 
            filename = "figures/pathway_enrichment/Fig1SC_Male_SuperPathwayByTissue", 
            plot_w = 4.5, plot_h = 2.5, leg_w = 5, leg_h = 5)
save_me_pls(ggobj = fS1d, 
            filename = "figures/pathway_enrichment/Fig1SD_Female_SuperPathwayByTissue", 
            plot_w = 4.5, plot_h = 2.5, leg_w = 5, leg_h = 5)



for (superpath in unique(enrichment_treatment$Superpathway)){
  a <- barplot_subpath_summary(filter(enrichment_treatment, Superpathway==superpath),
                          xlab = "Summated Subpathway Enrichment Poly(I:C)/Control", 
                          plot_title = " ")
  b <- barplot_subpath_summary(filter(enrichment_M_treatment, Superpathway==superpath),
                           xlab = "Summated Subpathway Enrichment Poly(I:C)/Control",
                           plot_title = paste0(superpath, " Superpathway Enrichment in Male Treatment Group"))
  c <- barplot_subpath_summary(filter(enrichment_F_treatment, Superpathway==superpath),
                           xlab = "Summated Subpathway Enrichment Poly(I:C)/Control",
                           plot_title = paste0(superpath, " Superpathway Enrichment in Female Treatment Group"))

  save_me_pls(ggobj = a, 
              filename = paste0("figures/pathway_enrichment/",superpath, "_all_Subpathway"), 
              plot_w = 4.5, plot_h = 2.5, leg_w = 5, leg_h = 5)
  save_me_pls(ggobj = b, 
              filename = paste0("figures/pathway_enrichment/",superpath, "_male_Subpathway"), 
              plot_w = 4.5, plot_h = 2.5, leg_w = 5, leg_h = 5)
  save_me_pls(ggobj = c, 
              filename = paste0("figures/pathway_enrichment/",superpath, "_female_Subpathway"), 
              plot_w = 4.5, plot_h = 2.5, leg_w = 5, leg_h = 5)

}

FIG_2A <- barplot_subpath(filter(enrichment_treatment, tissue_source=="Plasma"), xlab = "Enrichment Score", plot_title = " ")
FIG_3A <- barplot_subpath(filter(enrichment_treatment, tissue_source=="CSF"), xlab = "Enrichment Score", plot_title = " ")
FIG_4A <- barplot_subpath(filter(enrichment_treatment, tissue_source=="Feces"), xlab = "Enrichment Score", plot_title = " ")
FIG_5A <- barplot_subpath(filter(enrichment_treatment, tissue_source=="Ileum"), xlab = "Enrichment Score", plot_title = " ")
FIG_5D <- barplot_subpath(filter(enrichment_treatment, tissue_source=="Jejunum"), xlab = "Enrichment Score", plot_title = " ")
FIG_S9A <- barplot_subpath(filter(enrichment_treatment, tissue_source=="Colon"), xlab = "Enrichment Score", plot_title = " ")

 
ggsave(FIG_3A, filename = paste0("figures/pathway_enrichment/FIG_3A.svg"),
       height = 2.5 , width = 7)
ggsave(FIG_4A, filename = paste0("figures/pathway_enrichment/FIG_4A.svg"),
       height = 2.5 , width = 6.5)
ggsave(FIG_5A, filename = paste0("figures/pathway_enrichment/FIG_5A.svg"),
       height = 3.5 , width = 8)
ggsave(FIG_5D, filename = paste0("figures/pathway_enrichment/FIG_5D.svg"),
       height = 2.5 , width = 7.5)
ggsave(FIG_S9A, filename = paste0("figures/pathway_enrichment/FIG_S9A.svg"),
       height = 3 , width = 8)

# Subpathways summarized across tissues
FIG_1C <- barplot_subpath_all_tissues(
  enrichment_treatment, xlab = "Sum mated Enrichment Score Poly(I:C)/Control", plot_title = " ")
FIG_S1A <- barplot_subpath_all_tissues(
  enrichment_M_treatment, xlab = "Summated Enrichment Score Poly(I:C)/Control", 
  plot_title = "Metabolic Pathway Enrichment in Male Treatment Condition")
FIG_S1B <- barplot_subpath_all_tissues(
  enrichment_F_treatment, xlab = "Summated Enrichment Score Poly(I:C)/Control", 
  plot_title = "Metabolic Pathway Enrichment in Female Treatment Condition")

ggsave(FIG_1C, filename = paste0("figures/pathway_enrichment/FIG_1C.svg"),
       height = 5 , width = 8)
ggsave(FIG_S1A, filename = paste0("figures/pathway_enrichment/FIG_S1A.svg"),
       height = 5 , width = 8)
ggsave(FIG_S1B, filename = paste0("figures/pathway_enrichment/FIG_S1B.svg"),
       height = 5 , width = 8)


