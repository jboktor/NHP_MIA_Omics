
source("src/_load_packages.R")
source("src/_misc_functions.R")

# This script generates tSNE and PCA ordination plots for all samples 

plasmaInd <- read_xlsx('input_files/NHP Plasma HD4 Metabolon.xlsx', sheet = 'ScaledImpData')
plasmaInd <- plasmaInd[-2:-3]
plasma_labs <- c("BIOCHEMICAL", paste0(colnames(plasmaInd[2:24]), " Plasma")) 
colnames(plasmaInd) <- plasma_labs

csfInd <- read_xlsx('input_files/NHP CSF Metabolon.xlsx', sheet = 'ScaledImpData')
csfInd <- csfInd[-2:-3]
csf_labs <- c("BIOCHEMICAL", paste0(colnames(csfInd[2:22]), " CSF")) 
colnames(csfInd) <- csf_labs

jejunumInd <- read_xlsx('input_files/NHP Jejunum Metabolon.xlsx', sheet = 'ScaledImpData')
jejunumInd <- jejunumInd[-2:-3]
Jej_labs <- c("BIOCHEMICAL", paste0(colnames(jejunumInd[2:22]), " Jejunum")) 
colnames(jejunumInd) <- Jej_labs

ileumInd <- read_xlsx('input_files/NHP Ileum Metabolon.xlsx', sheet = 'ScaledImpData')
ileumInd <- ileumInd[-2:-3]
Ile_labs <- c("BIOCHEMICAL", paste0(colnames(ileumInd[2:22]), " Ileum")) 
colnames(ileumInd) <- Ile_labs

colonInd <- read_xlsx('input_files/NHP Colon Metabolon.xlsx', sheet = 'ScaledImpData')
colonInd <- colonInd[-2:-3]
Col_labs <- c("BIOCHEMICAL", paste0(colnames(colonInd[2:22]), " Colon")) 
colnames(colonInd) <- Col_labs

fecesInd <- read_xlsx('input_files/NHP Feces Metabolon.xlsx', sheet = 'ScaledImpData')
fecesInd  <-  fecesInd[-2:-3]
Fec_labs <- c("BIOCHEMICAL", paste0(colnames(fecesInd[2:22]), " Feces")) 
colnames(fecesInd) <- Fec_labs

allMetabolites <- plasmaInd
allMetabolites <- allMetabolites %>% full_join(csfInd)
allMetabolites <- allMetabolites %>% full_join(jejunumInd)
allMetabolites <- allMetabolites %>% full_join(ileumInd)
allMetabolites <- allMetabolites %>% full_join(colonInd)
allMetabolites <- allMetabolites %>% full_join(fecesInd)
samples <- colnames(allMetabolites[-1])
allMetabolites <- t(allMetabolites)

# Making First Row into Header
colnames(allMetabolites) <- allMetabolites[1,]
allMetabolites <- allMetabolites[-1,]

# Replacing all NA's with 0
allMetabolites[is.na(allMetabolites)] <- 0
allMetabolites <- as.data.frame(allMetabolites)

colnames(allMetabolites)[0] <- c("Sample")
allMet <-
  mutate(allMetabolites, tissue = if_else(
    grepl("Plasma", rownames(allMetabolites)),
    "Plasma",
    if_else(
      grepl("CSF", rownames(allMetabolites)),
      "CSF",
      if_else(
        grepl("Jejunum", rownames(allMetabolites)),
        "Jejunum",
        if_else(
          grepl("Ileum", rownames(allMetabolites)),
          "Ileum",
          if_else(
            grepl("Colon", rownames(allMetabolites)),
            "Colon",
            if_else(grepl("Feces", rownames(allMetabolites)), "Feces", "N/A")
          )
        )
      )
    )
  ))

# Make data frame Numeric
df <- allMet %>% 
  as.data.framez() %>% 
  select(-tissue) %>%
  mutate_if(is.character, as.numeric) %>%
  as.matrix()

#____________________________________________
#                   tSNE             ----
# ____________________________________________

## Executing the algorithm on curated data
set.seed(42)
tsne_out <- Rtsne(df, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
exeTimeTsne<- system.time(Rtsne(df, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500))
df.tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], col = tissue)

tsne <-
  df.tsne_plot %>%
  ggplot() + 
  geom_point(aes(x = x, y = y, color = col)) +
  scale_color_manual("", values = tissue_colors) +
  theme_classic() +
  labs(x = "tSNE dimension 1", y = "tSNE dimension 2")
tsne

ggsave(tsne, filename = "figures/ordination/tSNE_all_tissues.svg", 
       width = 5, height = 4)

#____________________________________________
#                   PCA             ----
# ____________________________________________

pca <- prcomp(df, center = TRUE, scale. = TRUE)

df_out <- as.data.frame(pca$x)
df_out$donor <- sapply(strsplit(as.character(row.names(df)), "_"), "[[", 1 )

pca.plot <-
  pca$x %>%
  as.data.frame() %>%
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  scale_color_manual("", values = tissue_colors) +
  theme_classic()
pca.plot

  
