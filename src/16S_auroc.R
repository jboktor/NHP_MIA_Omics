# Joe Boktor
# Caltech - Mazmanian Lab
# 2020

source("src/_load_packages.R")
source("src/_misc_functions.R")
source("src/_16S_functions.R")
library(DESeq2)
library(DEFormats)
library(limma)
library(edgeR)

ps_obj <- readRDS(file = "data/16S/phyloseq.rds")

aucs <- tibble()
start_time <- Sys.time()
for (tissue_origin in c("Stool", "Jejunum", "Ileum")) {
  ps_dat <- ps_obj %>%
    subset_samples(tissue == tissue_origin)
  ps_dat_list <- list(
    "ASVs" = ps_dat,
    "Phylum" = aggregate_taxa(ps_dat, "Phylum"),
    "Genus" = aggregate_taxa(ps_dat, "Genus"),
    "Species" = aggregate_taxa(ps_dat, "Species")
  )
  for (level in names(ps_dat_list)) {

    # _____________________________________________________________________________
    #####                    Limma - VOOM Normalization                     #####
    # _____________________________________________________________________________

    ps_object <- ps_dat_list[[level]]
    metadat <- meta(ps_object)
    print(ps_object)

    voom_matrix <- model.matrix(~ 0 + treatment + gender, data = metadat)

    dge <-
      ps_object %>%
      phyloseq_to_deseq2(~ treatment + gender) %>%
      as.DGEList()

    dge.norm <- dge %>% # dge.filtered %>%
      calcNormFactors(method = "TMM") %>%
      voom(
        design = voom_matrix, block = metadat$study, plot = TRUE,
        save.plot = TRUE, normalize.method = "none"
      )
    print(dim(t(dge.norm$E)))
    otu_table(ps_object) <- otu_table(dge.norm$E, taxa_are_rows = T)

    dat.polyic <- subset_samples(ps_object, treatment == "Poly(I:C)") %>%
      abundances()
    dat.control <- subset_samples(ps_object, treatment == "Saline") %>%
      abundances()

    require(foreach)
    require(doParallel)
    require(progress)

    iterations <- length(rownames(dat.polyic))
    pb <- progress_bar$new(
      format = "  calculating auroc [:bar] :current/:total :percent eta: :eta",
      total = iterations,
      width = 60
    )

    # allows progress bar to print within foreach
    progress_report <- function() {
      pb$tick()
    }
    opts <- list(progress = progress_report)

    cores <- detectCores()
    cl <- makeCluster(cores[1] - 2)
    registerDoParallel(cl)

    roc2add <-
      foreach(
        feat = rownames(dat.polyic),
        .combine = "rbind",
        .packages = c("magrittr", "pROC"),
        .options.snow = opts
      ) %dopar% {
        x <- dat.polyic[feat, ] %>% t()
        y <- dat.control[feat, ] %>% t()

        # AUROC
        rocdata <-
          c(roc(
            controls = y,
            cases = x,
            direction = "<",
            ci = TRUE,
            auc = TRUE
          )$ci)
        data.frame(
          "feature" = feat,
          "data_level" = level,
          "tissue" = tissue_origin,
          "ci_lower" = rocdata[1],
          "auroc" = rocdata[2],
          "ci_upper" = rocdata[3]
        )
      }

    stopCluster(cl)
    aucs <- rbind(aucs, roc2add)
  }
}


end_time <- Sys.time()
cat(
  "AUROCs calculated in : ",
  end_time - start_time, attr(end_time - start_time, "units"), "\n"
)

saveRDS(aucs, file = "data/16S/differential-abundance/feature_AUROCs_treatment.rds")

openxlsx::write.xlsx(aucs,
  file = "data/16S/differential-abundance/feature_AUROCs_treatment.xlsx",
  overwrite = T
)
