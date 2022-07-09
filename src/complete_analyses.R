


#' Running this script will generate nearly all figures, tables,
#' and statistical summaries in our analysis.


# _______________________________________________________________________________
# Load and process data ----
# _______________________________________________________________________________

source("src/load-data-metabolites.R")
source("src/load-data-cytokines.R")
source("src/load-data-behaviors.R")

# _____________________________________________________________________________
# Exploration of metabolome profiles ----
# _____________________________________________________________________________

source("src/tsne_pca.R")
source("src/summary_vis.R")

# _____________________________________________________________________________
# Visualizing metabolome alterations ----
# _____________________________________________________________________________

source("src/volcano_plots.R")
source("src/boxplots.R")
source("src/metabolon_pathway_enrichment_vis.R")

# _____________________________________________________________________________
# Cytokine and Behavioral Correlations ----
# _____________________________________________________________________________

source("src/correlations.R")
source("src/correlation_vis.R")

# _____________________________________________________________________________
# 16S Analyses ----
# _____________________________________________________________________________

source("src/16S_metadata-prep.R")
source("src/16S_load-data.R")
source("src/16S_diversity-analysis.R")
source("src/16S_auroc.R")
source("src/16S_differential-abundance.R")
source("src/16S_differential-abundance-plots.R")
source("src/16S_correlations.R")
source("src/16S_correlations-plots.R")
