# Global Metabolic Profiles in a Non-human Primate Maternal Immune Activation Model of Autism Spectrum Disorder
This R project integrates small molecule untargeted metabolomics, 16S sequencing, PBMC bioplex profiles, and behavioral data to produce analyses and publication figures.


## Layout

* data - analysis output files 
* src - R scripts analysis
* input_files - data for visualization and statistical analysis
* figures - data visualization


```
NHP_MIA_Omics Project 
 
├── NHP_MIA_Omics.Rproj
├── README.md
├── data
│   ├── correlations_behavior_metabolite
│   │   └── Behavior_Metabolite_Correlations.xlsx
│   ├── correlations_cytokine_metabolite
│   │   ├── Blood_Bioplex_Correlations.xlsx
│   │   ├── Brain_tissue_Bioplex_Correlations.xlsx
│   │   ├── GI_Bioplex_Correlations.xlsx
│   │   └── cytokine_correlation_pathway_enrichment_analysis.xlsx
│   └── correlations_metabolite-to-metabolite
│       ├── all_metabolite_corrs.xlsx
│       ├── subpathway_corrs.xlsx
│       └── superpathway_corrs.xlsx
├── figures
│   ├── boxplots
│   │   ├── ....
│   ├── correlations
│   │   ├── heatmaps
│   │   │   ├── ....
│   │   ├── heatmaps_by_tissue
│   │   │   ├── ....
│   │   ├── manual_scatterplots
│   │   │   ├── ....
│   │   ├── pathway_enrichment
│   │   │   ├──....
│   │   ├── scatterplots
│   │   │   ├── ....
│   │   └── scatterplots_tissue_specific
│   │       ├── ....
│   ├── ordination
│   │   ├── ....
│   ├── pathway_enrichment
│   │   ├── ....
│   └── volcano_plots
│       ├── ....
│      
├── input_files
│   ├── Behavior_Immune_data
│   │   ├── AMA\ 27\ GI\ tissues\ ConA\ 072513.xlsx
│   │   ├── AMA\ 27\ GI\ tissues\ PMA\ iono\ 073013.xlsx
│   │   ├── AMA\ 27\ GI\ tissues\ Poly\ IC\ 073013.xlsx
│   │   ├── AMA\ 27\ GI\ tissues\ media\ 072513.xlsx
│   │   ├── AMA27\ -\ brain\ cytokines.xlsx
│   │   ├── LongStereo\ -\ Year\ four\ stereotyped\ behaviors.xlsx
│   │   └── Year\ 4\ Blood\ Cytokine\ data.xlsx
│   ├── NHP\ CSF\ Metabolon.xlsx
│   ├── NHP\ Colon\ Metabolon.xlsx
│   ├── NHP\ Feces\ Metabolon.xlsx
│   ├── NHP\ Ileum\ Metabolon.xlsx
│   ├── NHP\ Jejunum\ Metabolon.xlsx
│   ├── NHP\ Pathway\ Enrichment\ Heatmaps.xlsx
│   ├── NHP\ Plasma\ HD4\ Metabolon.xlsx
│   ├── NHP\ Plasma\ Metabolon.xlsx
│   ├── Pathway\ Enrichment.xlsx
│   ├── TissueSharedMetabolites.xlsx
│   ├── pathway_map.xlsx
│   ├── sample_keys.xlsx
├── renv.lock
└── src
    ├── _load_packages.R
    ├── _misc_functions.R
    ├── boxplots.R
    ├── correlation_vis.R
    ├── correlations.R
    ├── pathway_enrichment.R
    ├── pathway_enrichment_simulations.R
    ├── summary_vis.R
    ├── tissue_shared_metabolites.R
    ├── tsne_pca.R
    └── volcano_plots.R
```


