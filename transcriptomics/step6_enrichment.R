# the R script below is adapted from Daniel Beiting and modified by me

# Load packages ----
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots

# Carry out GO enrichment using gProfiler2 ----
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits <- topTable(ebFit, adjust ="none", coef=1, number=10990, sort.by="logFC")

# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
gost.res <- gost(rownames(myTopHits), organism = "mmusculus", correction_method = "fdr")

# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.res, interactive = TRUE, capped = TRUE)

#set interactive=FALSE to get plot for publications
mygostplot <- gostplot(gost.res, interactive = FALSE, capped = TRUE)

# produce a publication quality static manhattan plot with specific GO terms highlighted
publish_gostplot(
    mygostplot, #your static gostplot from above
    highlight_terms = c(),
    filename = NULL,
    width = NA,
    height = NA)

# ==============================================================================

# Perform GSEA using clusterProfiler ----
# option: use the msigdb package to access up-to-date collections
# this option has the additional advantage of providing access to species-specific collections
msigdbr_species()
mus_gsea <- msigdbr(species = "Mus musculus") #gets all collections/signatures with mouse gene IDs

#take a look at the categories and subcategories of signatures available to you
mus_gsea %>% 
    dplyr::distinct(gs_collection, gs_subcollection) %>% 
    dplyr::arrange(gs_collection, gs_subcollection)

# choose a specific msigdb collection/subcollection
mus_gsea_c2 <- msigdbr(species = "Mus musculus", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
    dplyr::select(gs_name, gene_symbol) 
#just get the columns corresponding to signature name and gene symbols of genes in each signature

# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)

# construct a named vector
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
set.seed(123) #set a random seed so that we can reproducible ordering for our GSEA results below

#could replace C2CP with mus_gsea_c2 object you retrieved from msigdb above
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=mus_gsea_c2, verbose=FALSE)

myGSEA.df <- as_tibble(myGSEA.res@result)

# view results as an interactive table
datatable(myGSEA.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'RNA-Seq: Liver Fibrosis & Healthy Groups', 
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100"))) %>% 
    formatRound(columns=c(2:10), digits=2)

# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res,
          geneSetID = c(25, 100, 101, 86, 102, 139, 106, 149, 
                        178, 157, 210, 221, 196, 291, 337, 271, 
                        342, 417, 441, 437, 343), 
          pvalue_table = FALSE) #can set this to FALSE for a cleaner plot
          #title = myGSEA.res$Description[1]
          #can also turn off this title

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df[c(25, 100, 101, 86, 102, 139, 106, 149, 
                  178, 157, 210, 221, 196, 291, 337, 271, 
                  342, 417, 441, 437, 343),] %>% 
    mutate(phenotype = case_when(
        NES > 0 ~ "Liver Fibrosis", 
        NES < 0 ~ "Healthy"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df[1:21,], aes(x=phenotype, y=ID)) +
    geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) + 
    scale_color_gradient(low="blue", high="red") + 
    theme_bw()
