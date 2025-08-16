# the R script below is adapted from Daniel Beiting and modified by me

# Load packages -----
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(gt)
library(DT)
library(plotly)
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
targets <- read_tsv("studydesign.txt")
sampleLabels <- targets$sample
load(file = "../myDGEList1011FilterNorm")
log2.cpm.filtered.norm <- edgeR::cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

# ==============================================================================

# Set up your design matrix ----
group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# ==============================================================================

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)

# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# ==============================================================================

# Contrast matrix ----
contrast.matrix <- makeContrasts(infection = disease - healthy,
                                 levels=design)

# ==============================================================================

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)

#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
write.fit(ebFit, file="../lmfit1011_results.txt")

# ==============================================================================

# TopTable to view DEGs -----
#myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=10, sort.by="logFC")
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=10990, sort.by="logFC")

# convert to a tibble
myTopHits.df <- myTopHits %>% 
    as_tibble(rownames = "geneID")

temp <- left_join(log2.cpm.filtered.norm.df, myTopHits.df)
write.csv(temp, "../myDGEList1011log2_filtered_normalized_stat.csv")

# ==============================================================================

# Volcano Plots ----
vplot <- ggplot(myTopHits.df) + 
    aes(y=-log10(P.Value), x=logFC, text = paste("Symbol:", geneID)) + 
    geom_point(size=2) + 
    geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=0.8) + 
    geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=0.8) + 
    geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=0.8) + 
    annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.05), ymax = 6.5, alpha=.2, fill="#BE684D") + 
    annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.05), ymax = 6.5, alpha=.2, fill="#2C467A") + 
    labs(title="RNA-Seq: Liver Fibrosis & Healthy Groups") + 
    theme_bw()

EnhancedVolcano(myTopHits.df, 
                x="logFC", 
                y="P.Value", 
                xlim = c(min(-12), max(12)), 
                ylim = c(0, max(-log10(myTopHits.df[['P.Value']]), na.rm = TRUE) + 0.1), 
                xlab = bquote(~Log[2] ~ "Fold Change"), 
                ylab = bquote(~-Log[10] ~italic(P)~"-values"), axisLabSize = 18, 
                #title = "", 
                subtitle = "RNA-Seq: Liver Fibrosis & Healthy Groups", 
                lab = myTopHits.df$geneID, 
                pCutoff = 5e-2, 
                FCcutoff = 1)

# Now make the volcano plot above interactive with plotly
#ggplotly(vplot)

# ==============================================================================

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="none", p.value=0.05, lfc=1)

# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include=c("up", "down"))

# ==============================================================================

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes)

#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

#write your DEGs to a file
write_csv(diffGenes.df, "../DiffGenes1011.csv")

# ==============================================================================

# create interactive tables to display your DEGs ----
datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in Liver Fibrosis', 
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100"))) %>% 
    formatRound(columns=c(2:11), digits=2)
