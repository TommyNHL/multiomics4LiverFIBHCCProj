# the R script below is adapted from Daniel Beiting and modified by me

# Load packages ------
library(tidyverse) # you're familiar with this fromt the past two lectures
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables

# ==============================================================================

# Identify variables of interest in study design file ----
targets <- read_tsv("studydesign.txt")
targets
group <- targets$group
group <- factor(group)
sampleLabels <- targets$sample

# ==============================================================================

# Prepare your data -------
# for this part of the class you'll use your normalized and filtered data in log2 cpm
log2.cpm.filtered.norm.df <- read_csv("../myDGEList1011log2_filtered_normalized.csv")
log2.cpm.filtered.norm.df

load(file = "../myDGEList1011FilterNorm")
log2.cpm.filtered.norm <- edgeR::cpm(myDGEList.filtered.norm, log=TRUE)

# ==============================================================================

# Hierarchical clustering ---------------
distance <- dist(t(log2.cpm.filtered.norm), method = "euclidean") # "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
clusters <- hclust(distance, method = "complete") # "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
plot(clusters, labels=sampleLabels)

# ==============================================================================

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)

#look at the PCA result (pca.res) that you just created
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')

#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

# ==============================================================================

# Visualize your PCA result ------------------
pca.res.df <- as_tibble(pca.res$x)

ggplot(pca.res.df) + 
    aes(x=PC1, y=PC2, label=sampleLabels, color=group) + 
    geom_point(size=4) + 
    #geom_label() + 
    stat_ellipse() + 
    xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
    ylab(paste0("PC2 (",pc.per[2],"%",")")) + 
    labs(title="PCA plot", 
         #caption=paste0("produced on ", Sys.time())) + 
         caption="RNA-Seq: Liver Fibrosis & Healthy Groups") +
    coord_fixed() + 
    theme_bw()

# ==============================================================================

# Create a PCA 'small multiples' chart ----
pca.res.df <- pca.res$x[,1:4] %>%  # 80%, Magrittr package
    as_tibble() %>% 
    add_column(sample = sampleLabels, 
               group = group)
  
pca.pivot <- pivot_longer(
    pca.res.df, # dataframe to be pivoted
    cols = PC1:PC4, # column names to be stored as a SINGLE variable
    names_to = "PC", # name of that new variable (column)
    values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) + 
    aes(x=sample, y=loadings, fill=group) + 
    geom_bar(stat="identity") + 
    facet_wrap(~PC) + 
    labs(title="PCA 'small multiples' plot", 
         #caption=paste0("produced on ", Sys.time())) + 
         caption="RNA-Seq: Liver Fibrosis & Healthy Groups") + 
    theme_bw() + 
    coord_flip()

# ==============================================================================

# Use dplyr 'verbs' to modify our dataframe ----
# use dplyr 'mutate' function to add new columns based on existing data
mydata.df <- log2.cpm.filtered.norm.df %>% 
    mutate(healthy.AVG = (M1001R3 + M1004R1 + M1005R1 + M1010R1 + M1011R1 + M1012R1)/6, 
           disease.AVG = (M1101R1 + M1104R1 + M1105R1 + M1110R1 + M1111R1 + M1112R1)/6, 
           LogFC = (disease.AVG - healthy.AVG)) %>% 
    mutate_if(is.numeric, round, 2)

#now look at this modified data table
mydata.df

# Use dplyr 'arrange' and 'select' to sort your dataframe based on any variable
mydata.sort <- mydata.df %>% 
    dplyr::arrange(desc(LogFC)) %>% 
    dplyr::select(geneID, LogFC)

# Use dplyr "filter" and "select" functions to pick out genes of interest 
mydata.filter <- mydata.df %>% 
    dplyr::filter(geneID=="Gm10053" | geneID=="Rft1" | geneID=="App" | geneID=="Ano10" | geneID=="Tpt1-ps3"
                  | geneID=="Rps13-ps2" | geneID=="Grk6" | geneID=="Gpnmb" | geneID=="Mt2" ) %>% 
    dplyr::select(geneID, healthy.AVG, disease.AVG, LogFC) %>% 
    dplyr::arrange(desc(LogFC))

# you can also filter based on any regular expression
mydata.grep <- mydata.df %>% 
    dplyr::filter(grepl('Gm|Rft', geneID)) %>% 
    dplyr::select(geneID, healthy.AVG, disease.AVG, LogFC) %>% 
    dplyr::arrange(desc(geneID))

# ==============================================================================

# Produce publication-quality tables using the gt package ----
#gt(mydata.filter)

## now with a few more options
#mydata.filter %>% 
#    gt() %>% 
#    fmt_number(columns=2:4, decimals = 1) %>% 
#    tab_header(title = md("**Regulators of skin pathogenesis**"), 
#               subtitle = md("*during cutaneous leishmaniasis*")) %>% 
#    tab_footnote(
#        footnote = "Deletion or blockaid ameliorates disease in mice", 
#        locations = cells_body(
#            columns = geneID, 
#            rows = c(6, 7))) %>% 
#    tab_footnote(
#        footnote = "Associated with treatment failure in multiple studies", 
#        locations = cells_body(
#            columns = geneID, 
#            rows = c(2:9))) %>% 
#    tab_footnote(
#        footnote = "Implicated in parasite control", 
#        locations = cells_body(
#            columns = geneID, 
#            rows = c(2))) %>% 
#    tab_source_note(
#        source_note = md("Reference: Amorim *et al*., (2019). DOI: 10.1126/scitranslmed.aar3619"))

# ==============================================================================

# Make an interactive table using the DT package ----
#datatable(mydata.df[,c(1,12:14)], 
#          extensions = c('KeyTable', "FixedHeader"), 
#          filter = 'top', 
#          options = list(keys = TRUE, 
#                         searchHighlight = TRUE, 
#                         pageLength = 10, 
#                         lengthMenu = c("10", "25", "50", "100")))

# ==============================================================================

# Make an interactive scatter plot with plotly -----
## begin by storing your ggplot object
#myplot <- ggplot(mydata.df) + 
#    aes(x=healthy.AVG, y=disease.AVG) + 
#    geom_point(shape=16, size=1) + 
#    ggtitle("disease vs. healthy") + 
#    theme_bw()#

##now use the ggplotly function from the plotly package to convert this ggplot object into an interactive plot
#ggplotly(myplot)#

##let's customize this graphic by adding a more informative mouseover tooltip
#myplot <- ggplot(mydata.df) + 
#    aes(x=healthy.AVG, y=disease.AVG, 
#        text = paste("Symbol:", geneID)) + 
#    geom_point(shape=16, size=1) + 
#    ggtitle("disease vs. healthy") + 
#    theme_bw()#

#ggplotly(myplot)
