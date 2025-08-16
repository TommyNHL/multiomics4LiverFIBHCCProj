# the R script below is adapted from Daniel Beiting and modified by me

# Load packages -----
library(tidyverse) # already know about this from Step 1 script
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure

# ==============================================================================

# Examine your data up to this point ----
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
colSums(myTPM)
colSums(myCounts)

# capture sample labels from the study design file that you worked with and saved as 'targets' in step 1
targets
sampleLabels <- targets$sample

# ==============================================================================

# Generate summary stats for your data ----
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))

# look at what you created
head(myTPM.stats)

# ==============================================================================

# Create your first plot using ggplot2 ----
# produce a scatter plot of the transformed data
ggplot(myTPM.stats) + 
    aes(x = SD, y = MED) +
    #geom_hex(shape=25, size=2) + 
    geom_point(shape=25, size=1) + 
    theme_bw()

# Let's expand on the plot above a bit more and take a look at each 'layer' of the ggplot code
p1 <- ggplot(myTPM.stats) + 
          aes(x = SD, y = MED) + 
          geom_point(shape=16, size=2) + 
          geom_smooth(method=lm) + 
          geom_hex(show.legend = FALSE) + 
          labs(y="Median", x = "Standard deviation", 
               title="Transcripts per million (TPM)", 
               subtitle="unfiltered, non-normalized data", 
               caption="RNA-Seq: Liver Fibrosis & Healthy Groups") + 
          theme_classic() + 
          theme_dark() + 
          theme_bw()
p1
# ==============================================================================

# Make a DGElist from your counts, and plot ----
myDGEList <- DGEList(myCounts)

# take a look at the DGEList object 
myDGEList

#DEGList objects are a good R data file to consider saving to you working directory
save(myDGEList, file = "../myDGEList1011")

#Saved DGEList objects can be easily shared and loaded into an R environment
load(file = "../myDGEList1011")

# use the 'cpm' function from EdgeR to get counts per million
cpm <- edgeR::cpm(myDGEList) 
colSums(cpm)
log2.cpm <- edgeR::cpm(myDGEList, log=TRUE)

# 'coerce' your data matrix to a dataframe so that you can use tidyverse tools on it
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
log2.cpm.df

# add your sample names to this dataframe (we lost these when we read our data in with tximport)
colnames(log2.cpm.df) <- c("geneID", sampleLabels)

# use the tidy package to 'pivot' your dataframe (from wide to long)
log2.cpm.df.pivot <- pivot_longer(
    log2.cpm.df, # dataframe to be pivoted
    cols = M1001R3:M1112R1, # column names to be stored as a SINGLE variable
    names_to = "samples", # name of that new variable (column)
    values_to = "expression") # name of new variable (column) storing all the values (data)

# let's look at the impact of pivoting the data
log2.cpm.df.pivot

# not it is easy to plot this pivoted data
p2 <- ggplot(log2.cpm.df.pivot) +
          aes(x=samples, y=expression, fill=samples) +
          geom_violin(trim = FALSE, show.legend = FALSE) +
          stat_summary(fun = "median", 
                       geom = "point", 
                       shape = 95, 
                       size = 10, 
                       color = "black", 
                       show.legend = FALSE) +
          labs(y="log2 expression", x = "sample",
               title="Log2 Counts per Million (CPM)",
               subtitle="unfiltered, non-normalized",
               #caption=paste0("produced on ", Sys.time())) + 
               caption="RNA-Seq: Liver Fibrosis & Healthy Groups") +
               theme_bw() + 
               theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
p2

# ==============================================================================

# Filter your data ----
#first, take a look at how many genes or transcripts have no read counts at all
table(rowSums(myDGEList$counts==0)==6)

# The line below is important! This is where the filtering starts
# Be sure to adjust this cutoff for the number of samples in the smallest group of comparison.
keepers <- rowSums(cpm>1)>=6

# now use base R's simple subsetting method to filter your DGEList based on the logical produced above
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered)

log2.cpm.filtered <- edgeR::cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)

# pivot this FILTERED data, just as you did earlier
log2.cpm.filtered.df.pivot <- pivot_longer(
    log2.cpm.filtered.df, # dataframe to be pivoted
    cols = -1, # column names to be stored as a SINGLE variable
    names_to = "samples", # name of that new variable (column)
    values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.df.pivot) + 
          aes(x=samples, y=expression, fill=samples) + 
          geom_violin(trim = FALSE, show.legend = FALSE) + 
          stat_summary(fun = "median", 
                       geom = "point", 
                       shape = 95, 
                       size = 10, 
                       color = "black", 
                       show.legend = FALSE) + 
          labs(y="log2 expression", x = "sample", 
               title="Log2 Counts per Million (CPM)", 
               subtitle="filtered, non-normalized", 
               #caption=paste0("produced on ", Sys.time())) + 
               caption="RNA-Seq: Liver Fibrosis & Healthy Groups") + 
               theme_bw() + 
               theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
p3

# ==============================================================================

# Normalize your data ----
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2.cpm.filtered.norm <- edgeR::cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

# pivot this NORMALIZED data, just as you did earlier
log2.cpm.filtered.norm.df.pivot <- pivot_longer(
    log2.cpm.filtered.norm.df, # dataframe to be pivoted
    cols = -1, # column names to be stored as a SINGLE variable
    names_to = "samples", # name of that new variable (column)
    values_to = "expression") # name of new variable (column) storing all the values (data)

p4 <- ggplot(log2.cpm.filtered.norm.df.pivot) + 
          aes(x=samples, y=expression, fill=samples) + 
          geom_violin(trim = FALSE, show.legend = FALSE) + 
          stat_summary(fun = "median", 
                       geom = "point", 
                       shape = 95, 
                       size = 10, 
                       color = "black", 
                       show.legend = FALSE) + 
          labs(y="log2 expression", x = "sample", 
               title="Log2 Counts per Million (CPM)", 
               subtitle="filtered, TMM normalized", 
               #caption=paste0("produced on ", Sys.time())) + 
               caption="RNA-Seq: Liver Fibrosis & Healthy Groups") + 
               theme_bw() + 
               theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
p4

# we'll use the 'plot_grid' function from the cowplot package to put these together in a figure
plot_grid(p1, p2, p3, p4, labels = c('A', 'B', 'C', 'D'), label_size = 12)

# save
write.csv(myDGEList@.Data[1], "../myDGEList1011.csv")
write.csv(log2.cpm.df, "../myDGEList1011log2.csv")

save(myDGEList.filtered, file = "../myDGEList1011Filter")
write.csv(myDGEList.filtered@.Data[1], "../myDGEList1011_filtered.csv")
write.csv(log2.cpm.filtered.df, "../myDGEList1011log2_filtered.csv")

save(myDGEList.filtered.norm, file = "../myDGEList1011FilterNorm")
write.csv(myDGEList.filtered.norm@.Data[1], "../myDGEList1011_filtered_normalized.csv")
write.csv(log2.cpm.filtered.norm.df, "../myDGEList1011log2_filtered_normalized.csv")

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
cpm.filtered.norm <- edgeR::cpm(myDGEList.filtered.norm, log=FALSE)
cpm.filtered.norm.df <- as_tibble(cpm.filtered.norm, rownames = "geneID")
colnames(cpm.filtered.norm.df) <- c("geneID", sampleLabels)
write.csv(cpm.filtered.norm.df, "../myDGEList1011log2FALSE_filtered_normalized.csv")
