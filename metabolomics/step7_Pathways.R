library(conflicted)
library(dplyr)
library(massprocesser)
library(tidyverse)
library(xcms)
library(MSnbase)
library(mzR)
library(tidymass)
library(tidyverse)

dataSource = "E:/LC_test/data/"
setwd(dataSource)

load("./finalT3Object_1011")

dir.create(path = "pathway_enrichment", showWarnings = FALSE)

diff_metabolites <-
  finalT3Object_0010 %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  metid::filter(p_value_adjust < 0.05) %>% 
  extract_variable_info()

head(diff_metabolites)

# ---Load KEGG human pathway database---
data("kegg_hsa_pathway", package = "metpath")
kegg_hsa_pathway

get_pathway_class(kegg_hsa_pathway)

# ---Remove the disease pathways---
pathway_class = metpath::pathway_class(kegg_hsa_pathway)
head(pathway_class)

remain_idx = pathway_class %>% 
    unlist() %>%
    stringr::str_detect("Disease") %>%
    `!`() %>%
    which()

remain_idx

pathway_database = kegg_hsa_pathway[remain_idx]
pathway_database

KEGGpath_table <- pathway_database@pathway_name
write.csv(KEGGpath_table, "KEGGpath_tableTable.csv")

kegg_id <- diff_metabolites$KEGG.ID 
kegg_id <- kegg_id[!is.na(kegg_id)]
kegg_id

result <- enrich_kegg(query_id = kegg_id, 
                      query_type = "compound", 
                      id_type = "KEGG",
                      pathway_database = pathway_database, 
                      p_cutoff = 0.05, 
                      p_adjust_method = "none", 
                      threads = 3)

enrich_bar_plot(object = result,
                x_axis = "p_value",
                cutoff = 0.05)

enrich_scatter_plot(object = result, y_axis = "p_value", y_axis_cutoff = 0.05)

enrich_network(object = result, 
               point_size = "p_value", 
               p_cutoff = 0.05, 
               only_significant_pathway = TRUE)
