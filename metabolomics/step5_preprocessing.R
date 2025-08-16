
dataSource = "E:/LC_test/data/t3Neg/progressionStudy/AMCFIB/"
setwd(dataSource)

load("./object_alignT3NegBatch")

library(massprocesser)
library(tidyverse)
library(xcms)
library(MSnbase)
library(mzR)

##-------------------------Process data--step one-------------------------------

posObject <- new_object
##-------------------- MS1 identification---------------------------------------
library(metid)

data("snyder_database_rplc0.0.3")
snyder_database_rplc0.0.3

# MS1 identification
posObject <-
  annotate_metabolites_mass_dataset(object = new_object,
                                    ms1.match.ppm = 10,
                                    rt.match.tol = 10000000,  # dysfunction by assigning a value of 10,000,000
                                    polarity = "negative",
                                    database = snyder_database_rplc0.0.3)

anotationTable <- posObject@annotation_table

write.csv(anotationTable, "./t31234NegAnnotationTable.csv")

##---------------------add MS2 to object----------------------------------------
library(massdataset)

posObject = mutate_ms2(object = posObject, 
                       column = "hilic",  # hilic
                       polarity = "positive",
                       ms1.ms2.match.mz.tol = 10,
                       ms1.ms2.match.rt.tol = 15,
                       path = "./MS2/")
#path = "../MS2_1/")

###----------------Remove noise features----------------------------------------
posObject %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::count(group)

# posObject <- 
#     posObject %>% 
#     activate_mass_dataset(what = "annotation_table") %>% 
#     group_by(Compound.name) %>% 
#     metid::filter(Level == min(Level)) %>% 
#     metid::filter(SS == max(SS)) %>% 
#     slice_head(n = 1)
# 
# posObject <- 
#     posObject %>% 
#     activate_mass_dataset(what = "annotation_table") %>% 
#     group_by(variable_id) %>% 
#     metid::filter(Level == min(Level)) %>% 
#     metid::filter(SS == max(SS)) %>% 
#     slice_head(n = 1)

qc_id = posObject %>% 
  activate_mass_dataset(what = "sample_info") %>%
  metid::filter(class == "QC") %>%
  pull(sample_id)

control_id = posObject %>%
  activate_mass_dataset(what = "sample_info") %>%
  metid::filter(group == "Control") %>%
  pull(sample_id)

case_id = posObject %>%
  activate_mass_dataset(what = "sample_info") %>%
  metid::filter(group == "Case") %>%
  pull(sample_id)

# case2_id =
#   object_pos1 %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   filter(group == "case2") %>%
#   pull(sample_id)

posObject = posObject %>%
  mutate_variable_na_freq(according_to_samples = qc_id) %>%
  mutate_variable_na_freq(according_to_samples = control_id) %>%
  mutate_variable_na_freq(according_to_samples = case_id)

head(extract_variable_info(posObject))

##-------------------------------remove variables-------------------------------
posObject <- posObject %>% 
  activate_mass_dataset(what = "variable_info") %>%
  metid::filter(na_freq < 0.2 & (na_freq.1 < 0.5 | na_freq.2 < 0.5))

# object_pos <- object_pos %>% 
#     activate_mass_dataset(what = "variable_info") %>%
#     filter(na_freq.3 < 0.5)

posObject

##------------------------------------Detect outlier----------------------------
library(MetNormalizer)
library(masscleaner)

outlier_samples = posObject %>%
  `+`(1) %>% 
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

outlier_table <- extract_outlier_table(outlier_samples)

outlier_table %>% head()

outlier_table %>% apply(1, function(x){sum(x)}) %>% `>`(0) %>% which()

##------------------------Missing value---------------------------------
get_mv_number(posObject)
# amidePos, 2687, 8542, 10893
# amideNeg, 4239, 13767, 
# t3Pos, 3247, 10502
# t3Neg, 3859, 21842

tryObject <- impute_mv(object = posObject, method = "knn")

get_mv_number(tryObject)

##---------------------normalization------------------------------------ 
tryObject <- normalize_data(tryObject, method = "svr")

object_pos_norm <- integrate_data(tryObject, method = "subject_median")
#object_pos_norm <- integrate_data(tryObject, method = "qc_mean")
#object_pos_norm <- tryObject

object_pos_norm %>% `+`(1) %>% log(2) %>% massqc::massqc_pca(color_by = "group", line = FALSE)


object_final4Merge <- object_pos_norm

save(object_final4Merge, file = "./object_final4Merge_T3Neg")
