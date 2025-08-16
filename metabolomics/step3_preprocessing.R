path <- getwd()
setwd(path)
print(path)
dataSource = "E:/LC_test/data/t3Neg/case234ControlStudy/"
dataSource = "E:/LC_test/data/t3Neg/progressionStudy/AMCFIB/AMCFIB234/"
setwd(dataSource)

library(massprocesser)
library(tidyverse)
library(xcms)
library(MSnbase)
library(mzR)

##-------------------------Process data--step one-------------------------------
load("./Result/object")

posObject <- object

###----------------Remove noise features----------------------------------------
posObject %>%
    activate_mass_dataset(what = "sample_info") %>%
    dplyr::count(group)

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

##------------------------Missing value---------------------------------
library(MetNormalizer)
library(masscleaner)

get_mv_number(posObject)
# amidePos, 2687, 8542
# amideNeg, 4236, 12610 
# t3Pos, 3257, 10429
# t3Neg, 3690, 21933

tryObject <- impute_mv(object = posObject, method = "knn")

get_mv_number(tryObject)

##---------------------normalization------------------------------------ 
tryObject <- normalize_data(tryObject, method = "svr")

#object_pos_norm <- integrate_data(tryObject, method = "subject_median")
object_pos_norm <- integrate_data(tryObject, method = "qc_mean")

object_pos_norm %>% `+`(1) %>% log(2) %>% massqc::massqc_pca(color_by = "group", line = FALSE)

object_norm234 <- object_pos_norm

save(object_norm234, file = "./Result/object_norm234")
