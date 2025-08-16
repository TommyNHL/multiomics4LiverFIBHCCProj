
dataSource = "E:/LC_test/ms2/QCs/neg/"
setwd(dataSource)

load("./Result/object")

library(massprocesser)
library(tidyverse)
library(xcms)
library(MSnbase)
library(mzR)

##-------------------------Process data--step one-------------------------------

posObject <- object
##-------------------- MS1 identification---------------------------------------
library(metid)

data("snyder_database_rplc0.0.3")
snyder_database_rplc0.0.3

# MS1 identification
posObject <-
  annotate_metabolites_mass_dataset(object = object,
                                    ms1.match.ppm = 10,
                                    rt.match.tol = 10000000,  # dysfunction by assigning a value of 10,000,000
                                    polarity = "negative",
                                    database = snyder_database_rplc0.0.3)

anotationTable <- posObject@annotation_table

write.csv(anotationTable, "./t3QCMS1NegAnnotationTable.csv")

##---------------------add MS2 to object----------------------------------------
library(massdataset)

posObject = mutate_ms2(object = posObject, 
                       column = "rp",  # hilic
                       polarity = "negative",
                       ms1.ms2.match.mz.tol = 10,
                       ms1.ms2.match.rt.tol = 15,
                       path = "../ms2/neg/")
#path = "../MS2_1/")

variable_info_pos <- 
  extract_variable_info(object = posObject)
head(variable_info_pos)

extract_ms2_data(posObject)

load("database_12FQE_20240502")
T3_database_12FQE_20240502

STD_table <- T3_database_12FQE_20240502@spectra.info
write.csv(STD_table, "STD_table.csv")

posObject2 <-
  annotate_metabolites_mass_dataset(object = posObject, 
                                    ms1.match.ppm = 10, 
                                    rt.match.tol = 20, 
                                    polarity = "negative",
                                    database = T3_database_12FQE_20240502)

anotationMS2_table <- posObject2@annotation_table

write.csv(anotationMS2_table, "t3QCMS2NegAnnotationTable.csv")

save(posObject2, file = "./t3NegObject")


ms2_plot_mass_dataset(object = posObject2, 
                      variable_id = "M191T97_1_NEG", 
                      database = T3_database_12FQE_20240502)
