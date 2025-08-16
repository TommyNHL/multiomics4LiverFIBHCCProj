
dataSource = "G:/LC_test/ms2/QCs/pos/"
setwd(dataSource)

load("./Result/object")

load("SKLP.database")

#library(massprocesser)
library(tidyverse)
library(xcms)
library(MSnbase)
library(mzR)

##-------------------------Process data--step one-------------------------------

posObject <- object
##-------------------- MS1 identification---------------------------------------
library(metid)

#data("snyder_database_rplc0.0.3")
#snyder_database_rplc0.0.3

# MS1 identification
posObject <-
  annotate_metabolites_mass_dataset(object = object,
                                    ms1.match.ppm = 10,
                                    rt.match.tol = 120,  # dysfunction by assigning a value of 10,000,000
                                    polarity = "positive",
                                    database = SKLP.database)

anotationTable <- posObject@annotation_table

write.csv(anotationTable, "./t3QCMS1PosAnnotationTable_new.csv")

##---------------------add MS2 to object----------------------------------------
library(massdataset)

posObject = mutate_ms2(object = posObject, 
                       column = "rp",  # hilic
                       polarity = "positive",
                       ms1.ms2.match.mz.tol = 10,
                       ms1.ms2.match.rt.tol = 20,
                       path = "../ms2/pos/")
#path = "../MS2_1/")

variable_info_pos <- 
  extract_variable_info(object = posObject)
head(variable_info_pos)

extract_ms2_data(posObject)


posObject <-
  annotate_metabolites_mass_dataset(object = posObject, 
                                    ms1.match.ppm = 10, 
                                    rt.match.tol = 20, 
                                    polarity = "positive", 
                                    database = SKLP.database)


anotationMS2_table <- posObject@annotation_table

write.csv(anotationMS2_table, "annotate_pos_new.csv")
#write.csv(posObject@annotation_table,"annotate_neg_new.csv")


# ms2_plot_mass_dataset(object = posObject, 
#                       variable_id = "M191T97_1_NEG", 
#                       database = T3_database_12FQE_20240502)


save(posObject, file = "./t3PosObject_new")

library(readr)
library(officer)
library(rvg)
library(ggplot2)
library(grid)
library(rsvg)
library(eoffice)
library(metid)

data <- read.csv("annotate_pos_new.csv")
ppt <- read_pptx() 

# 遍历每个variable_id
for (variable_id in unique(data$variable_id)) {
  # 生成图表
  plot <- ms2_plot_mass_dataset(object = posObject, polarity="positive",variable_id = variable_id, database = SKLP.database)
  
  # 将图表保存为临时文件
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, width =4400, height =3300, res = 600)  # 设置更高的分辨率
  print(plot)
  dev.off()
  
  # 添加新幻灯片
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
  
  # 将图表插入PPT
  ppt <- ph_with(ppt, external_img(temp_file), location = ph_location_fullsize())
}

output_file <- "plots_presentation_pos.pptx" 
print(ppt, target = output_file)
