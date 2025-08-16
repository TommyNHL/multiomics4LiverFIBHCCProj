if(!require(devtools)){
  install.packages("devtools")
}

library(masscleaner)
library(tidyverse)
library(devtools)

path <- getwd()
setwd(path)
print(path)
dataSource1 = "E:/LC_test/data/t3Neg/progressionStudy/AMCFIB/AMCFIB1/"
setwd(dataSource1)
#load("./HKBUSKLP/.data/posDemoData/Result/object")  # on VS Code
#load("./.data/posDemoData/Result/object")  # on RStudio
load("./Result/object_norm1")  # on RStudio
b1 = object_norm1

dataSource2 = "E:/LC_test/data/t3Neg/progressionStudy/AMCFIB/AMCFIB234/"
setwd(dataSource2)
load("./Result/object_norm234") 
b234 = object_norm234

dataSource = "E:/LC_test/data/t3Neg/progressionStudy/AMCFIB/"
setwd(dataSource)

match_result = align_batch(x = b1, 
                           y = b234, 
                           combine.mz.tol = 15,
                           combine.rt.tol = 15,
                           use.int.tol = FALSE,
                           return_index = TRUE)

head(match_result)

new_object = align_batch(x = b1, 
                         y = b234, 
                         combine.mz.tol = 15,
                         combine.rt.tol = 15,
                         use.int.tol = FALSE,
                         return_index = FALSE)

save(new_object, file = "./object_alignT3NegBatch")
