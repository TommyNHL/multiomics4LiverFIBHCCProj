path <- getwd()
setwd(path)
print(path)
dataSource = "E:/LC_test/data/amidePos/caseControlStudy/"
dataSource = "E:/LC_test/ms2/QCs/neg/"
setwd(dataSource)

library(conflicted)
library(dplyr)
library(massprocesser)
library(tidyverse)
library(xcms)
library(MSnbase)
library(mzR)

# AmidePos, 0.50 min_fract AMCHCC 5-120, 0.95 min_fract AMCCIR 5-115, AMCFIB 5-95
# 0.50 min_fract 00-11fib 5-100, 11-21cirr 5-105, 21-31hcc 5-100
# 0.95 min_fract 00-10fibAMC 5-100, 10-20cirrAMC 5-105, 20-30hccAMC 5-100


# AmideNeg, 0.95 min_fract AMCHCC 5-115, AMCCIR 5-95, AMCFIB 5-45
# 0.50 min_fract 00-11fib 5-80, 11-21cirr 5-95, 21-31hcc 5-75
# 0.95 min_fract 00-10fibAMC 5-100,  

# T3Pos, 0.50 min_fract AMCHCC 5-50, AMCCIR 5-45, AMCFIB 5-40
# 0.95 min_fract 00-11fib 5-45, 11-21cirr 5-45, 21-31hcc 5-40
# 00-10fibAMC 5-35, 10-20cirrAMC 5-35, 20-30hccAMC 5-35
# 00-1121fibcirr 5-45
# 00-1020fibcirr 5-35

# T3Neg, 0.50 min_fract AMCHCC 5-40, AMCCIR 5-55, AMCFIB 5-4045
# 0.95 min_fract 00-11fib 5-45, 11-21cirr 5-45, 21-31hcc 5-30
# 00-10fibAMC 5-40, 10-20cirrAMC 5-40, 20-30hccAMC 5-40

massprocesser::process_data(#path = "./HKBUSKLP/.data/posDemoData/", 
                            #path = "./.data/posDemoData/",
                            path = dataSource, 
                            polarity = "negative",
                            ppm = 10,
                            peakwidth = c(5, 45),  # (5, 30)
                            snthresh = 4,
                            noise = 500,
                            threads = 6, 
                            output_tic = TRUE,
                            output_bpc = TRUE,
                            output_rt_correction_plot = TRUE,
                            min_fraction = 0.50,
                            fill_peaks = FALSE,
                            group_for_figure = "QC")

#load("./HKBUSKLP/.data/posDemoData/Result/object")  # on VS Code
#load("./.data/posDemoData/Result/object")  # on RStudio
load("./Result/object")  # on RStudio
object

load("./Result/intermediate_data/xdata2")
plot_rt = plot_adjusted_rt(object = xdata2, 
                           group_for_figure = "QC", 
                           interactive = TRUE)
plot_rt
