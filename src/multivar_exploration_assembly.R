############################
#### Libraries and setup ###
############################

library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

setwd(dir = "/home/sid/Dev/MVDataExploration/data/")
options(scipen=1000000)

# principal dataset (filtered)
quast <- read.csv(file = "quast_all_metrics_reduced.csv", header=TRUE)
