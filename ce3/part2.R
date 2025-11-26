## --------- Load libraries, source and data ---------
library(ctsmr)
library(splines)
load("Exercise3.RData")

Hour <- as.numeric(strftime(AllDat$date, format="%H"))

