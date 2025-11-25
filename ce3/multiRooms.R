## --------- Load libraries, source and data ---------
library(ctsmr)
library(splines)
load("Exercise3.RData")

Hour <- as.numeric(strftime(AllDat$date, format="%H"))

## Use only part of data, set n to 3111 to use all data
n <- 1000
source("sde4Rooms.R")
fit <- sde4Rooms(AllDat[1:n,]) 
Pred <- predict(fit, covariance = TRUE)

summary(fit, extended = TRUE)

##
y_preds <- matrix(c(Pred[[1]]$state$pred$T1[1:n], Pred[[1]]$state$pred$T3[1:n], Pred[[1]]$state$pred$T3[1:n], Pred[[1]]$state$pred$T4[1:n]), nrow=n, ncol=4)  # predicted measurement
y_trues <- matrix(c(AllDat$yTi1[1:n], AllDat$yTi2[1:n], AllDat$yTi3[1:n], AllDat$yTi4[1:n]), nrow=n, ncol=4)  # measured values
error <- y_preds - y_trues

# Basic residual plot with improved style
plot(Hour[1:n], error[1:n, 1],
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Residual (Prediction - True)",
     main = "Prediction Residuals",
     cex = 0.8)
# Horizontal zero line
abline(h = 0, col = "black", lwd = 2)