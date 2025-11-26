## --------- Load libraries, source and data ---------
library(ctsmr)
library(splines)
load("Exercise3.RData")

Hour <- as.numeric(strftime(AllDat$date, format="%H"))

## --------- EDA (Exploratory Data Analysis) ---------
plot(Hour, AllDat$yTi1,
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Temperature",
     main = "Indoor temperature 1",
     cex = 0.8)

plot(Hour, AllDat$Gv,
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Solar Radiation",
     main = "Global Horizontal Solar Radiation",
     cex = 0.8)

plot(Hour, AllDat$Ta,
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Ambient temperature",
     main = "Ambient temperature",
     cex = 0.8)

plot(Hour, AllDat$Ph1,
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Heating power",
     main = "Heating power in northern circuit",
     cex = 0.8)

## --------- Generate spline functions ---------
idx <- (Hour>8 & Hour < 23)
bs = bs(Hour[idx],df=5,intercept=TRUE)

bs1 <- bs2 <- bs3 <- bs4 <- bs5 <- numeric(dim(AllDat)[1])

bs1[idx] = bs[,1]
bs2[idx] = bs[,2]
bs3[idx] = bs[,3]
bs4[idx] = bs[,4]
bs5[idx] = bs[,5]

AllDat$bs1 = bs1
AllDat$bs2 = bs2
AllDat$bs3 = bs3
AllDat$bs4 = bs4
AllDat$bs5 = bs5

plot(bs1[15:38] ~ Hour[15:38], type='l')
lines(bs2[15:38] ~ Hour[15:38])
lines(bs3[15:38] ~ Hour[15:38])
lines(bs4[15:38] ~ Hour[15:38])
lines(bs5[15:38] ~ Hour[15:38])

## --------- Fit model to data and calculate predictions ---------
## Use only part of data, set n to 3111 to use all data
n <- 3111
source("sdeExRoom1.R")
fit <- sdeExRoom1(AllDat[1:n,],AllDat$yTi1[1:n],AllDat$Ph1[1:n]) 
Pred <- predict(fit, covariance = TRUE)

summary(fit, extended = TRUE)

## --------- Compute standard accuracy metrics and plot residuals ---------
y_pred <- Pred[[1]]$state$pred$Ti[1:n]   # predicted measurement
y_true <- AllDat$yTi1[1:n]            # measured values
error <- y_pred - y_true
var_pred <- Pred[[1]]$state$var[1,1,]
upper    <- y_pred + 2 * sqrt(var_pred)
lower    <- y_pred - 2 * sqrt(var_pred)

# Basic residual plot
plot(Hour[1:n], error,
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Residual (Prediction - True)",
     main = "Prediction Residuals",
     cex = 0.8)
# Horizontal zero line
abline(h = 0, col = "black", lwd = 2)

t <- AllDat$date[1:n]
# Start plot with true values as X marks
plot(t, y_true,
     pch = 4,                 # X marker
     col = "black",
     xlab = "Time",
     ylab = "Temperature [°C]",
     main = "True vs Predicted Temperature with Confidence Interval")
# Add predicted values as open circles
points(t, y_pred,
       pch = 1,               # open circle
       col = "blue")
# Add the upper and lower CI as lines
lines(t, upper, col = "blue", lwd = 1)
lines(t, lower, col = "blue", lwd = 1)
# Optional: Add a legend
legend("topright",
       legend = c("True (X)", "Predicted (O)", "95% CI"),
       pch    = c(4, 1, NA),
       lty    = c(0, 0, 1),
       col    = c("black", "blue", "blue"),
       pt.cex = 1.2,
       bty = "n")

plot(Hour[1:n], var_pred, 
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Variance",
     main = "Prediction Variance",
     cex = 0.8)

  
MSE   <- mean(error^2, na.rm = TRUE)
MAE    <- mean(abs(error), na.rm = TRUE)
R2     <- 1 - sum(error^2) / sum((y_true - mean(y_true))^2)