## --------- Load libraries, source and data ---------
library(ctsmr)
library(splines)
load("Exercise3.RData")

Hour <- as.numeric(strftime(AllDat$date, format="%H"))

## ---------- Part 2a ----------
## ------ Original model ------
# Fit model and predict
source("sdeTiTm.R")
fit_orig <- sdeTiTm(AllDat,AllDat$yTi1,AllDat$Ph1)
summary(fit_orig,extended=TRUE)
Pred_orig <- predict(fit_orig)

# Plot residuals
plot(Pred_orig[[1]]$state$pred$Ti - AllDat$yTi1 ~ Hour, 
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Residuals",
     cex = 0.8)
abline(h = 0, col = "black", lwd = 2)

## ----- Splines model -----
# Add splines
idx <- (Hour>8 & Hour < 23) # It is impossible to fit a window area for the hours without any sun, so we limit the window area estimation to the hours with sun.
bs = bs(Hour[idx],df=5,intercept=TRUE) 

bs1 <- bs2 <- bs3 <- bs4 <- bs5 <- bs6 <- numeric(dim(AllDat)[1])

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

# Plot what the splines look like
plot(bs1[15:38] ~ Hour[15:38], type='l', 
     xlab="Hour",
     ylab="")
lines(bs2[15:38] ~ Hour[15:38])
lines(bs3[15:38] ~ Hour[15:38])
lines(bs4[15:38] ~ Hour[15:38])
lines(bs5[15:38] ~ Hour[15:38])

# Fit model and predict
source("sdeTiTmAv.R")
fit_splines <- sdeTiTmAv(AllDat,AllDat$yTi1,AllDat$Ph1)
summary(fit_splines, extended = TRUE)
Pred_splines <- predict(fit_splines)

# Plot fitted splines
plot(Hour, bs1*fit_splines$xm[3]+bs2*fit_splines$xm[4]+bs3*fit_splines$xm[5]+bs4*fit_splines$xm[6]+bs5*fit_splines$xm[7],
     type='l', 
     xlab="Hour of day", 
     ylab="")

# Plot residulas
Pred2 <- predict(fit2)
plot(Pred2[[1]]$state$pred$Ti - AllDat$yTi4 ~ Hour, 
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Residuals",
     cex = 0.8)
abline(h = 0, col = "black", lwd = 2)

## ----- Compare the two models -----
# Residuals
res_orig <- Pred_orig[[1]]$state$pred$Ti - AllDat$yTi1
res_splines <- Pred_splines[[1]]$state$pred$Ti  - AllDat$yTi1

# Plot residuals in the same plot
plot(Hour, res_orig,
     pch = 16,                 # filled circles
     col = "red",
     xlab = "Hour of Day",
     ylab = "Residuals",
     cex = 0.8)
abline(h = 0, col = "black", lwd = 2)
points(Hour, res_splines,
       pch = 16,              # open circle
       col = "blue",
       cex = 0.8)
legend("bottomleft",
       legend = c("Residuals original model", "Residuals with splines"),
       col = c("red", "blue"),
       pch = c(16, 16),
       pt.cex = c(0.8, 0.7),
       bty = "n")

# Compute standard metrics
MSE_orig <- mean(res_orig^2, na.rm = TRUE)
MSE_splines <- mean(res_splines^2, na.rm = TRUE)

MAE_orig <- mean(abs(res_orig), na.rm = TRUE)
MAE_splines <- mean(abs(res_splines), na.rm = TRUE)

R2_orig <- 1 - sum(res_orig^2) / sum((AllDat$yTi1 - mean(AllDat$yTi1))^2)
R2_splines <- 1 - sum(res_splines^2) / sum((AllDat$yTi1 - mean(AllDat$yTi1))^2)


## ----------- Part 2b -----------
# ----- EDA (Exploratory Data Analysis) ----- 
plot(Hour, AllDat$yTi1,
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Temperature",
     cex = 0.8)
plot(Hour, AllDat$Ph1,
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Heating power",
     cex = 0.8)
# There seems to be some delay in the heating power compared to the indoor temperature

