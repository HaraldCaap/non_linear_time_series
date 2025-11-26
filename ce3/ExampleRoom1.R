library(ctsmr)
library(splines)
source("sdeTiTm.R")
load("Exercise3.RData")
fit1 <- sdeTiTm(AllDat,AllDat$yTi1,AllDat$Ph1)

summary(fit1,extended=TRUE)

Hour <- as.numeric(strftime(AllDat$date, format="%H"))

Pred <- predict(fit1)
plot(Pred[[1]]$state$pred$Ti - AllDat$yTi1 ~ Hour, 
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Residuals",
     cex = 0.8)
abline(h = 0, col = "black", lwd = 2)
# What is going on 10 AM?
# Try to fit a varying effective window area


plot(AllDat$Gv ~ Hour)


idx <- (Hour>8 & Hour < 23) # It is impossible to fit a window area for the hours without any sun, so we limit the window area estimation to the hours with sun.
bs = bs(Hour[idx],df=5,intercept=TRUE) 

# What does the splines look like?
plot(bs[14:27,1],type='l')
lines(bs[ 14:27,2])
lines(bs[ 14:27,3])
lines(bs[ 14:27,4])
lines(bs[ 14:27,5])

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

plot(bs1[15:38] ~ Hour[15:38], type='l', 
     xlab="Hour",
     ylab="")
lines(bs2[15:38] ~ Hour[15:38])
lines(bs3[15:38] ~ Hour[15:38])
lines(bs4[15:38] ~ Hour[15:38])
lines(bs5[15:38] ~ Hour[15:38])


### You will have to implement sdeTiTmAv ###
source("sdeTiTmAv.R")
fit2 <- sdeTiTmAv(AllDat,AllDat$yTi1,AllDat$Ph1)
summary(fit2, extended = TRUE)

plot(bs[14:27,1]*fit2$xm[3]+bs[14:27,2]*fit2$xm[4]+bs[14:27,3]*fit2$xm[5]+bs[14:27,4]*fit2$xm[6]+bs[14:27,5]*fit2$xm[7],type='l')
plot(bs1*fit2$xm[3]+bs2*fit2$xm[4]+bs3*fit2$xm[5]+bs4*fit2$xm[6]+bs5*fit2$xm[7] ~ Hour,type='l', xlab="Hour of day", ylab="")


Pred2 <- predict(fit2)
plot(Pred2[[1]]$state$pred$Ti - AllDat$yTi4 ~ Hour, 
     pch = 1,                 # open circle
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Residuals",
     cex = 0.8)
abline(h = 0, col = "black", lwd = 2)

# Residuals
res2 <- Pred2[[1]]$state$pred$Ti - AllDat$yTi1
res1 <- Pred[[1]]$state$pred$Ti  - AllDat$yTi1

# Base plot: Pred residuals (open circles)
plot(Hour, res1,
     pch = 16,                 # filled circles
     col = "red",
     xlab = "Hour of Day",
     ylab = "Residuals",
     cex = 0.8)

# Add reference line
abline(h = 0, col = "black", lwd = 2)

# Add Pred2 residuals as additional points (filled circles or different color)
points(Hour, res2,
       pch = 16,              # open circle
       col = "blue",
       cex = 0.8)

# Add legend
legend("bottomleft",
       legend = c("Residuals original model", "Residuals with splines"),
       col = c("red", "blue"),
       pch = c(16, 16),
       pt.cex = c(0.8, 0.7),
       bty = "n")

MSE1   <- mean(res1^2, na.rm = TRUE)
MSE2   <- mean(res2^2, na.rm = TRUE)

## --------- Compute standard accuracy metrics and plot residuals ---------
y_pred <- Pred2[[1]]$state$pred$Ti   # predicted measurement
y_true <- AllDat$yTi1            # measured values
error <- y_pred - y_true
sd_pred <- Pred2[[1]]$state$sd$Ti
upper    <- y_pred + 2 * sd_pred
lower    <- y_pred - 2 * sd_pred

t <- AllDat$date[1:200]
# Start plot with true values as X marks
plot(t, y_true[1:200],
     pch = 4,                 # X marker
     col = "black",
     xlab = "Time",
     ylab = "Temperature [°C]",
     main = "True vs Predicted Temperature with Confidence Interval")
# Add predicted values as open circles
points(t, y_pred[1:200],
       pch = 1,               # open circle
       col = "blue")
# Add the upper and lower CI as lines
lines(t, upper[1:200], col = "blue", lwd = 1)
lines(t, lower[1:200], col = "blue", lwd = 1)
# Optional: Add a legend
legend("topright",
       legend = c("True (X)", "Predicted (O)", "95% CI"),
       pch    = c(4, 1, NA),
       lty    = c(0, 0, 1),
       col    = c("black", "blue", "blue"),
       pt.cex = 1.2,
       bty = "n")
