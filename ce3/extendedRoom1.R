# --------- Load libraries, source and data ---------
library(ctsmr)
library(splines)
load("Exercise3.RData")
source("sdeExRoom1.R")

# --------- Generate spline functions ---------
Hour <- as.numeric(strftime(AllDat$date, format="%H"))
idx <- (Hour>8 & Hour < 23) # It is impossible to fit a window area for the hours without any sun, so we limit the window area estimation to the hours with sun.
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

# --------- Fit model to data and calculate predictions ---------
fit <- sdeExRoom1(AllDat,AllDat$yTi1,AllDat$Ph1)
Pred <- predict(fit)
plot(Pred[[1]]$state$pred$Ti - AllDat$yTi1 ~ Hour)