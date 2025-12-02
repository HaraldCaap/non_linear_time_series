library(ctsmr)
library(splines)

source("sde_models.R")

load("Exercise3.RData")
Hour <- as.numeric(strftime(AllDat$date, format="%H"))

############### Task 2a ###############


########### Original Fit 

fit_orig <- original_model(AllDat,AllDat$yTi1,AllDat$Ph1)
summary(fit_orig,extended=TRUE)
Pred_orig <- predict(fit_orig)

res_orig <- Pred_orig[[1]]$state$pred$Ti - AllDat$yTi1
RMSE_orig <- sqrt(mean(res_orig^2, na.rm = TRUE))
MAE_orig <- mean(abs(res_orig), na.rm = TRUE)
R2_orig <- 1 - sum(res_orig^2) / sum((AllDat$yTi1 - mean(AllDat$yTi1))^2)

# Plot residuals
plot(Pred_orig[[1]]$state$pred$Ti - AllDat$yTi1 ~ Hour, 
     pch = 1,
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Residuals",
     cex = 0.8)
abline(h = 0, col = "black", lwd = 2)


############### Splines

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
fit_splines <- splines_model(AllDat,AllDat$yTi1,AllDat$Ph1)

summary(fit_splines, extended = TRUE)
Pred_splines <- predict(fit_splines)


res_splines <- Pred_splines[[1]]$state$pred$Ti  - AllDat$yTi1

# Compute standard metrics
RMSE_splines <- sqrt(mean(res_splines^2, na.rm = TRUE))
MAE_splines <- mean(abs(res_splines), na.rm = TRUE)
R2_splines <- 1 - sum(res_splines^2) / sum((AllDat$yTi1 - mean(AllDat$yTi1))^2)

# Plot fitted splines
plot(Hour, bs1*fit_splines$xm[3]+bs2*fit_splines$xm[4]+bs3*fit_splines$xm[5]+bs4*fit_splines$xm[6]+bs5*fit_splines$xm[7],
     type='l', 
     xlab="Hour of day", 
     ylab="")

# Plot residulas
plot(Pred_splines[[1]]$state$pred$Ti - AllDat$yTi1 ~ Hour, 
     pch = 1,
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Residuals",
     cex = 0.8)
abline(h = 0, col = "black", lwd = 2)

## ----- Compare the two models -----

# Plot residuals in the same plot
plot(Hour, res_orig,
     pch = 16,
     col = "red",
     xlab = "Hour of Day",
     ylab = "Residuals",
     cex = 1)
abline(h = 0, col = "black", lwd = 2)
points(Hour, res_splines,
       pch = 16,
       col = "blue",
       cex = 0.7)
legend("bottomleft",
       legend = c("Residuals original model", "Residuals with splines"),
       col = c("red", "blue"),
       pch = c(16, 16),
       pt.cex = c(1, 0.7),
       bty = "n")

############### Task 2b ###############

plot(Hour, AllDat$yTi1,
     pch = 1,
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Temperature",
     cex = 0.8)
plot(Hour, AllDat$Ph1,
     pch = 1,
     col = "darkred",
     xlab = "Hour of Day",
     ylab = "Heating power",
     cex = 0.8)

########## Modelling with temp in room 2 
# Fit model and predict
fit_splines_neighbor <- splines_neighbor_model(AllDat,AllDat$yTi1,AllDat$Ph1, AllDat$yTi2)

summary(fit_splines_neighbor, extended = TRUE)
Pred_splines_neighbor <- predict(fit_splines_neighbor)


res_splines_neighbor <- Pred_splines_neighbor[[1]]$state$pred$Ti  - AllDat$yTi1

# Compute standard metrics
RMSE_splines_neighbor <- sqrt(mean(res_splines_neighbor^2, na.rm = TRUE))
MAE_splines_neighbor <- mean(abs(res_splines_neighbor), na.rm = TRUE)
R2_splines_neighbor <- 1 - sum(res_splines_neighbor^2) / sum((AllDat$yTi1 - mean(AllDat$yTi1))^2)



############ Heating Delay 

fit_splines_heating_delay <- splines_heating_delay_model(AllDat,AllDat$yTi1,AllDat$Ph1)

summary(fit_splines_heating_delay, extended = TRUE)
Pred_splines_heating_delay <- predict(fit_splines_heating_delay)


res_splines_heating_delay <- Pred_splines_heating_delay[[1]]$state$pred$Ti  - AllDat$yTi1

# Compute standard metrics
RMSE_splines_heating_delay <- sqrt(mean(res_splines_heating_delay^2, na.rm = TRUE))
MAE_splines_heating_delay <- mean(abs(res_splines_heating_delay), na.rm = TRUE)
R2_splines_heating_delay <- 1 - sum(res_splines_heating_delay^2) / sum((AllDat$yTi1 - mean(AllDat$yTi1))^2)

########### Solar Delay

fit_splines_solar_delay <- splines_solar_delay_model(AllDat,AllDat$yTi1,AllDat$Ph1)

summary(fit_splines_solar_delay, extended = TRUE)
Pred_splines_solar_delay <- predict(fit_splines_solar_delay)


res_splines_solar_delay <- Pred_splines_solar_delay[[1]]$state$pred$Ti  - AllDat$yTi1

# Compute standard metrics
RMSE_splines_solar_delay <- sqrt(mean(res_splines_solar_delay^2, na.rm = TRUE))
MAE_splines_solar_delay <- mean(abs(res_splines_solar_delay), na.rm = TRUE)
R2_splines_solar_delay <- 1 - sum(res_splines_solar_delay^2) / sum((AllDat$yTi1 - mean(AllDat$yTi1))^2)


###########  Radiator Temp
fit_splines_rad_temp <- splines_radiatorTemp_model(AllDat,AllDat$yTi1,AllDat$Ph1)

summary(fit_splines_rad_temp, extended = TRUE)
Pred_splines_rad_temp <- predict(fit_splines_rad_temp)

res_splines_rad_temp <- Pred_splines_rad_temp[[1]]$state$pred$Ti  - AllDat$yTi1

# Compute standard metrics
RMSE_splines_rad_temp <- sqrt(mean(res_splines_rad_temp^2, na.rm = TRUE))
MAE_splines_rad_temp <- mean(abs(res_splines_rad_temp), na.rm = TRUE)
R2_splines_rad_temp <- 1 - sum(res_splines_rad_temp^2) / sum((AllDat$yTi1 - mean(AllDat$yTi1))^2)


########## Plotting results

### --- Collect all residuals in one data frame ---

res_df <- data.frame(
  Hour = Hour,
  Spline = res_splines,
  Neighbor = res_splines_neighbor,
  HeatDelay = res_splines_heating_delay,
  SolarDelay = res_splines_solar_delay,
  RadTemp = res_splines_rad_temp
)

### --- Plot residuals of each model vs spline residuals ---

par(mfrow=c(2,2))

models <- c("Neighbor", "HeatDelay", "SolarDelay", "RadTemp")

for (m in models) {
  plot(res_df$Hour, res_df$Spline,
       pch = 16, col = "red",
       main = paste("Model:", m),
       xlab = "Hour", ylab = "Residuals")
  points(res_df$Hour, res_df[[m]], pch=16, col="blue", cex=0.8)
  legend("bottomleft", legend=c("Splines model", m), col=c("red", "blue"),pch=c(16,16), bty="n")
  #abline(h = 0, col = "black", lwd = 2)
}

## =========================
## 2) Sum of spline functions per model
## =========================

## =========================
## 2) Sum of spline functions per model
## =========================

## Design matrix for splines
B_spline <- cbind(bs1, bs2, bs3, bs4, bs5)

## Extract Aw1..Aw5 by fixed index
get_Aw <- function(fit) {
  fit$xm[3:7]
}

Aw_splines_neighbor       <- get_Aw(fit_splines_neighbor)
Aw_splines_heating        <- get_Aw(fit_splines_heating_delay)
Aw_splines_solar          <- get_Aw(fit_splines_solar_delay)
Aw_splines_radiator_temp  <- get_Aw(fit_splines_rad_temp)

## Sum of spline functions (effective window area term)
sum_spline_neighbor        <- as.numeric(B_spline %*% Aw_splines_neighbor)
sum_spline_heating_delay   <- as.numeric(B_spline %*% Aw_splines_heating)
sum_spline_solar_delay     <- as.numeric(B_spline %*% Aw_splines_solar)
sum_spline_radiator_temp   <- as.numeric(B_spline %*% Aw_splines_radiator_temp)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

plot(Hour, sum_spline_neighbor,        type="l", main="Neighbor model")
plot(Hour, sum_spline_heating_delay,   type="l", main="Heating delay model")
plot(Hour, sum_spline_solar_delay,     type="l", main="Solar delay model")
plot(Hour, sum_spline_radiator_temp,   type="l", main="Radiator temp model")


########## Extended part B model

fit_full <- combined_spline_delay_neighbor_model(AllDat,AllDat$yTi1,AllDat$Ph1, AllDat$yTi2)

summary(fit_full, extended = TRUE)
Pred_full <- predict(fit_full)


res_full <- Pred_full[[1]]$state$pred$Ti  - AllDat$yTi1

# Compute standard metrics
RMSE_full <- sqrt(mean(res_full^2, na.rm = TRUE))
MAE_full <- mean(abs(res_full), na.rm = TRUE)
R2_full <- 1 - sum(res_full^2) / sum((AllDat$yTi1 - mean(AllDat$yTi1))^2)

par(mfrow = c(1, 1))
plot(Hour, res_splines,
     pch = 16, col = "red",
     xlab = "Hour", ylab = "Residuals")
points(Hour, res_full, pch=16, col="blue", cex=0.8)
legend("bottomleft", legend=c("Splines model", "Extended model"), col=c("red", "blue"),pch=c(16,16), bty="n")

fit_full2 <- combined_spline_temp_delay_neighbor_model(AllDat,AllDat$yTi1,AllDat$Ph1, AllDat$yTi2)

summary(fit_full2, extended = TRUE)
Pred_full2 <- predict(fit_full2)


res_full2 <- Pred_full2[[1]]$state$pred$Ti  - AllDat$yTi1

# Compute standard metrics
RMSE_full2 <- sqrt(mean(res_full2^2, na.rm = TRUE))
MAE_full2 <- mean(abs(res_full2), na.rm = TRUE)
R2_full2 <- 1 - sum(res_full2^2) / sum((AllDat$yTi1 - mean(AllDat$yTi1))^2)


############### Task 2c ###############

#For part 2c we limit the amount of data for faster compute
nl <- 1
nu <- 3111

# Original 
fit_4rooms <- sde4Rooms_orig(AllDat[nl:nu,])

summary(fit_4rooms, extended = TRUE)

Pred_4rooms <- predict(fit_4rooms, covariance = TRUE)

# Residuals 
res_4rooms <- Pred_4rooms[[1]]$state$pred$T1 - AllDat$yTi1[nl:nu]

# Compute standard metrics
RMSE_4rooms <- sqrt(mean(res_4rooms^2, na.rm = TRUE))
MAE_4rooms <- mean(abs(res_4rooms), na.rm = TRUE)
R2_4rooms <- 1 - sum(res_4rooms^2) / sum((AllDat$yTi1[nl:nu] - mean(AllDat$yTi1[nl:nu]))^2)


res_4rooms_all <- c(
  Pred_4rooms[[1]]$state$pred$T1 - AllDat$yTi1[nl:nu],
  Pred_4rooms[[1]]$state$pred$T2 - AllDat$yTi2[nl:nu],
  Pred_4rooms[[1]]$state$pred$T3 - AllDat$yTi3[nl:nu],
  Pred_4rooms[[1]]$state$pred$T4 - AllDat$yTi4[nl:nu]
)

# Compute standard metrics
RMSE_4rooms_all <- sqrt(mean(res_4rooms_all^2, na.rm = TRUE))
MAE_4rooms_all <- mean(abs(res_4rooms_all), na.rm = TRUE)
R2_4rooms_all <- 1 - sum(res_4rooms_all^2) / sum((AllDat$yTi1[nl:nu] - mean(AllDat$yTi1[nl:nu]))^2 + (AllDat$yTi2[nl:nu] - mean(AllDat$yTi2[nl:nu]))^2 + (AllDat$yTi3[nl:nu] - mean(AllDat$yTi3[nl:nu]))^2+ (AllDat$yTi4[nl:nu] - mean(AllDat$yTi4[nl:nu]))^2)

# Individual residuals
res_T1 <- Pred_4rooms[[1]]$state$pred$T1 - AllDat$yTi1[nl:nu]
res_T2 <- Pred_4rooms[[1]]$state$pred$T2 - AllDat$yTi2[nl:nu]
res_T3 <- Pred_4rooms[[1]]$state$pred$T3 - AllDat$yTi3[nl:nu]
res_T4 <- Pred_4rooms[[1]]$state$pred$T4 - AllDat$yTi4[nl:nu]

# Put into long format for ggplot
N <- length(res_T1)
res_df <- data.frame(
  time = rep(seq_len(N), 4),
  residual = c(res_T1, res_T2, res_T3, res_T4),
  room = factor(rep(c("Room 1", "Room 2", "Room 3", "Room 4"), each = N))
)

# Plot
ggplot(res_df, aes(x = time, y = residual)) +
  geom_line(color = "steelblue") +
  facet_wrap(~ room, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Residuals in the 4 Rooms",
    x = "Time index",
    y = "Residual (Observed - Predicted)"
  )

#Floorplan

fit_4floorplan <- sde4Rooms_floorplan(AllDat[nl:nu,])

summary(fit_4floorplan, extended = TRUE)

Pred_4floorplan <- predict(fit_4floorplan, covariance = TRUE)

# Residuals 
room_names <- c("T1","T2","T3","T4")
obs_names  <- c("yTi1","yTi2","yTi3","yTi4")

res_list <- lapply(1:4, function(i){
  Pred_4floorplan[[1]]$state$pred[[ room_names[i] ]] -
    AllDat[[ obs_names[i] ]][nl:nu]
})

names(res_list) <- room_names

res_4floorplan_all <- c(
  Pred_4floorplan[[1]]$state$pred$T1 - AllDat$yTi1[nl:nu],
  Pred_4floorplan[[1]]$state$pred$T2 - AllDat$yTi2[nl:nu],
  Pred_4floorplan[[1]]$state$pred$T3 - AllDat$yTi3[nl:nu],
  Pred_4floorplan[[1]]$state$pred$T4 - AllDat$yTi4[nl:nu]
)

# Compute standard metrics
RMSE_4floorplan_all <- sqrt(mean(res_4floorplan_all^2, na.rm = TRUE))
MAE_4floorplan_all <- mean(abs(res_4floorplan_all), na.rm = TRUE)
R2_4floorplan_all <- 1 - sum(res_4floorplan_all^2) / sum((AllDat$yTi1[nl:nu] - mean(AllDat$yTi1[nl:nu]))^2 + (AllDat$yTi2[nl:nu] - mean(AllDat$yTi2[nl:nu]))^2 + (AllDat$yTi3[nl:nu] - mean(AllDat$yTi3[nl:nu]))^2+ (AllDat$yTi4[nl:nu] - mean(AllDat$yTi4[nl:nu]))^2)

####### Floorplan Heating Delay #######
fit_4floorplan_heatdelay <- sde4Rooms_floorplan_heatdelay(AllDat[nl:nu,])

summary(fit_4floorplan_heatdelay, extended = TRUE)

Pred_4floorplan_heatdelay <- predict(fit_4floorplan_heatdelay, covariance = TRUE)

# Residuals 
room_names <- c("T1","T2","T3","T4")
obs_names  <- c("yTi1","yTi2","yTi3","yTi4")

res_list <- lapply(1:4, function(i){
  Pred_4floorplan_heatdelay[[1]]$state$pred[[ room_names[i] ]] -
    AllDat[[ obs_names[i] ]][nl:nu]
})

names(res_list) <- room_names

res_4floorplan_heatdelay_all <- c(
  Pred_4floorplan_heatdelay[[1]]$state$pred$T1 - AllDat$yTi1[nl:nu],
  Pred_4floorplan_heatdelay[[1]]$state$pred$T2 - AllDat$yTi2[nl:nu],
  Pred_4floorplan_heatdelay[[1]]$state$pred$T3 - AllDat$yTi3[nl:nu],
  Pred_4floorplan_heatdelay[[1]]$state$pred$T4 - AllDat$yTi4[nl:nu]
)

# Compute standard metrics
RMSE_4floorplan__heatdelay_all <- sqrt(mean(res_4floorplan_heatdelay_all^2, na.rm = TRUE))
MAE_4floorplan__heatdelay_all <- mean(abs(res_4floorplan_heatdelay_all), na.rm = TRUE)
R2_4floorplan__heatdelay_all <- 1 - sum(res_4floorplan_heatdelay_all^2) / sum((AllDat$yTi1[nl:nu] - mean(AllDat$yTi1[nl:nu]))^2 + (AllDat$yTi2[nl:nu] - mean(AllDat$yTi2[nl:nu]))^2 + (AllDat$yTi3[nl:nu] - mean(AllDat$yTi3[nl:nu]))^2+ (AllDat$yTi4[nl:nu] - mean(AllDat$yTi4[nl:nu]))^2)

##### Original with heatdelay

fit_4rooms_orig_heatdelay <- sde4Rooms_orig_heatdelay(AllDat[nl:nu, ])

summary(fit_4rooms_orig_heatdelay, extended = TRUE)

Pred_4rooms_orig_heatdelay <- predict(fit_4rooms_orig_heatdelay, covariance = TRUE)

room_names <- c("T1","T2","T3","T4")
obs_names  <- c("yTi1","yTi2","yTi3","yTi4")

res_list <- lapply(1:4, function(i){
  Pred_4rooms_orig_heatdelay[[1]]$state$pred[[ room_names[i] ]] -
    AllDat[[ obs_names[i] ]][nl:nu]
})

names(res_list) <- room_names

res_4rooms_orig_heatdelay_all <- c(
  Pred_4rooms_orig_heatdelay[[1]]$state$pred$T1 - AllDat$yTi1[nl:nu],
  Pred_4rooms_orig_heatdelay[[1]]$state$pred$T2 - AllDat$yTi2[nl:nu],
  Pred_4rooms_orig_heatdelay[[1]]$state$pred$T3 - AllDat$yTi3[nl:nu],
  Pred_4rooms_orig_heatdelay[[1]]$state$pred$T4 - AllDat$yTi4[nl:nu]
)

RMSE_4rooms_orig_heatdelay <- sqrt(mean(res_4rooms_orig_heatdelay_all^2, na.rm = TRUE))

MAE_4rooms_orig_heatdelay <- mean(abs(res_4rooms_orig_heatdelay_all), na.rm = TRUE)

# Total variance across all 4 rooms for proper R²
total_variance <- sum(
  (AllDat$yTi1[nl:nu] - mean(AllDat$yTi1[nl:nu]))^2 +
    (AllDat$yTi2[nl:nu] - mean(AllDat$yTi2[nl:nu]))^2 +
    (AllDat$yTi3[nl:nu] - mean(AllDat$yTi3[nl:nu]))^2 +
    (AllDat$yTi4[nl:nu] - mean(AllDat$yTi4[nl:nu]))^2
)

R2_4rooms_orig_heatdelay <- 1 - sum(res_4rooms_orig_heatdelay_all^2) / total_variance

################ Plot results of 2c

res_T1_orig <- Pred_4rooms[[1]]$state$pred$T1 - AllDat$yTi1[nl:nu]
res_T2_orig <- Pred_4rooms[[1]]$state$pred$T2 - AllDat$yTi2[nl:nu]
res_T3_orig <- Pred_4rooms[[1]]$state$pred$T3 - AllDat$yTi3[nl:nu]
res_T4_orig <- Pred_4rooms[[1]]$state$pred$T4 - AllDat$yTi4[nl:nu]

res_T1_fp <- Pred_4floorplan_heatdelay[[1]]$state$pred$T1 - AllDat$yTi1[nl:nu]
res_T2_fp <- Pred_4floorplan_heatdelay[[1]]$state$pred$T2 - AllDat$yTi2[nl:nu]
res_T3_fp <- Pred_4floorplan_heatdelay[[1]]$state$pred$T3 - AllDat$yTi3[nl:nu]
res_T4_fp <- Pred_4floorplan_heatdelay[[1]]$state$pred$T4 - AllDat$yTi4[nl:nu]

res_df_compare <- data.frame(
  Hour = Hour[nl:nu],
  Orig_R1 = res_T1_orig,
  Orig_R2 = res_T2_orig,
  Orig_R3 = res_T3_orig,
  Orig_R4 = res_T4_orig,
  FP_R1 = res_T1_fp,
  FP_R2 = res_T2_fp,
  FP_R3 = res_T3_fp,
  FP_R4 = res_T4_fp
)

par(mfrow = c(2,2), mar = c(4,4,3,1))

rooms <- c("R1","R2","R3","R4")

par(mfrow = c(2,2), mar = c(4,4,3,1))

rooms <- c("R1","R2","R3","R4")

for (r in rooms) {
  
  orig_col <- paste0("Orig_", r)
  hd_col   <- paste0("FP_", r)
  
  # Plot ORIGINAL first (blue)
  plot(res_df_compare$Hour,
       res_df_compare[[orig_col]],
       pch = 16, col = "blue", cex = 1,
       xlab = "Hour of Day",
       ylab = "Residual",
       main = paste("Residuals – Room", substr(r,2,2)),
       ylim = range(
         c(res_df_compare[[orig_col]], 
           res_df_compare[[hd_col]]),
         na.rm = TRUE
       )
  )
  
  # Add HEAT-DELAY results second (red)
  points(res_df_compare$Hour,
         res_df_compare[[hd_col]],
         pch = 16, col = "red", cex = 0.7)

  
  legend("bottomleft",
         legend = c("Original 4 room model", "Area heat delay model"),
         col = c("blue", "red"), pch = 16, bty = "n")
}

