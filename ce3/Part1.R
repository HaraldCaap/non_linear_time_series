# Part 1: Bonhoeffer–Van der Pol SDE

library(ggplot2)
library(MASS)
library(dplyr)
library(viridisLite)

# Parameters
theta1 <- 0.7; theta2 <- 0.8; theta3 <- 3.0; theta4 <- -0.34
dt <- 2^-9; Tmax <- 100
N <- as.integer(Tmax / dt)
tgrid <- seq(0, by = dt, length.out = N + 1)

y10 <- -1.9; y20 <- 1.2
sigma_vals <- c(0, 0.1, 0.2, 0.3, 0.4)

set.seed(42)
Z <- rnorm(N)

# Euler–Maruyama
simulate_system <- function(sigma) {
  y1 <- numeric(N+1); y2 <- numeric(N+1)
  y1[1] <- y10; y2[1] <- y20
  sdt <- sqrt(dt)
  
  for (i in 1:N) {
    drift1 <- theta3 * (y1[i] + y2[i] - (1/3)*y1[i]^3 + theta4)
    drift2 <- -(1/theta3) * (y1[i] + theta2*y2[i] - theta1)
    y1[i+1] <- y1[i] + drift1*dt + sigma*sdt*Z[i]
    y2[i+1] <- y2[i] + drift2*dt
  }
  data.frame(t = tgrid, Y1 = y1, Y2 = y2, sigma = factor(sigma))
}

sim_df <- do.call(rbind, lapply(sigma_vals, simulate_system))
sim_df$idx <- ave(seq_len(nrow(sim_df)), sim_df$sigma, FUN = seq_along)
plot_df <- sim_df[sim_df$idx %% 10 == 1, ]

# ------------------
# Part 1a: Plots
# ------------------

ggplot(plot_df, aes(t, Y1)) +
  geom_line(size = 0.25) +
  facet_wrap(~ sigma, ncol = 1, scales = "free_y") +
  labs(title = "Y1 over time", x = "t", y = "Y1") +
  theme_minimal()

ggplot(plot_df, aes(t, Y2)) +
  geom_line(size = 0.25) +
  facet_wrap(~ sigma, ncol = 1, scales = "free_y") +
  labs(title = "Y2 over time", x = "t", y = "Y2") +
  theme_minimal()

ggplot(plot_df, aes(Y1, Y2)) +
  geom_path(size = 0.25) +
  facet_wrap(~ sigma, ncol = 3) +
  labs(title = "Phase plot", x = "Y1", y = "Y2") +
  theme_minimal()

# ------------------
# Part 1b: Densities
# ------------------

sigma_list <- c("0.1", "0.2", "0.3", "0.4")
sim_df$sigma <- as.character(sim_df$sigma)

global_x <- range(sim_df$Y1)
global_y <- range(sim_df$Y2)

dens_all <- do.call(rbind, lapply(sigma_list, function(s) {
  df_s <- sim_df[sim_df$sigma == s, ]
  kd <- kde2d(df_s$Y1, df_s$Y2, n = 200,
              lims = c(global_x[1], global_x[2],
                       global_y[1], global_y[2]))
  expand.grid(Y1 = kd$x, Y2 = kd$y, sigma = s) |>
    mutate(z = as.vector(kd$z))
}))

ggplot(dens_all, aes(Y1, Y2, fill = z)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  coord_fixed(xlim = global_x, ylim = global_y) +
  facet_wrap(~ sigma, ncol = 2) +
  labs(title = "Smoothed density for different σ", x = "Y1", y = "Y2") +
  theme_minimal()


# 3D density plot for sigma = 0.1

df_s <- sim_df[sim_df$sigma == "0.1", ]

kd3 <- kde2d(df_s$Y1, df_s$Y2, n = 150,
             lims = c(global_x[1], global_x[2],
                      global_y[1], global_y[2]))

persp(kd3$x, kd3$y, kd3$z,
      theta = 40, phi = 35,
      col = "lightblue",
      shade = 0.6,
      xlab = "Y1", ylab = "Y2", zlab = "Density",
      main = "3D Density Surface (σ = 0.1)")