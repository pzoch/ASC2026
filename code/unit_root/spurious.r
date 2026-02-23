rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
})

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

# Random walk with drift:
# y_t = mu + y_(t-1) + e_t
rho <- 1
sigma2 <- 1
mu <- c(0, 0)
nobs <- 100

simulate_rw_pair <- function(n, rho, mu, sigma2) {
  rw <- matrix(0, nrow = n, ncol = 2)
  for (k in 1:2) {
    for (t in 2:n) {
      rw[t, k] <- mu[k] + rho * rw[t - 1, k] + rnorm(1, sd = sqrt(sigma2))
    }
  }
  rw
}

set.seed(234)
rw_sim <- simulate_rw_pair(n = nobs, rho = rho, mu = mu, sigma2 = sigma2)
rw_df <- data.frame(
  obs = seq_len(nobs),
  x1 = rw_sim[, 1],
  x2 = rw_sim[, 2]
)

# Time series plot -----------------------------------------------------------
rw_long <- rbind(
  data.frame(obs = rw_df$obs, value = rw_df$x1, series = "x1"),
  data.frame(obs = rw_df$obs, value = rw_df$x2, series = "x2")
)

print(
  ggplot(rw_long, aes(x = obs, y = value, color = series)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c("x1" = "#1f77b4", "x2" = "#e76f51")) +
    labs(
      title = "Random walk simulation",
      x = "Observation",
      y = NULL,
      color = "Series"
    ) +
    plot_theme
)

# Spurious regression --------------------------------------------------------
print(
  ggplot(rw_df, aes(x = x1, y = x2)) +
    geom_point(alpha = 0.65, color = "#264653") +
    geom_smooth(method = "lm", se = FALSE, color = "#e76f51", linewidth = 0.9) +
    labs(
      title = "Scatter plot of two independent random walks",
      x = "x1",
      y = "x2"
    ) +
    plot_theme
)

spurious_model <- lm(x1 ~ 0 + x2, data = rw_df)
print(summary(spurious_model))

rw_df$residual <- resid(spurious_model)
print(
  ggplot(rw_df, aes(x = obs, y = residual)) +
    geom_line(linewidth = 0.9, color = "#2a9d8f") +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", color = "gray40") +
    labs(title = "Residuals from OLS fit", x = "Observation", y = "Residual") +
    plot_theme
)

# Monte Carlo distribution of spurious t-statistics -------------------------
set.seed(234)
nexper <- 5000
t_values <- numeric(nexper)

for (n in 1:nexper) {
  rw_sim <- simulate_rw_pair(n = nobs, rho = rho, mu = mu, sigma2 = sigma2)
  rw_df <- data.frame(x1 = rw_sim[, 1], x2 = rw_sim[, 2])
  spurious_model <- lm(x1 ~ 0 + x2, data = rw_df)
  t_values[n] <- coef(summary(spurious_model))[1, "t value"]
}

t_df <- data.frame(tstat = t_values)
print(
  ggplot(t_df, aes(x = tstat)) +
    geom_density(linewidth = 1, fill = "#1f77b4", alpha = 0.25, color = "#1f77b4") +
    stat_function(
      fun = dnorm,
      args = list(mean = 0, sd = 1),
      color = "#e76f51",
      linetype = "dotted",
      linewidth = 1
    ) +
    geom_vline(xintercept = mean(t_values), color = "#2a9d8f", linetype = "dashed", linewidth = 0.9) +
    geom_vline(xintercept = c(-1.96, 1.96), color = "gray30", linetype = "dashed", linewidth = 0.7) +
    coord_cartesian(xlim = c(-10, 10)) +
    labs(
      title = "Distribution of t-statistics from spurious regressions",
      x = "t-statistic",
      y = "Density"
    ) +
    plot_theme
)
