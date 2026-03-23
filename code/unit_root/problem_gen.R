rm(list = ls())

suppressPackageStartupMessages({
  library(dynlm)
  library(forecast)
  library(ggplot2)
  library(lmtest)
  library(urca)
  library(zoo)
})

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

plot_series <- function(x, title, y_label = "Value", color = "#1f77b4") {
  df <- data.frame(t = seq_along(x), value = as.numeric(x))
  ggplot(df, aes(x = t, y = value)) +
    geom_line(linewidth = 0.9, color = color) +
    labs(title = title, x = "t", y = y_label) +
    plot_theme
}

run_bg_loop <- function(model, orders) {
  bg_results <- lapply(orders, function(ord) bgtest(model, order = ord))
  bg_summary <- data.frame(
    order = orders,
    statistic = sapply(bg_results, function(res) as.numeric(res$statistic)),
    p_value = sapply(bg_results, function(res) res$p.value)
  )

  cat("\nBreusch-Godfrey LM tests\n")
  print(bg_summary, row.names = FALSE)
  invisible(bg_results)
}

adf_bg_table <- function(y_zoo, type = "none", max_lag = 6, bg_order = 4) {
  dy <- diff(y_zoo)
  results <- lapply(0:max_lag, function(p) {
    if (p == 0) {
      lag_part <- ""
    } else {
      lag_part <- paste("+", paste0("L(dy, ", 1:p, ")", collapse = " + "))
    }
    fml_str <- switch(type,
      none  = paste("dy ~ L(y_zoo, 1)", lag_part, "- 1"),
      drift = paste("dy ~ L(y_zoo, 1)", lag_part),
      trend = paste("dy ~ L(y_zoo, 1) + trend(dy)", lag_part)
    )
    reg <- dynlm(as.formula(fml_str))
    bg <- bgtest(reg, order = bg_order)
    data.frame(
      adf_lags = p,
      BG_stat = round(as.numeric(bg$statistic), 3),
      BG_pval = round(bg$p.value, 4)
    )
  })
  do.call(rbind, results)
}

simulate_ar2 <- function(phi, n, sigma2 = 1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  x <- numeric(n)
  x[1] <- rnorm(1, sd = sqrt(sigma2))
  x[2] <- phi[1] * x[1] + rnorm(1, sd = sqrt(sigma2))
  for (t in 3:n) {
    x[t] <- phi[1] * x[t - 1] + phi[2] * x[t - 2] + rnorm(1, sd = sqrt(sigma2))
  }
  x
}

# Near-unit-root AR(2) example ----------------------------------------------
phi <- c(0.99, -0.1)
sigma2 <- 1
nobs <- 100
x <- simulate_ar2(phi = phi, n = nobs, sigma2 = sigma2, seed = 220)
dx <- diff(x)

print(plot_series(x, "Near-unit-root AR(2) simulation", "x_t"))
print(plot_series(dx, "First difference", "Delta x_t", "#e76f51"))

print(ggAcf(x, lag.max = 24) + labs(title = "Near-unit-root series - ACF") + plot_theme)
print(ggPacf(x, lag.max = 24) + labs(title = "Near-unit-root series - PACF") + plot_theme)

# DF tests ------------------------------------------------------------------
df_test1 <- ur.df(x, type = "none", lags = 0)
df_test2 <- ur.df(dx, type = "none", lags = 0)
df_test3 <- ur.df(x, type = "drift", lags = 0)
df_test4 <- ur.df(x, type = "trend", lags = 0)

print(plot_series(df_test1@res, "DF(no deterministic terms) residuals", "Residual", "#264653"))
print(ggAcf(df_test1@res, lag.max = 24) + labs(title = "DF residuals - ACF") + plot_theme)
print(ggPacf(df_test1@res, lag.max = 24) + labs(title = "DF residuals - PACF") + plot_theme)

cat("\nDF test summaries\n")
print(summary(df_test1))
print(summary(df_test2))
print(summary(df_test3))
print(summary(df_test4))

# BG tests across ADF lag lengths -------------------------------------------
x_zoo <- zoo(x)
cat("\nBG test (order 4) for ADF(none) with varying lag lengths\n")
print(adf_bg_table(x_zoo, type = "none", max_lag = 6, bg_order = 4))

# KPSS tests ----------------------------------------------------------------
kpss_test1 <- ur.kpss(x, type = "mu")
kpss_test2 <- ur.kpss(x, type = "tau")

cat("\nKPSS test summaries\n")
print(summary(kpss_test1))
print(summary(kpss_test2))

# Monte Carlo: rejection rates ----------------------------------------------
# DGP: y_t = (1-rho)*alpha + rho*y_{t-1} + e_t, e_t ~ N(0, sigma2)
set.seed(777)
T <- 4 * 25
N <- 100
rho <- 0.5
sigma2 <- 1.0
alpha <- 1.0

count <- matrix(NA_real_, nrow = N, ncol = 2)
for (n in 1:N) {
  e <- rnorm(T, sd = sqrt(sigma2))
  y <- numeric(T)

  y[1] <- alpha + e[1] / sqrt(1 - rho ^ 2)
  for (t in 2:T) {
    y[t] <- (1 - rho) * alpha + rho * y[t - 1] + e[t]
  }

  adf <- ur.df(y, type = "drift", lags = 0)
  count[n, 1] <- adf@teststat[1] < adf@cval[1, 2]

  kpss <- ur.kpss(y, type = "mu")
  count[n, 2] <- kpss@teststat < kpss@cval[2]
}

rejection_rates <- colMeans(count) * 100
names(rejection_rates) <- c("ADF", "KPSS")

cat("\nEstimated rejection rates (%)\n")
print(rejection_rates)
