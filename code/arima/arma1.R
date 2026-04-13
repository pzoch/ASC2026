rm(list = ls())
set.seed(13)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forecast)
  library(lmtest)
  library(tseries)
})

plot_theme <- theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot", panel.grid.minor = element_blank())

sim_proc <- function(ar = numeric(0), ma = numeric(0), n = 1000, sd = 1) {
  arima.sim(model = list(ar = ar, ma = ma), n = n, sd = sd)
}

to_long <- function(series_list) {
  rbindlist(lapply(names(series_list), function(nm) {
    data.table(t = seq_along(series_list[[nm]]),
               y = as.numeric(series_list[[nm]]),
               series = nm)
  }))
}

coef_ar <- function(m) {
  p <- arimaorder(m)["p"]
  if (p > 0) as.numeric(coef(m)[paste0("ar", seq_len(p))]) else numeric(0)
}

coef_ma <- function(m) {
  q <- arimaorder(m)["q"]
  if (q > 0) as.numeric(coef(m)[paste0("ma", seq_len(q))]) else numeric(0)
}

irf_from <- function(m, h = 24) {
  psi <- ARMAtoMA(ar = coef_ar(m), ma = coef_ma(m), lag.max = h)
  ts(c(1, psi), start = 0)
}

# MA(4) processes -----------------------------------------------------------
ma_sets <- list(
  "MA vec 1" = c(1, 1, 1, 1),
  "MA vec 2" = c(-1, 2, -3, 4),
  "MA vec 3" = c(2, 1, 0.5, 0.2)
)
ma_series <- lapply(ma_sets, function(v) sim_proc(ma = v, n = 1000))

print(
  ggplot(to_long(ma_series), aes(x = t, y = y, color = series)) +
    geom_line(linewidth = 0.4) +
    labs(title = "Simulated MA(4) processes", x = NULL, y = NULL, color = NULL) +
    plot_theme
)

for (nm in names(ma_series)) {
  print(ggAcf(ma_series[[nm]], lag.max = 12) +
          labs(title = paste("ACF —", nm)) + plot_theme)
  print(ggPacf(ma_series[[nm]], lag.max = 12) +
          labs(title = paste("PACF —", nm)) + plot_theme)
}

# AR(4) processes -----------------------------------------------------------
ar_sets <- list(
  "AR vec 1" = c(0.5, 0, 0, 0),
  "AR vec 2" = c(0.4, 0.2, 0.1, 0.05),
  "AR vec 3" = c(-0.3, 0.2, -0.1, 0.05)
)
ar_series <- lapply(ar_sets, function(v) sim_proc(ar = v, n = 1000))

print(
  ggplot(to_long(ar_series), aes(x = t, y = y, color = series)) +
    geom_line(linewidth = 0.4) +
    labs(title = "Simulated AR(4) processes", x = NULL, y = NULL, color = NULL) +
    plot_theme
)

for (nm in names(ar_series)) {
  print(ggAcf(ar_series[[nm]], lag.max = 12) +
          labs(title = paste("ACF —", nm)) + plot_theme)
  print(ggPacf(ar_series[[nm]], lag.max = 12) +
          labs(title = paste("PACF —", nm)) + plot_theme)
}

# ARMA example: fit and compare specifications -----------------------------
ar_vec <- c(0.9, -0.3)
ma_vec <- c(0.3, 3, 0.1)
y_arma <- sim_proc(ar = ar_vec, ma = ma_vec, n = 100)

print(
  autoplot(y_arma) +
    labs(title = "Simulated ARMA(2,3) process", x = NULL, y = NULL) +
    plot_theme
)
print(ggAcf(y_arma,  lag.max = 12) + labs(title = "ACF")  + plot_theme)
print(ggPacf(y_arma, lag.max = 12) + labs(title = "PACF") + plot_theme)

model_23 <- Arima(y_arma, order = c(2, 0, 3), method = "ML")
model_32 <- Arima(y_arma, order = c(3, 0, 2), method = "ML")

cat("\nARMA(2,3) coefficient tests:\n")
print(coeftest(model_23))
cat("\nARMA(3,2) coefficient tests:\n")
print(coeftest(model_32))

# Residual diagnostics ------------------------------------------------------
resids <- residuals(model_23)
print(ggAcf(resids,  lag.max = 12) + labs(title = "Residual ACF")  + plot_theme)
print(ggPacf(resids, lag.max = 12) + labs(title = "Residual PACF") + plot_theme)

cat("\nLjung-Box on residuals (ARMA(2,3)):\n")
print(Box.test(resids, lag = 12, type = "Ljung-Box", fitdf = 2 + 3))
cat("\nJarque-Bera normality test on residuals:\n")
print(jarque.bera.test(resids))
checkresiduals(model_23)

# Information criteria search -----------------------------------------------
model_aic <- auto.arima(y_arma, max.d = 0, max.p = 4, max.q = 4, ic = "aic",
                        seasonal = FALSE, stepwise = FALSE, approximation = FALSE)
model_bic <- auto.arima(y_arma, max.d = 0, max.p = 4, max.q = 4, ic = "bic",
                        seasonal = FALSE, stepwise = FALSE, approximation = FALSE)
cat("\nauto.arima (AIC) selected order:", arimaorder(model_aic), "\n")
cat("auto.arima (BIC) selected order:", arimaorder(model_bic), "\n")

# Impulse response: compare estimated vs true -------------------------------
h_max <- 24
irf_df <- rbind(
  data.table(h = 0:h_max, psi = as.numeric(ts(c(1, ARMAtoMA(ar = ar_vec, ma = ma_vec, lag.max = h_max)))), series = "True"),
  data.table(h = 0:h_max, psi = as.numeric(irf_from(model_aic, h_max)), series = "AIC"),
  data.table(h = 0:h_max, psi = as.numeric(irf_from(model_bic, h_max)), series = "BIC"),
  data.table(h = 0:h_max, psi = as.numeric(irf_from(model_23,  h_max)), series = "ARMA(2,3)")
)

print(
  ggplot(irf_df, aes(x = h, y = psi, color = series)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    labs(title = "Impulse response functions",
         x = "h", y = expression(psi[h]), color = NULL) +
    plot_theme
)
