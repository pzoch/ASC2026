rm(list = ls())
set.seed(13)

suppressPackageStartupMessages({
  library(forecast)
  library(ggplot2)
  library(urca)
  library(uroot)
})

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

plot_ts <- function(x, title, y_label = "Value", color = "#1f77b4") {
  df <- data.frame(time = as.numeric(time(x)), value = as.numeric(x))
  ggplot(df, aes(x = time, y = value)) +
    geom_line(linewidth = 0.9, color = color) +
    labs(title = title, x = "Time", y = y_label) +
    plot_theme
}

plot_periodogram <- function(x, title, span = NULL) {
  spec <- spec.pgram(as.numeric(x), log = "no", span = span, plot = FALSE)
  df <- data.frame(freq = spec$freq, spec = spec$spec)
  ggplot(df, aes(x = freq, y = spec)) +
    geom_line(linewidth = 0.8, color = "#264653") +
    labs(title = title, x = "Frequency", y = "Spectral density") +
    plot_theme
}

plot_corr <- function(x, label, lag_max = 24) {
  print(
    ggAcf(as.numeric(x), lag.max = lag_max) +
      labs(title = paste(label, "- ACF"), x = "Lag", y = "ACF") +
      plot_theme
  )
  print(
    ggPacf(as.numeric(x), lag.max = lag_max) +
      labs(title = paste(label, "- PACF"), x = "Lag", y = "PACF") +
      plot_theme
  )
}

# Simulated quarterly seasonality -------------------------------------------
no_periods <- 20 * 4
noise_scale <- 0.1

white_noise <- arima.sim(model = list(order = c(0, 0, 0)), n = no_periods, sd = noise_scale)
x_det <- ts(rep(c(-1, 0, 1, 2), no_periods / 4) + white_noise, start = c(1980, 1), frequency = 4)

model_ar <- Arima(
  ts(rnorm(no_periods, sd = noise_scale), frequency = 4),
  order = c(0, 0, 0),
  seasonal = list(order = c(1, 0, 0), period = 4),
  fixed = c(0.5, 0)
)
x_ar <- ts(simulate(model_ar, nsim = no_periods), start = c(1980, 1), frequency = 4)

model_ur <- Arima(
  ts(rnorm(no_periods, sd = noise_scale), frequency = 4),
  order = c(0, 0, 0),
  seasonal = list(order = c(0, 1, 0), period = 4)
)
x_ur <- ts(simulate(model_ur, nsim = no_periods), start = c(1980, 1), frequency = 4)

horizon <- 12
quarters <- rep(1:4, length.out = no_periods + horizon)
dummies <- ts(
  cbind(
    as.integer(quarters == 1),
    as.integer(quarters == 2),
    as.integer(quarters == 3)
  ),
  start = c(1980, 1),
  frequency = 4
)
dummies_reg <- window(dummies, end = end(x_det))
dummies_fcast <- window(dummies, start = c(2000, 1))

model_det <- Arima(x_det, order = c(0, 0, 0), xreg = dummies_reg, method = "ML")
det_fcst <- forecast(model_det, h = NROW(dummies_fcast), xreg = dummies_fcast)
print(
  autoplot(det_fcst) +
    labs(title = "Deterministic seasonality: forecast with quarterly dummies") +
    plot_theme
)

model_ur_fit <- Arima(
  x_ur,
  order = c(0, 0, 0),
  seasonal = list(order = c(0, 1, 0), period = 4)
)
ur_fcst <- forecast(model_ur_fit, h = horizon)
print(
  autoplot(ur_fcst) +
    labs(title = "Seasonal unit root: forecast from ARIMA(0,0,0)(0,1,0)[4]") +
    plot_theme
)

sim_df <- rbind(
  data.frame(time = as.numeric(time(x_det)), value = as.numeric(x_det), series = "Deterministic"),
  data.frame(time = as.numeric(time(x_ar)), value = as.numeric(x_ar), series = "Seasonal AR"),
  data.frame(time = as.numeric(time(x_ur)), value = as.numeric(x_ur), series = "Seasonal unit root")
)

print(
  ggplot(sim_df, aes(x = time, y = value, color = series)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c("Deterministic" = "#1f77b4", "Seasonal AR" = "#2a9d8f", "Seasonal unit root" = "#e76f51")) +
    labs(title = "Three types of quarterly seasonality", x = "Time", y = "Value", color = "Process") +
    plot_theme
)

sim_df_ar_ur <- subset(sim_df, series %in% c("Seasonal AR", "Seasonal unit root"))
print(
  ggplot(sim_df_ar_ur, aes(x = time, y = value, color = series)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c("Seasonal AR" = "#2a9d8f", "Seasonal unit root" = "#e76f51")) +
    labs(title = "Seasonal AR vs seasonal unit root", x = "Time", y = "Value", color = "Process") +
    plot_theme
)

print(plot_periodogram(x_det, "Deterministic seasonality - periodogram"))
print(plot_periodogram(x_ar, "Seasonal AR - periodogram"))
print(plot_periodogram(x_ur, "Seasonal unit root - periodogram"))

plot_corr(x_det, "Deterministic seasonality")
plot_corr(x_ar, "Seasonal AR")
plot_corr(x_ur, "Seasonal unit root")

# AirPassengers: Box-Jenkins workflow ---------------------------------------
print(plot_ts(AirPassengers, "AirPassengers", "Passengers"))

selected_years <- c(1951, 1953, 1955)
ap_year_df <- do.call(
  rbind,
  lapply(selected_years, function(year) {
    year_window <- window(AirPassengers, start = c(year, 1), end = c(year, 12))
    data.frame(year = factor(year), month = c(cycle(year_window)), value = c(year_window))
  })
)

print(
  ggplot(ap_year_df, aes(x = month, y = value)) +
    geom_line(linewidth = 0.9, color = "#264653") +
    geom_point(size = 1.5, color = "#264653") +
    facet_wrap(~ year, nrow = 1) +
    labs(title = "AirPassengers seasonality by selected years", x = "Month", y = "Passengers") +
    plot_theme
)

lAP <- log(AirPassengers)
print(plot_ts(lAP, "Log AirPassengers", "log(passengers)", "#2a9d8f"))

lap_year_df <- do.call(
  rbind,
  lapply(selected_years, function(year) {
    year_window <- window(lAP, start = c(year, 1), end = c(year, 12))
    data.frame(year = factor(year), month = c(cycle(year_window)), value = c(year_window))
  })
)

print(
  ggplot(lap_year_df, aes(x = month, y = value)) +
    geom_line(linewidth = 0.9, color = "#2a9d8f") +
    geom_point(size = 1.5, color = "#2a9d8f") +
    facet_wrap(~ year, nrow = 1) +
    labs(title = "Log AirPassengers seasonality by selected years", x = "Month", y = "log(passengers)") +
    plot_theme
)

plot_corr(lAP, "Log AirPassengers")
print(plot_periodogram(lAP, "Log AirPassengers - periodogram", span = 5))

df_test_none <- ur.df(lAP, type = "none", lags = 0)
df_test_trend <- ur.df(lAP, type = "trend", lags = 12)
cat("\nADF tests for log AirPassengers\n")
print(summary(df_test_none))
print(summary(df_test_trend))

dlAP12 <- diff(lAP, lag = 12)
print(plot_ts(dlAP12, "Seasonal difference: Delta_12 log(AirPassengers)", "Value", "#e76f51"))
plot_corr(dlAP12, "Seasonal difference")
print(plot_periodogram(dlAP12, "Seasonal difference - periodogram", span = 5))
print(summary(ur.df(dlAP12, type = "none", lags = 12)))

HEGY_lAP <- hegy.test(lAP, deterministic = c(1, 1, 1), lag.method = "AIC", maxlag = 24)
CH_lAP <- ch.test(lAP, type = "trigonometric", sid = "all", lag1 = TRUE)
print(HEGY_lAP)
print(CH_lAP)

dlAP <- diff(lAP)
print(plot_ts(dlAP, "First difference: Delta log(AirPassengers)", "Value", "#f4a261"))
plot_corr(dlAP, "First difference")
print(plot_periodogram(dlAP, "First difference - periodogram", span = 5))

HEGY_dlAP <- hegy.test(dlAP, deterministic = c(1, 1, 1), lag.method = "AIC", maxlag = 24)
CH_dlAP <- ch.test(dlAP, type = "trigonometric", sid = "all", lag1 = TRUE)
print(HEGY_dlAP)
print(CH_dlAP)

print(summary(ur.df(dlAP, type = "drift", lags = 0)))
print(summary(ur.df(dlAP, type = "trend", lags = 0)))

dSlAP <- diff(lAP, lag = 12)
print(plot_ts(dSlAP, "Seasonal difference dS(log AP)", "Value", "#e76f51"))
plot_corr(dSlAP, "Seasonal difference dS(log AP)")
print(plot_periodogram(dSlAP, "Seasonal difference dS(log AP) - periodogram", span = 5))

HEGY_dSlAP <- hegy.test(dSlAP, deterministic = c(1, 1, 1), lag.method = "AIC", maxlag = 24)
CH_dSlAP <- ch.test(dSlAP, type = "trigonometric", sid = "all", lag1 = TRUE)
print(HEGY_dSlAP)
print(CH_dSlAP)

print(summary(ur.df(dSlAP, type = "trend", lags = 12)))

ddSlAP <- diff(dlAP, lag = 12)
print(plot_ts(ddSlAP, "Combined difference dS(d log AP)", "Value", "#f4a261"))
plot_corr(ddSlAP, "Combined difference dS(d log AP)")
print(plot_periodogram(ddSlAP, "Combined difference dS(d log AP) - periodogram", span = 5))

HEGY_ddSlAP <- hegy.test(ddSlAP, deterministic = c(1, 1, 1), lag.method = "AIC", maxlag = 24)
CH_ddSlAP <- ch.test(ddSlAP, type = "trigonometric", sid = "all", lag1 = TRUE)
print(HEGY_ddSlAP)
print(CH_ddSlAP)

model_s <- arima(
  lAP,
  order = c(1, 1, 1),
  seasonal = list(order = c(1, 1, 1), period = 12),
  method = "ML"
)
res_model_s <- residuals(model_s)
print(Box.test(res_model_s, lag = 1, type = "Ljung"))
print(Box.test(res_model_s, lag = 12, type = "Ljung"))
checkresiduals(model_s)

model_auto <- auto.arima(lAP, method = "ML", max.d = 1, max.p = 24, max.q = 24, ic = "aic")
checkresiduals(model_auto)
print(
  autoplot(forecast(model_auto)) +
    labs(title = "Forecast from auto.arima(log AirPassengers)") +
    plot_theme
)
