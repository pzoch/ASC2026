rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dynlm)
  library(forecast)
  library(ggplot2)
  library(lmtest)
  library(lubridate)
  library(urca)
  library(zoo)
})

# Locate the project root so the script works from any working directory.
find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "data")) &&
        dir.exists(file.path(current, "code", "unit_root"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Project root not found. Expected folders: data/ and code/unit_root/.")
    }
    current <- parent
  }
}

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

plot_series <- function(series, title, y_label = "Value", color = "#1f77b4") {
  df <- data.frame(
    date = as.POSIXct(zoo::index(series)),
    value = as.numeric(series)
  )

  ggplot(df, aes(x = date, y = value)) +
    geom_line(linewidth = 0.9, color = color) +
    labs(title = title, x = NULL, y = y_label) +
    plot_theme
}

plot_correlation <- function(series, label, lag_max = 12) {
  acf_plot <- ggAcf(as.numeric(series), lag.max = lag_max) +
    labs(title = paste(label, "- ACF"), x = "Lag", y = "ACF") +
    plot_theme

  pacf_plot <- ggPacf(as.numeric(series), lag.max = lag_max) +
    labs(title = paste(label, "- PACF"), x = "Lag", y = "PACF") +
    plot_theme

  print(acf_plot)
  print(pacf_plot)
}

run_bg_loop <- function(model, orders) {
  bg_results <- lapply(orders, function(ord) bgtest(model, order = ord))
  bg_summary <- data.frame(
    order = orders,
    statistic = sapply(bg_results, function(res) as.numeric(res$statistic)),
    p_value = sapply(bg_results, function(res) res$p.value)
  )

  cat("\nBreusch-Godfrey LM tests (ADF residual autocorrelation)\n")
  print(bg_summary, row.names = FALSE)
  invisible(bg_results)
}

# Data ---------------------------------------------------------------------
project_root <- find_project_root()
gdp_path <- file.path(project_root, "data", "GDP_POLAND.csv")
if (!file.exists(gdp_path)) {
  stop("File not found: ", gdp_path)
}

data_pl <- fread(gdp_path)
setnames(data_pl, c("date", "Y"))
data_pl[, date := as.POSIXct(strptime(date, "%m/%d/%Y", tz = "UTC"))]
data_pl <- data_pl[order(date)]

Y <- zoo(data_pl$Y, order.by = data_pl$date)
sample_end <- parse_date_time("01/01/2015", orders = "mdy", tz = "UTC")
Y_short <- window(Y, end = sample_end)

# Plots: levels, logs, differences ----------------------------------------
print(plot_series(Y, "Poland GDP level", "GDP"))
print(plot_series(Y_short, "Poland GDP level (sample ending 2015-01-01)", "GDP"))

y <- log(Y)
y_short <- log(Y_short)
print(plot_series(y, "Log GDP", "log(GDP)", "#2a9d8f"))
print(plot_series(y_short, "Log GDP (sample ending 2015-01-01)", "log(GDP)", "#2a9d8f"))

dy <- diff(y)
dy_short <- diff(y_short)
print(plot_series(dy, "Quarterly log difference", "Delta log(GDP)", "#e76f51"))
print(plot_series(dy_short, "Quarterly log difference (sample ending 2015-01-01)", "Delta log(GDP)", "#e76f51"))

# Deterministic detrending --------------------------------------------------
trend <- seq_along(y)
trend_sq <- trend ^ 2

reg_detr_linear <- lm(as.numeric(y) ~ trend)
y_detr_linear <- zoo(resid(reg_detr_linear), order.by = index(y))
print(plot_series(y_detr_linear, "Detrended log GDP (linear trend)", "Residual"))

reg_detr_quadratic <- lm(as.numeric(y) ~ trend + trend_sq)
y_detr_quadratic <- zoo(resid(reg_detr_quadratic), order.by = index(y))
print(plot_series(y_detr_quadratic, "Detrended log GDP (quadratic trend)", "Residual"))

# Correlation diagnostics ---------------------------------------------------
plot_correlation(y, "Log GDP")
plot_correlation(y_short, "Log GDP (short sample)")
plot_correlation(dy, "Delta log GDP")
plot_correlation(dy_short, "Delta log GDP (short sample)")

# ADF tests ----------------------------------------------------------------
cat("\nADF tests\n")
print(summary(ur.df(as.numeric(y), type = "none", lags = 4)))
print(summary(ur.df(as.numeric(y), type = "drift", lags = 4)))
print(summary(ur.df(as.numeric(y), type = "trend", lags = 4)))
print(summary(ur.df(as.numeric(dy), type = "drift", lags = 4)))
print(summary(ur.df(as.numeric(y_short), type = "none", lags = 4)))
print(summary(ur.df(as.numeric(y_short), type = "drift", lags = 4)))
print(summary(ur.df(as.numeric(y_short), type = "trend", lags = 4)))
print(summary(ur.df(as.numeric(dy_short), type = "drift", lags = 4)))

# BG test for the no-drift DF regression -----------------------------------
reg_df <- dynlm(dy ~ L(y, 1) - 1)
print(summary(reg_df))

df_resid <- residuals(reg_df)
print(
  ggAcf(as.numeric(df_resid), lag.max = 12) +
    labs(title = "DF regression residuals - ACF", x = "Lag", y = "ACF") +
    plot_theme
)
print(
  ggPacf(as.numeric(df_resid), lag.max = 12) +
    labs(title = "DF regression residuals - PACF", x = "Lag", y = "PACF") +
    plot_theme
)

run_bg_loop(reg_df, orders = 1:4)

# KPSS tests ---------------------------------------------------------------
cat("\nKPSS tests\n")
print(summary(ur.kpss(as.numeric(y), type = "tau", lags = "short")))
print(summary(ur.kpss(as.numeric(dy), type = "mu", lags = "short")))
print(summary(ur.kpss(as.numeric(y_short), type = "tau", lags = "short")))
print(summary(ur.kpss(as.numeric(dy_short), type = "mu", lags = "short")))
