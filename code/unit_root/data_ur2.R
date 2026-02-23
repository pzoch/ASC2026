rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dynlm)
  library(forecast)
  library(ggplot2)
  library(lmtest)
  library(urca)
  library(zoo)
})

# Locate the project root so this script can be run from any folder.
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

project_root <- find_project_root()
agg_path <- file.path(project_root, "data", "agg_data.csv")
if (!file.exists(agg_path)) {
  stop("File not found: ", agg_path)
}

# Data prep ----------------------------------------------------------------
agg_data <- fread(agg_path)
forward_rate_cols <- c(paste0("TIPSF0", 2:9), paste0("TIPSF", 10:20))
hqm_spread_cols <- paste0("hqm_TD_spread", c(2, 5, 10, 20))

agg_data[, leverage := liquid_effective_leverage_agg]
keep_cols <- c(
  "date", "real_rate_1year", forward_rate_cols,
  "CP_TB_spread", hqm_spread_cols, "leverage"
)
agg_data <- na.omit(agg_data[, ..keep_cols])
agg_data[, date := as.Date(date)]
setorder(agg_data, date)

y <- agg_data$leverage
y_zoo <- zoo(y, order.by = agg_data$date)
dy_zoo <- diff(y_zoo)

# Plot leverage -------------------------------------------------------------
lev_df <- data.frame(
  date = agg_data$date,
  leverage = y
)

print(
  ggplot(lev_df, aes(x = date, y = leverage)) +
    geom_line(linewidth = 0.9, color = "#264653") +
    labs(
      title = "Aggregate effective leverage",
      x = NULL,
      y = "Leverage"
    ) +
    plot_theme
)

print(
  ggAcf(y, lag.max = 12) +
    labs(title = "Leverage - ACF", x = "Lag", y = "ACF") +
    plot_theme
)
print(
  ggPacf(y, lag.max = 12) +
    labs(title = "Leverage - PACF", x = "Lag", y = "PACF") +
    plot_theme
)

# ADF tests ----------------------------------------------------------------
cat("\nADF tests\n")
print(summary(ur.df(y, type = "none", lags = 4)))
print(summary(ur.df(y, type = "drift", lags = 5)))
print(summary(ur.df(y, type = "trend", lags = 5)))
print(summary(ur.df(as.numeric(dy_zoo), type = "drift", lags = 4)))

# BG test for the trend-ADF specification ----------------------------------
trend_zoo <- zoo(seq_along(y), order.by = agg_data$date)
reg_df <- dynlm(
  dy_zoo ~ L(y_zoo, 1) + trend_zoo +
    L(dy_zoo, 1) + L(dy_zoo, 2) + L(dy_zoo, 3) + L(dy_zoo, 4) + L(dy_zoo, 5)
)
print(summary(reg_df))

reg_resid <- residuals(reg_df)
print(
  ggAcf(as.numeric(reg_resid), lag.max = 12) +
    labs(title = "Trend-ADF residuals - ACF", x = "Lag", y = "ACF") +
    plot_theme
)
print(
  ggPacf(as.numeric(reg_resid), lag.max = 12) +
    labs(title = "Trend-ADF residuals - PACF", x = "Lag", y = "PACF") +
    plot_theme
)

run_bg_loop(reg_df, orders = c(1, 2, 3, 4, 12))

# KPSS tests ---------------------------------------------------------------
cat("\nKPSS tests\n")
print(summary(ur.kpss(y, type = "tau", lags = "short")))
print(summary(ur.kpss(as.numeric(dy_zoo), type = "mu", lags = "short")))
