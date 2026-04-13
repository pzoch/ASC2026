rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forecast)
  library(lmtest)
  library(urca)
})

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "data")) &&
        dir.exists(file.path(current, "code"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) stop("Project root not found.")
    current <- parent
  }
}

plot_theme <- theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot", panel.grid.minor = element_blank())

project_root <- find_project_root()

# Data: US CPI monthly inflation rate ---------------------------------------
cpi <- fread(file.path(project_root, "data", "CPI_US_NSA.csv"))
setnames(cpi, c("date", "pi"))
cpi[, date := as.Date(date, format = "%m/%d/%Y")]

pi_ts <- ts(cpi$pi, start = c(1960, 1), frequency = 12)
pi_ts <- window(pi_ts, start = c(2002, 1), end = c(2019, 12))

train <- window(pi_ts, end = c(2018, 12))
test  <- window(pi_ts, start = c(2019, 1))

print(
  autoplot(pi_ts) +
    labs(title = "US CPI monthly inflation rate (%)", x = NULL, y = "%") +
    plot_theme
)

# Box-Jenkins identification ------------------------------------------------
print(ggAcf(train,  lag.max = 24) + labs(title = "ACF")  + plot_theme)
print(ggPacf(train, lag.max = 24) + labs(title = "PACF") + plot_theme)

cat("\nADF test on levels (type = none, lags = 4):\n")
print(summary(ur.df(train, type = "none", lags = 4)))

dtrain <- diff(train)
print(ggAcf(dtrain,  lag.max = 24) + labs(title = "ACF of first difference")  + plot_theme)
print(ggPacf(dtrain, lag.max = 24) + labs(title = "PACF of first difference") + plot_theme)

# Initial model and alternatives --------------------------------------------
model         <- Arima(train, order = c(1, 1, 2), method = "ML")
model_plus_ar <- Arima(train, order = c(2, 1, 2), method = "ML")
model_plus_ma <- Arima(train, order = c(1, 1, 3), method = "ML")

cat("\nARIMA(1,1,2):\n");  print(coeftest(model))
cat("\nARIMA(2,1,2):\n");  print(coeftest(model_plus_ar))
cat("\nARIMA(1,1,3):\n");  print(coeftest(model_plus_ma))

cat("\nLjung-Box on ARIMA(1,1,2) residuals:\n")
print(Box.test(residuals(model), lag = 24, type = "Ljung-Box", fitdf = 1 + 2))
checkresiduals(model)

# Information criteria ------------------------------------------------------
model_aic <- auto.arima(train, max.d = 1, max.p = 12, max.q = 12, ic = "aic",
                        seasonal = FALSE, stepwise = FALSE, approximation = FALSE)
model_bic <- auto.arima(train, max.d = 1, max.p = 12, max.q = 12, ic = "bic",
                        seasonal = FALSE, stepwise = FALSE, approximation = FALSE)

cat("\nauto.arima (AIC) selected order:", arimaorder(model_aic), "\n")
cat("auto.arima (BIC) selected order:", arimaorder(model_bic), "\n")
checkresiduals(model_aic)
checkresiduals(model_bic)

# General-to-specific -------------------------------------------------------
model_general  <- Arima(train, order = c(1, 1, 2), method = "ML")
model_specific <- Arima(train, order = c(1, 1, 0), method = "ML")

lr_stat <- -2 * as.numeric(logLik(model_specific) - logLik(model_general))
df_diff <- length(coef(model_general)) - length(coef(model_specific))
p_value <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)

cat("\nLR test (general vs specific):\n")
cat("  statistic =", round(lr_stat, 3),
    "  df =", df_diff,
    "  p-value =", signif(p_value, 4), "\n")

# Forecasts -----------------------------------------------------------------
h <- length(test)

print(
  autoplot(forecast(model, h = h)) +
    autolayer(pi_ts, series = "Actual") +
    labs(title = "ARIMA(1,1,2) forecast", x = NULL, y = "%", colour = NULL) +
    plot_theme
)
print(
  autoplot(forecast(model_aic, h = h)) +
    autolayer(pi_ts, series = "Actual") +
    labs(title = "auto.arima (AIC) forecast", x = NULL, y = "%", colour = NULL) +
    plot_theme
)
print(
  autoplot(forecast(model_bic, h = h)) +
    autolayer(pi_ts, series = "Actual") +
    labs(title = "auto.arima (BIC) forecast", x = NULL, y = "%", colour = NULL) +
    plot_theme
)

# Seasonal ARIMA ------------------------------------------------------------
dtrain12 <- diff(train, lag = 12)
print(ggAcf(dtrain12,  lag.max = 24) + labs(title = "ACF of seasonal difference")  + plot_theme)
print(ggPacf(dtrain12, lag.max = 24) + labs(title = "PACF of seasonal difference") + plot_theme)

model_s <- Arima(train, order = c(2, 0, 0),
                 seasonal = list(order = c(1, 1, 0), period = 12),
                 method = "ML")
cat("\nSARIMA(2,0,0)(1,1,0)[12]:\n")
print(coeftest(model_s))
checkresiduals(model_s)

model_s_aic <- auto.arima(train, max.d = 1, max.D = 1,
                          max.p = 3, max.q = 3,
                          max.P = 2, max.Q = 2,
                          ic = "aic", stepwise = FALSE, approximation = FALSE)
model_s_bic <- auto.arima(train, max.d = 1, max.D = 1,
                          max.p = 3, max.q = 3,
                          max.P = 2, max.Q = 2,
                          ic = "bic", stepwise = FALSE, approximation = FALSE)

cat("\nauto.arima SARIMA (AIC):", arimaorder(model_s_aic), "\n")
cat("auto.arima SARIMA (BIC):", arimaorder(model_s_bic), "\n")
checkresiduals(model_s_aic)
checkresiduals(model_s_bic)

print(
  autoplot(forecast(model_s, h = h)) +
    autolayer(pi_ts, series = "Actual") +
    labs(title = "SARIMA forecast", x = NULL, y = "%", colour = NULL) +
    plot_theme
)
print(
  autoplot(forecast(model_s_aic, h = h)) +
    autolayer(pi_ts, series = "Actual") +
    labs(title = "auto.arima SARIMA (AIC) forecast",
         x = NULL, y = "%", colour = NULL) +
    plot_theme
)
print(
  autoplot(forecast(model_s_bic, h = h)) +
    autolayer(pi_ts, series = "Actual") +
    labs(title = "auto.arima SARIMA (BIC) forecast",
         x = NULL, y = "%", colour = NULL) +
    plot_theme
)
