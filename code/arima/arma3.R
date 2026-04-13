rm(list = ls())
set.seed(13)

suppressPackageStartupMessages({
  library(ggplot2)
  library(forecast)
  library(lmtest)
  library(urca)
})

plot_theme <- theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot", panel.grid.minor = element_blank())

# Simulate AR(3) ------------------------------------------------------------
n_obs <- 240
ar_vec <- c(0.6, 0.2, 0.1)
y <- ts(arima.sim(model = list(ar = ar_vec), n = n_obs, sd = 1),
        frequency = 12, start = c(2000, 1))

print(
  autoplot(y) +
    labs(title = "Simulated AR(3) process", x = NULL, y = NULL) +
    plot_theme
)
print(ggAcf(y,  lag.max = 24) + labs(title = "ACF")  + plot_theme)
print(ggPacf(y, lag.max = 24) + labs(title = "PACF") + plot_theme)

# Train/test split ----------------------------------------------------------
train <- window(y, end = c(2018, 12))
test  <- window(y, start = c(2019, 1))

# Box-Jenkins identification ------------------------------------------------
print(ggAcf(train,  lag.max = 24) + labs(title = "ACF (train)")  + plot_theme)
print(ggPacf(train, lag.max = 24) + labs(title = "PACF (train)") + plot_theme)

cat("\nADF test (type = none, lags = 4):\n")
print(summary(ur.df(train, type = "none", lags = 4)))

# Initial model and alternatives --------------------------------------------
model         <- Arima(train, order = c(2, 0, 0), method = "ML")
model_plus_ar <- Arima(train, order = c(3, 0, 0), method = "ML")
model_plus_ma <- Arima(train, order = c(2, 0, 1), method = "ML")

cat("\nAR(2):\n");      print(coeftest(model))
cat("\nAR(3):\n");      print(coeftest(model_plus_ar))
cat("\nARMA(2,1):\n");  print(coeftest(model_plus_ma))

cat("\nLjung-Box on AR(2) residuals:\n")
print(Box.test(residuals(model), lag = 24, type = "Ljung-Box", fitdf = 2))
checkresiduals(model)

# Information criteria ------------------------------------------------------
model_aic <- auto.arima(train, max.d = 0, max.p = 24, max.q = 24, ic = "aic",
                        seasonal = FALSE, stepwise = FALSE, approximation = FALSE)
model_bic <- auto.arima(train, max.d = 0, max.p = 24, max.q = 24, ic = "bic",
                        seasonal = FALSE, stepwise = FALSE, approximation = FALSE)

cat("\nauto.arima (AIC) selected order:", arimaorder(model_aic), "\n")
cat("auto.arima (BIC) selected order:", arimaorder(model_bic), "\n")
checkresiduals(model_aic)
checkresiduals(model_bic)

# General-to-specific -------------------------------------------------------
model_general  <- Arima(train, order = c(3, 0, 1), method = "ML")
model_specific <- Arima(train, order = c(2, 0, 0), method = "ML")

lr_stat <- -2 * as.numeric(logLik(model_specific) - logLik(model_general))
df_diff <- length(coef(model_general)) - length(coef(model_specific))
p_value <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)

cat("\nLR test (general vs specific):\n")
cat("  statistic =", round(lr_stat, 3),
    "  df =", df_diff,
    "  p-value =", signif(p_value, 4), "\n")

# Impulse response ----------------------------------------------------------
ar_coef <- as.numeric(coef(model)[paste0("ar", 1:2)])
irf <- ts(c(1, ARMAtoMA(ar = ar_coef, ma = numeric(0), lag.max = 24)), start = 0)

print(
  autoplot(irf) +
    labs(title = "Impulse response — AR(2)",
         x = "h", y = expression(psi[h])) +
    plot_theme
)

# Forecasts -----------------------------------------------------------------
h <- length(test)

print(
  autoplot(forecast(model, h = h)) +
    autolayer(y, series = "Actual") +
    labs(title = "AR(2) forecast", x = NULL, y = NULL, colour = NULL) +
    plot_theme
)
print(
  autoplot(forecast(model_aic, h = h)) +
    autolayer(y, series = "Actual") +
    labs(title = "auto.arima (AIC) forecast", x = NULL, y = NULL, colour = NULL) +
    plot_theme
)
print(
  autoplot(forecast(model_bic, h = h)) +
    autolayer(y, series = "Actual") +
    labs(title = "auto.arima (BIC) forecast", x = NULL, y = NULL, colour = NULL) +
    plot_theme
)
