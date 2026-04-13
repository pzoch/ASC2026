rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forecast)
  library(lmtest)
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

# Data: Polish GDP (SA) quarterly growth rate -------------------------------
gdp <- fread(file.path(project_root, "data", "GDP_POLAND.csv"))
setnames(gdp, c("date", "Y"))
gdp[, date := as.Date(date, format = "%m/%d/%Y")]

logY <- ts(log(gdp$Y), start = c(1995, 1), frequency = 4)
dlogY <- diff(logY) * 100

# Window and train/test split ------------------------------------------------
dlogY <- window(dlogY, start = c(2004, 4), end = c(2019, 4))
train <- window(dlogY, end = c(2018, 4))
test  <- window(dlogY, start = c(2019, 1))

print(
  autoplot(train) +
    labs(title = "Polish GDP growth (SA, quarterly, %)", x = NULL, y = "%") +
    plot_theme
)

# Box-Jenkins identification: ACF and PACF ----------------------------------
print(ggAcf(train,  lag.max = 8) + labs(title = "ACF")  + plot_theme)
print(ggPacf(train, lag.max = 8) + labs(title = "PACF") + plot_theme)

# Initial model and alternatives --------------------------------------------
model         <- Arima(train, order = c(3, 0, 0), method = "ML")
model_plus_ar <- Arima(train, order = c(4, 0, 0), method = "ML")
model_plus_ma <- Arima(train, order = c(3, 0, 1), method = "ML")

cat("\nAR(3):\n");      print(coeftest(model))
cat("\nAR(4):\n");      print(coeftest(model_plus_ar))
cat("\nARMA(3,1):\n");  print(coeftest(model_plus_ma))

cat("\nLjung-Box on AR(3) residuals:\n")
print(Box.test(residuals(model), lag = 8, type = "Ljung-Box", fitdf = 3))
checkresiduals(model)

# Information criteria ------------------------------------------------------
model_aic <- auto.arima(train, max.d = 0, max.p = 8, max.q = 8, ic = "aic",
                        seasonal = FALSE, stepwise = FALSE, approximation = FALSE)
model_bic <- auto.arima(train, max.d = 0, max.p = 8, max.q = 8, ic = "bic",
                        seasonal = FALSE, stepwise = FALSE, approximation = FALSE)

cat("\nauto.arima (AIC) selected order:", arimaorder(model_aic), "\n")
cat("auto.arima (BIC) selected order:", arimaorder(model_bic), "\n")
checkresiduals(model_aic)
checkresiduals(model_bic)

# General-to-specific -------------------------------------------------------
model_general  <- Arima(train, order = c(4, 0, 4), method = "ML")
model_specific <- Arima(train, order = c(2, 0, 3), method = "ML")

lr_stat <- -2 * as.numeric(logLik(model_specific) - logLik(model_general))
df_diff <- length(coef(model_general)) - length(coef(model_specific))
p_value <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)

cat("\nLR test (general vs specific):\n")
cat("  statistic =", round(lr_stat, 3),
    "  df =", df_diff,
    "  p-value =", signif(p_value, 4), "\n")
checkresiduals(model_specific)

# Impulse response (MA(∞) representation) -----------------------------------
ar_coef <- as.numeric(coef(model_specific)[paste0("ar", 1:2)])
ma_coef <- as.numeric(coef(model_specific)[paste0("ma", 1:3)])
irf <- ts(c(1, ARMAtoMA(ar = ar_coef, ma = ma_coef, lag.max = 12)), start = 0)

print(
  autoplot(irf) +
    labs(title = "Impulse response — ARMA(2,3)",
         x = "h", y = expression(psi[h])) +
    plot_theme
)

# Forecasts -----------------------------------------------------------------
h <- length(test)
fc <- forecast(model_specific, h = h)

print(
  autoplot(fc) +
    autolayer(test, series = "Actual") +
    labs(title = "ARMA(2,3) forecast vs actual",
         x = NULL, y = "%", colour = NULL) +
    plot_theme
)
