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

# Data: Polish GDP (NSA) quarterly growth rate ------------------------------
gdp <- fread(file.path(project_root, "data", "GDP_POLAND_NSA.csv"))
setnames(gdp, c("date", "Y"))
gdp[, date := as.Date(date, format = "%m/%d/%Y")]

logY <- ts(log(gdp$Y), start = c(1995, 1), frequency = 4)
dlogY <- diff(logY) * 100

dlogY <- window(dlogY, start = c(2004, 4), end = c(2019, 4))
train <- window(dlogY, end = c(2018, 4))
test  <- window(dlogY, start = c(2019, 1))

print(
  autoplot(dlogY) +
    labs(title = "Polish GDP growth (NSA, quarterly, %)", x = NULL, y = "%") +
    plot_theme
)

# Box-Jenkins identification ------------------------------------------------
print(ggAcf(train,  lag.max = 12) + labs(title = "ACF")  + plot_theme)
print(ggPacf(train, lag.max = 12) + labs(title = "PACF") + plot_theme)

# Seasonal difference
train_s4 <- diff(train, lag = 4)
print(ggAcf(train_s4,  lag.max = 12) + labs(title = "ACF of seasonal difference")  + plot_theme)
print(ggPacf(train_s4, lag.max = 12) + labs(title = "PACF of seasonal difference") + plot_theme)

# Seasonal ARMA model -------------------------------------------------------
model_s <- Arima(train, order = c(1, 0, 1),
                 seasonal = list(order = c(1, 0, 1), period = 4),
                 method = "ML")
cat("\nSARIMA(1,0,1)(1,0,1)[4]:\n")
print(coeftest(model_s))
checkresiduals(model_s)

resids_s <- residuals(model_s)
print(ggAcf(resids_s,  lag.max = 24) + labs(title = "Residual ACF")  + plot_theme)
print(ggPacf(resids_s, lag.max = 24) + labs(title = "Residual PACF") + plot_theme)

# auto.arima: no seasonal differencing --------------------------------------
model_aic <- auto.arima(train, max.d = 0, max.D = 0,
                        max.p = 4, max.q = 4,
                        ic = "aic", stepwise = FALSE, approximation = FALSE)
model_bic <- auto.arima(train, max.d = 0, max.D = 0,
                        max.p = 4, max.q = 4,
                        ic = "bic", stepwise = FALSE, approximation = FALSE)

cat("\nauto.arima (AIC, D=0):", arimaorder(model_aic), "\n")
cat("auto.arima (BIC, D=0):", arimaorder(model_bic), "\n")
checkresiduals(model_aic)
checkresiduals(model_bic)

# auto.arima: allowing seasonal differencing --------------------------------
model_s_aic <- auto.arima(train, max.d = 0, max.D = 1,
                          max.p = 4, max.q = 4,
                          ic = "aic", stepwise = FALSE, approximation = FALSE)
model_s_bic <- auto.arima(train, max.d = 0, max.D = 1,
                          max.p = 4, max.q = 4,
                          ic = "bic", stepwise = FALSE, approximation = FALSE)

cat("\nauto.arima (AIC, D<=1):", arimaorder(model_s_aic), "\n")
cat("auto.arima (BIC, D<=1):", arimaorder(model_s_bic), "\n")
checkresiduals(model_s_aic)
checkresiduals(model_s_bic)

# Non-seasonal benchmark ----------------------------------------------------
model_ns <- Arima(train, order = c(1, 0, 1), method = "ML")
cat("\nARMA(1,1) without seasonal terms:\n")
print(coeftest(model_ns))

# Forecasts -----------------------------------------------------------------
h <- length(test)

print(
  autoplot(forecast(model_s, h = h)) +
    autolayer(dlogY, series = "Actual") +
    labs(title = "SARIMA(1,0,1)(1,0,1)[4] forecast",
         x = NULL, y = "%", colour = NULL) +
    plot_theme
)
print(
  autoplot(forecast(model_s_aic, h = h)) +
    autolayer(dlogY, series = "Actual") +
    labs(title = "auto.arima seasonal (AIC) forecast",
         x = NULL, y = "%", colour = NULL) +
    plot_theme
)
print(
  autoplot(forecast(model_s_bic, h = h)) +
    autolayer(dlogY, series = "Actual") +
    labs(title = "auto.arima seasonal (BIC) forecast",
         x = NULL, y = "%", colour = NULL) +
    plot_theme
)
print(
  autoplot(forecast(model_ns, h = h)) +
    autolayer(dlogY, series = "Actual") +
    labs(title = "ARMA(1,1) forecast (no seasonal)",
         x = NULL, y = "%", colour = NULL) +
    plot_theme
)

cat("\nIn-sample accuracy of SARIMA(1,0,1)(1,0,1)[4]:\n")
print(accuracy(model_s))
