rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(forecast)
  library(ggplot2)
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

# Data: Polish GDP (SA and NSA) ----------------------------------------------
gdp_sa <- fread(file.path(project_root, "data", "GDP_POLAND.csv"))
gdp_nsa <- fread(file.path(project_root, "data", "GDP_POLAND_NSA.csv"))
setnames(gdp_sa, c("date", "Y"))
setnames(gdp_nsa, c("date", "Y"))

logY_SA <- ts(log(gdp_sa$Y), start = c(1995, 1), frequency = 4)
logY_NS <- ts(log(gdp_nsa$Y), start = c(1995, 1), frequency = 4)
dlogY_SA <- diff(logY_SA) * 100
dlogY_NS <- diff(logY_NS) * 100

print(
  autoplot(logY_SA, series = "SA") +
    autolayer(logY_NS, series = "NSA") +
    labs(title = "Log GDP Poland", x = NULL, y = "log(GDP)") +
    plot_theme
)

# SES on growth rates --------------------------------------------------------
train_dSA <- window(dlogY_SA, end = c(2014, 4))
test_dSA <- window(dlogY_SA, start = c(2015, 1))
train_dNS <- window(dlogY_NS, end = c(2014, 4))
test_dNS <- window(dlogY_NS, start = c(2015, 1))
h <- length(test_dSA)

f_ses_dSA <- ses(train_dSA, h = h, alpha = 0.3)
f_ses_dNS <- ses(train_dNS, h = h, alpha = 0.3)

print(
  autoplot(train_dNS) +
    autolayer(fitted(f_ses_dNS), series = "Fitted") +
    autolayer(test_dNS, series = "Test") +
    autolayer(f_ses_dNS$mean, series = "Forecast") +
    labs(title = "SES on GDP growth (NSA, alpha=0.3)", x = NULL, y = "%") +
    plot_theme
)

# SES on log levels ----------------------------------------------------------
train_lSA <- window(logY_SA, end = c(2014, 4))
test_lSA <- window(logY_SA, start = c(2015, 1))
train_lNS <- window(logY_NS, end = c(2014, 4))
test_lNS <- window(logY_NS, start = c(2015, 1))

f_ses_lSA <- ses(train_lSA, h = h, alpha = 0.3)
f_ses_lNS <- ses(train_lNS, h = h, alpha = 0.3)

print(
  autoplot(train_lSA) +
    autolayer(fitted(f_ses_lSA), series = "Fitted") +
    autolayer(test_lSA, series = "Test") +
    autolayer(f_ses_lSA$mean, series = "Forecast") +
    labs(title = "SES on log GDP (SA, alpha=0.3)", x = NULL, y = "log(GDP)") +
    plot_theme
)

print(
  autoplot(train_lNS) +
    autolayer(fitted(f_ses_lNS), series = "Fitted") +
    autolayer(test_lNS, series = "Test") +
    autolayer(f_ses_lNS$mean, series = "Forecast") +
    labs(title = "SES on log GDP (NSA, alpha=0.3)", x = NULL, y = "log(GDP)") +
    plot_theme
)

# Holt on log levels ---------------------------------------------------------
f_holt_SA <- holt(train_lSA, h = h, alpha = 0.3, beta = 0.3)
f_holt_NS <- holt(train_lNS, h = h, alpha = 0.3, beta = 0.3)

print(
  autoplot(train_lSA) +
    autolayer(fitted(f_holt_SA), series = "Fitted") +
    autolayer(test_lSA, series = "Test") +
    autolayer(f_holt_SA$mean, series = "Forecast") +
    labs(title = "Holt on log GDP (SA, alpha=0.3, beta=0.3)", x = NULL, y = "log(GDP)") +
    plot_theme
)

print(
  autoplot(train_lNS) +
    autolayer(fitted(f_holt_NS), series = "Fitted") +
    autolayer(test_lNS, series = "Test") +
    autolayer(f_holt_NS$mean, series = "Forecast") +
    labs(title = "Holt on log GDP (NSA, alpha=0.3, beta=0.3)", x = NULL, y = "log(GDP)") +
    plot_theme
)

# Holt-Winters on log levels (additive) -------------------------------------
f_hw_SA <- hw(train_lSA, h = h, alpha = 0.3, beta = 0.3, gamma = 0.1, seasonal = "additive")
f_hw_NS <- hw(train_lNS, h = h, alpha = 0.3, beta = 0.3, gamma = 0.1, seasonal = "additive")

print(
  autoplot(train_lSA) +
    autolayer(fitted(f_hw_SA), series = "Fitted") +
    autolayer(test_lSA, series = "Test") +
    autolayer(f_hw_SA$mean, series = "Forecast") +
    labs(title = "HW additive on log GDP (SA)", x = NULL, y = "log(GDP)") +
    plot_theme
)

print(
  autoplot(train_lNS) +
    autolayer(fitted(f_hw_NS), series = "Fitted") +
    autolayer(test_lNS, series = "Test") +
    autolayer(f_hw_NS$mean, series = "Forecast") +
    labs(title = "HW additive on log GDP (NSA)", x = NULL, y = "log(GDP)") +
    plot_theme
)

# HW decomposition states
plot(f_hw_NS$model$states, main = "HW states (NSA)")

# Varying HW parameters on NSA ----------------------------------------------
param_sets <- list(
  list(alpha = 0.1, beta = 0.1, gamma = 0.5, label = "a=0.1,b=0.1,g=0.5"),
  list(alpha = 0.1, beta = 0.1, gamma = 0.1, label = "a=0.1,b=0.1,g=0.1"),
  list(alpha = 0.2, beta = 0.2, gamma = 0.5, label = "a=0.2,b=0.2,g=0.5")
)

for (ps in param_sets) {
  f <- hw(train_lNS, h = h, alpha = ps$alpha, beta = ps$beta, gamma = ps$gamma, seasonal = "additive")
  print(
    autoplot(train_lNS) +
      autolayer(fitted(f), series = "Fitted") +
      autolayer(test_lNS, series = "Test") +
      autolayer(f$mean, series = "Forecast") +
      labs(title = paste("HW additive (NSA):", ps$label), x = NULL, y = "log(GDP)") +
      plot_theme
  )
}

# Optimal HW ----------------------------------------------------------------
f_hw_opt <- hw(train_lNS, h = h, seasonal = "additive")
cat("\nOptimal HW parameters:\n")
cat("  alpha =", f_hw_opt$model$par["alpha"], "\n")
cat("  beta  =", f_hw_opt$model$par["beta"], "\n")
cat("  gamma =", f_hw_opt$model$par["gamma"], "\n")

print(
  autoplot(train_lNS) +
    autolayer(fitted(f_hw_opt), series = "Fitted") +
    autolayer(test_lNS, series = "Test") +
    autolayer(f_hw_opt$mean, series = "Forecast") +
    labs(title = "HW additive (NSA, optimal parameters)", x = NULL, y = "log(GDP)") +
    plot_theme
)

cat("\nAccuracy comparison (NSA, log levels)\n")
print(accuracy(f_ses_lNS, test_lNS))
print(accuracy(f_holt_NS, test_lNS))
print(accuracy(f_hw_NS, test_lNS))
print(accuracy(f_hw_opt, test_lNS))
