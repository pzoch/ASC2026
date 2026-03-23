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

# Data: US CPI (NSA, monthly) -----------------------------------------------
cpi <- fread(file.path(project_root, "data", "CPI_US_NSA.csv"))
setnames(cpi, c("date", "Pi"))

Pi_full <- ts(cpi$Pi, start = c(1960, 1), frequency = 12)
Pi <- window(Pi_full, start = c(2000, 1))

print(
  autoplot(Pi) +
    labs(title = "US CPI inflation (NSA, monthly)", x = NULL, y = "%") +
    plot_theme
)

# Residual analysis: benchmark methods ---------------------------------------
f_naive <- naive(Pi)
f_snaive <- snaive(Pi)
f_meanf <- meanf(Pi)
f_ses <- ses(Pi, alpha = 0.25)

checkresiduals(f_naive)
checkresiduals(f_snaive)
checkresiduals(f_meanf)
checkresiduals(f_ses)

# Train/test split -----------------------------------------------------------
train <- window(Pi, end = c(2017, 12))
test <- window(Pi, start = c(2018, 1), end = c(2018, 12))
h <- 12

# Holt-Winters with various parameters --------------------------------------
hw_specs <- list(
  list(alpha = 0.45, beta = 0.1, gamma = 0.50, label = "a=0.45,b=0.1,g=0.5"),
  list(alpha = 0.25, beta = 0.1, gamma = 0.50, label = "a=0.25,b=0.1,g=0.5"),
  list(alpha = 0.15, beta = 0.1, gamma = 0.50, label = "a=0.15,b=0.1,g=0.5")
)

hw_fits <- lapply(hw_specs, function(s) {
  hw(train, h = h, alpha = s$alpha, beta = s$beta, gamma = s$gamma, seasonal = "additive")
})

for (i in seq_along(hw_fits)) {
  print(
    autoplot(hw_fits[[i]]) +
      autolayer(fitted(hw_fits[[i]]), series = "Fitted") +
      autolayer(test, series = "Test") +
      labs(title = paste("HW additive:", hw_specs[[i]]$label), x = NULL, y = "%") +
      plot_theme
  )
}

checkresiduals(hw_fits[[1]])

# Optimal HW (additive — multiplicative fails with negative CPI values) -----
f_hw_opt <- hw(train, h = h, seasonal = "additive")
checkresiduals(f_hw_opt)

cat("\nOptimal HW (additive) parameters:\n")
cat("  alpha =", f_hw_opt$model$par["alpha"], "\n")
cat("  beta  =", f_hw_opt$model$par["beta"], "\n")
cat("  gamma =", f_hw_opt$model$par["gamma"], "\n")

# Benchmarks for comparison
f_naive_h <- naive(train, h = h)
f_snaive_h <- snaive(train, h = h)

cat("\nAccuracy on test set (2018)\n")
for (i in seq_along(hw_fits)) {
  cat("\n", hw_specs[[i]]$label, "\n")
  print(accuracy(hw_fits[[i]], test))
}
cat("\nHW additive (optimal)\n")
print(accuracy(f_hw_opt, test))
cat("\nNaive\n")
print(accuracy(f_naive_h, test))
cat("\nSeasonal naive\n")
print(accuracy(f_snaive_h, test))

# Time series cross-validation -----------------------------------------------
cat("\n\nTime series cross-validation (h=1)\n")
e_naive_cv <- tsCV(train, forecastfunction = naive, h = 1)
e_naive_base <- residuals(naive(train, h = 1))

rmse_cv <- sqrt(mean(e_naive_cv^2, na.rm = TRUE))
rmse_base <- sqrt(mean(e_naive_base^2, na.rm = TRUE))
cat("Naive RMSE - CV:", round(rmse_cv, 4), " base:", round(rmse_base, 4), "\n")

cat("\nCross-validation comparison (h=12, evaluating at horizon 12)\n")
cv_specs <- list(
  list(fn = function(y, h, ...) hw(y, h = h, alpha = 0.45, beta = 0.1, gamma = 0.50, seasonal = "additive"),
       label = "HW(0.45,0.1,0.5)"),
  list(fn = function(y, h, ...) hw(y, h = h, alpha = 0.35, beta = 0.1, gamma = 0.10, seasonal = "additive"),
       label = "HW(0.35,0.1,0.1)"),
  list(fn = function(y, h, ...) hw(y, h = h, alpha = 0.25, beta = 0.05, gamma = 0.10, seasonal = "additive"),
       label = "HW(0.25,0.05,0.1)"),
  list(fn = snaive, label = "Seasonal naive")
)

for (s in cv_specs) {
  e <- tsCV(train, forecastfunction = s$fn, h = 12)
  rmse_h12 <- sqrt(mean(e[, 12]^2, na.rm = TRUE))
  cat(s$label, "- RMSE at h=12:", round(rmse_h12, 4), "\n")
}

# Bootstrap vs normal prediction intervals -----------------------------------
f_snaive_normal <- snaive(train, bootstrap = FALSE, h = h)
f_snaive_boot <- snaive(train, bootstrap = TRUE, h = h)

print(
  autoplot(f_snaive_normal) +
    labs(title = "Seasonal naive: normal intervals") +
    plot_theme
)
print(
  autoplot(f_snaive_boot) +
    labs(title = "Seasonal naive: bootstrap intervals") +
    plot_theme
)
