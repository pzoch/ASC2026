rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(forecast)
  library(ggplot2)
  library(zoo)
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

# Moving average forecast: y_t^f = mean(y_{t-k}, ..., y_{t-1})
ma_forecast <- function(y, k) {
  n <- length(y)
  yf <- rep(NA_real_, n)
  for (t in (k + 1):n) {
    yf[t] <- mean(y[(t - k):(t - 1)])
  }
  yf
}

ma_mse <- function(y, k) {
  yf <- ma_forecast(y, k)
  valid <- !is.na(yf)
  mean((y[valid] - yf[valid])^2)
}

project_root <- find_project_root()

# Data: Polish GDP (SA) -----------------------------------------------------
gdp <- fread(file.path(project_root, "data", "GDP_POLAND.csv"))
setnames(gdp, c("date", "Y"))
gdp[, date := as.Date(date, format = "%m/%d/%Y")]

logY <- ts(log(gdp$Y), start = c(1995, 1), frequency = 4)
dlogY <- diff(logY) * 100

print(
  autoplot(logY) +
    labs(title = "Log GDP Poland (SA)", x = NULL, y = "log(GDP)") +
    plot_theme
)
print(
  autoplot(dlogY) +
    labs(title = "GDP growth rate (q/q, %)", x = NULL, y = "%") +
    plot_theme
)

# Train/test split -----------------------------------------------------------
train <- window(dlogY, end = c(2014, 4))
test <- window(dlogY, start = c(2015, 1))
h <- length(test)

# Naive forecasts ------------------------------------------------------------
f_naive <- naive(train, h = h)
f_snaive <- snaive(train, h = h)
f_drift <- rwf(train, h = h, drift = TRUE)
f_mean <- meanf(train, h = h)

print(
  autoplot(train) +
    autolayer(test, series = "Test") +
    autolayer(f_naive$mean, series = "Naive") +
    autolayer(f_snaive$mean, series = "Seasonal naive") +
    autolayer(f_drift$mean, series = "Drift") +
    autolayer(f_mean$mean, series = "Mean") +
    labs(title = "Naive forecasts vs test data", x = NULL, y = "%") +
    plot_theme
)

cat("\nNaive forecast accuracy\n")
print(accuracy(f_naive, test))
print(accuracy(f_snaive, test))
print(accuracy(f_drift, test))
print(accuracy(f_mean, test))

# Moving average: smoothing --------------------------------------------------
# Centered MA with various window sizes (k = 2n+1)
y_vec <- as.numeric(dlogY)
time_idx <- as.Date(as.yearqtr(time(dlogY)))

smooth_df <- data.frame(date = time_idx, original = y_vec)
for (k in c(3, 5, 7, 9)) {
  smooth_df[[paste0("MA_", k)]] <- rollmean(y_vec, k = k, fill = NA, align = "center")
}

smooth_long <- reshape(
  smooth_df,
  direction = "long",
  varying = names(smooth_df)[-1],
  v.names = "value",
  timevar = "series",
  times = names(smooth_df)[-1]
)

print(
  ggplot(smooth_long, aes(x = date, y = value, color = series)) +
    geom_line(linewidth = 0.7) +
    labs(title = "MA smoothing (centered)", x = NULL, y = "%", color = NULL) +
    plot_theme
)

# Moving average: forecasting ------------------------------------------------
# y_t^f = (1/k) * sum(y_{t-k}, ..., y_{t-1})
train_vec <- as.numeric(train)

forecast_df <- data.frame(
  date = as.Date(as.yearqtr(time(train))),
  actual = train_vec
)
for (k in c(2, 4, 6, 8)) {
  forecast_df[[paste0("MA_", k)]] <- ma_forecast(train_vec, k)
}

forecast_long <- reshape(
  forecast_df,
  direction = "long",
  varying = names(forecast_df)[-1],
  v.names = "value",
  timevar = "series",
  times = names(forecast_df)[-1]
)

print(
  ggplot(forecast_long, aes(x = date, y = value, color = series)) +
    geom_line(linewidth = 0.7) +
    labs(title = "MA forecasts (trailing, lagged)", x = NULL, y = "%", color = NULL) +
    plot_theme
)

# Optimal k via MSE ----------------------------------------------------------
k_grid <- 1:12
mse_vals <- sapply(k_grid, function(k) ma_mse(train_vec, k))
cat("\nMA forecast MSE by k:\n")
print(data.frame(k = k_grid, MSE = round(mse_vals, 4)))
cat("Optimal k:", k_grid[which.min(mse_vals)], "\n")

print(
  ggplot(data.frame(k = k_grid, MSE = mse_vals), aes(x = k, y = MSE)) +
    geom_line(linewidth = 0.9, color = "#264653") +
    geom_point(size = 2, color = "#264653") +
    scale_x_continuous(breaks = k_grid) +
    labs(title = "MSE by MA window size k", x = "k", y = "MSE") +
    plot_theme
)

# Simple exponential smoothing -----------------------------------------------
alphas <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.99)
ses_fits <- lapply(alphas, function(a) ses(train, h = h, alpha = a))
ses_opt <- ses(train, h = h)

cat("\nSES: optimal alpha =", ses_opt$model$par["alpha"], "\n")

print(
  autoplot(train) +
    autolayer(test, series = "Test") +
    autolayer(ses_fits[[2]]$mean, series = "alpha=0.1") +
    autolayer(ses_fits[[4]]$mean, series = "alpha=0.5") +
    autolayer(ses_fits[[6]]$mean, series = "alpha=0.99") +
    autolayer(ses_opt$mean, series = "alpha=optimal") +
    labs(title = "SES forecasts", x = NULL, y = "%") +
    plot_theme
)

cat("\nSES accuracy (alpha = optimal)\n")
print(accuracy(ses_opt, test))
