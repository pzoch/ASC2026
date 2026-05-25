rm(list = ls())

packages <- c("vars", "urca", "forecast", "ggplot2", "lmtest", "sandwich")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install missing packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(vars)
library(urca)
library(forecast)
library(ggplot2)
library(lmtest)
library(sandwich)

options(digits = 4)

###############################################################################
# 1. Spurious regression and cointegration: simulation
###############################################################################

set.seed(111)
nobs <- 200

random_walk <- function(n, sd = 1) {
  cumsum(rnorm(n, sd = sd))
}

# Two unrelated I(1) variables.
x_spurious <- ts(random_walk(nobs))
y_spurious <- ts(random_walk(nobs))
spurious_data <- cbind(x_spurious, y_spurious)
colnames(spurious_data) <- c("x_spurious", "y_spurious")

autoplot(spurious_data) +
  ggtitle("Two unrelated random walks")

summary(lm(y_spurious ~ x_spurious))
summary(lm(diff(y_spurious) ~ diff(x_spurious)))

# Two I(1) variables with one common stochastic trend.
z <- random_walk(nobs)
x <- ts(z + rnorm(nobs, sd = sqrt(0.25)))
y <- ts(z + rnorm(nobs, sd = sqrt(0.75)))
coint_data_sim <- cbind(x, y)

autoplot(coint_data_sim) +
  ggtitle("Cointegrated variables with a common stochastic trend")

long_run_sim <- lm(y ~ x)
summary(long_run_sim)

u_hat_sim <- ts(resid(long_run_sim))
autoplot(u_hat_sim) +
  ggtitle("Residual from long-run relation")

# Unit root checks for variables and residuals.
summary(ur.df(x, type = "drift", lags = 4, selectlags = "AIC"))
summary(ur.df(y, type = "drift", lags = 4, selectlags = "AIC"))
summary(ur.df(diff(x), type = "none", lags = 4, selectlags = "AIC"))
summary(ur.df(diff(y), type = "none", lags = 4, selectlags = "AIC"))

# Engle-Granger logic: residuals from the long-run regression should be I(0).
# Use cointegration critical values in applied work, not ordinary ADF values.
eg_residual_test <- ur.df(u_hat_sim, type = "none", lags = 4, selectlags = "AIC")
summary(eg_residual_test)

# Johansen test on the simulated system.
jo_sim <- ca.jo(coint_data_sim, type = "trace", ecdet = "const", K = 2)
summary(jo_sim)

###############################################################################
# 2. Cointegration and VECM: Finnish money demand data
###############################################################################

data("finland", package = "urca")
finland_ts <- ts(finland, start = c(1958, 1), frequency = 4)

# lrm1: log real money M1
# lny: log real income
money_income <- finland_ts[, c("lrm1", "lny")]
colnames(money_income) <- c("m", "y")

autoplot(money_income) +
  ggtitle("Finland: log real money and log real income")

# Unit root tests.
summary(ur.df(money_income[, "m"], type = "drift", lags = 8, selectlags = "AIC"))
summary(ur.df(money_income[, "y"], type = "drift", lags = 8, selectlags = "AIC"))
summary(ur.df(diff(money_income[, "m"]), type = "none", lags = 8, selectlags = "AIC"))
summary(ur.df(diff(money_income[, "y"]), type = "none", lags = 8, selectlags = "AIC"))

# Engle-Granger two-step check.
long_run_fin <- lm(m ~ y, data = as.data.frame(money_income))
summary(long_run_fin)

u_hat_fin <- ts(
  resid(long_run_fin),
  start = start(money_income),
  frequency = frequency(money_income)
)
autoplot(u_hat_fin) +
  ggtitle("Finland: residual from long-run money demand relation")

summary(ur.df(u_hat_fin, type = "none", lags = 8, selectlags = "AIC"))

# Johansen test. K is the lag order in levels; VECM has K - 1 lagged differences.
VARselect(money_income, lag.max = 8, type = "const")

K_vecm <- 3
jo_fin <- ca.jo(
  money_income,
  type = "trace",
  ecdet = "const",
  K = K_vecm,
  season = 4
)
summary(jo_fin)

# Estimate restricted VECM for rank r = 1.
vecm_fin <- cajorls(jo_fin, r = 1)
vecm_fin

# Convert VECM to a VAR representation, useful for forecasts and IRFs.
var_from_vecm_fin <- vec2var(jo_fin, r = 1)
summary(var_from_vecm_fin)
serial.test(var_from_vecm_fin, lags.bg = 8, type = "BG")

irf_vecm_fin <- irf(
  var_from_vecm_fin,
  impulse = "y",
  response = "m",
  n.ahead = 12,
  boot = TRUE,
  runs = 500,
  seed = 123
)
plot(irf_vecm_fin)

predict(var_from_vecm_fin, n.ahead = 8)

###############################################################################
# 3. Manual ECM from the two-step relation
###############################################################################

delta_m <- diff(money_income[, "m"])
delta_y <- diff(money_income[, "y"])
ec_term_lag <- u_hat_fin[-length(u_hat_fin)]

ecm_data <- na.omit(data.frame(
  delta_m = as.numeric(delta_m),
  delta_y = as.numeric(delta_y),
  delta_m_lag1 = c(NA, head(as.numeric(delta_m), -1)),
  delta_y_lag1 = c(NA, head(as.numeric(delta_y), -1)),
  ec_term_lag = as.numeric(ec_term_lag)
))

ecm_fin <- lm(delta_m ~ delta_y + delta_m_lag1 + delta_y_lag1 + ec_term_lag, data = ecm_data)
summary(ecm_fin)
coeftest(ecm_fin, vcov = vcovHAC(ecm_fin))

adjustment <- coef(ecm_fin)["ec_term_lag"]
half_life <- -log(2) / log(1 + adjustment)
half_life
