rm(list = ls())

packages <- c("vars", "forecast", "ggplot2")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install missing packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(vars)
library(forecast)
library(ggplot2)

options(digits = 4)

###############################################################################
# VAR: example on built-in macro data
###############################################################################

data("Canada", package = "vars")
canada <- Canada

autoplot(canada) +
  ggtitle("Canada data from vars package")

# Lag length selection.
VARselect(canada, lag.max = 8, type = "const")

p_var <- 2
var_canada <- VAR(canada, p = p_var, type = "const")
summary(var_canada)

# Coefficients in compact matrix form.
Bcoef(var_canada)
Acoef(var_canada)

###############################################################################
# Residual diagnostics
###############################################################################

serial.test(var_canada, lags.bg = 8, type = "BG")
normality.test(var_canada)
arch.test(var_canada, lags.multi = 4)

###############################################################################
# Stability
###############################################################################

# Values returned by roots() should be below 1 in modulus.
roots(var_canada, modulus = TRUE)
plot(roots(var_canada), ylab = "modulus")
abline(h = 1, lty = 2, col = "red")

###############################################################################
# Forecast
###############################################################################

h <- 8
var_forecast <- predict(var_canada, n.ahead = h, ci = 0.95)
var_forecast
fanchart(var_forecast)

###############################################################################
# Impulse responses and forecast error variance decomposition
###############################################################################

irf_canada <- irf(
  var_canada,
  impulse = "e",
  response = c("prod", "rw", "U"),
  n.ahead = 12,
  boot = TRUE,
  runs = 500,
  seed = 123
)
plot(irf_canada)

fevd_canada <- fevd(var_canada, n.ahead = 12)
plot(fevd_canada)
