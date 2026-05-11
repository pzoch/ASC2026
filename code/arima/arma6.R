rm(list = ls())
set.seed(13)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forecast)
  library(lmtest)
  library(urca)
})

plot_theme <- theme_minimal(base_size = 12) +
  theme(plot.title.position = "plot", panel.grid.minor = element_blank())

# Helpers -------------------------------------------------------------------
adf_brief <- function(adf_obj) {
  stat <- adf_obj@teststat[1, "tau2"]
  cval <- adf_obj@cval["tau2", ]
  data.frame(tau = round(stat, 3),
             cv1 = cval["1pct"], cv5 = cval["5pct"], cv10 = cval["10pct"],
             row.names = NULL)
}

adf_bg_check <- function(adf_obj, orders = 1:24) {
  reg     <- adf_obj@testreg
  max_ord <- length(residuals(reg)) - length(coef(reg)) - 1L
  orders  <- orders[orders <= max_ord]
  if (!length(orders)) {
    return(data.frame(order = integer(), BG_stat = numeric(), BG_pval = numeric()))
  }
  bg <- lapply(orders, function(k) bgtest(reg, order = k))
  data.frame(
    order   = orders,
    BG_stat = round(sapply(bg, function(b) as.numeric(b$statistic)), 3),
    BG_pval = round(sapply(bg, function(b) b$p.value), 4)
  )
}

fit_spec <- function(spec, y) {
  Arima(y, order = spec$order,
        seasonal = list(order = spec$seasonal, period = s),
        method = "ML")
}

label_spec <- function(spec) {
  sprintf("(%d,%d,%d)(%d,%d,%d)[%d]",
          spec$order[1], spec$order[2], spec$order[3],
          spec$seasonal[1], spec$seasonal[2], spec$seasonal[3], s)
}

sarima_row <- function(name, fit) {
  lb <- Box.test(residuals(fit), lag = 24, type = "Ljung",
                 fitdf = length(coef(fit)))
  data.frame(
    spec   = name,
    npar   = length(coef(fit)),
    logLik = round(as.numeric(logLik(fit)), 2),
    AIC    = round(AIC(fit), 2),
    BIC    = round(BIC(fit), 2),
    LB24_p = round(lb$p.value, 4),
    row.names = NULL
  )
}

lr_test <- function(general, specific) {
  lr <- -2 * as.numeric(logLik(specific) - logLik(general))
  df <- length(coef(general)) - length(coef(specific))
  data.frame(stat = round(lr, 3), df = df,
             p = round(pchisq(lr, df = df, lower.tail = FALSE), 4),
             row.names = NULL)
}

simulate_sarima <- function(order, seasonal, period, fixed, n,
                            sigma = 1, start = c(2000, 1)) {
  template <- Arima(
    ts(rnorm(n, sd = sigma), frequency = period),
    order    = order,
    seasonal = list(order = seasonal, period = period),
    fixed    = fixed
  )
  template$sigma2 <- sigma^2
  ts(simulate(template, nsim = n, future = FALSE),
     start = start, frequency = period)
}

# Simulate a SARIMA process -------------------------------------------------
# fixed-coef order: ar1..arp, ma1..maq, sar1..sarP, sma1..smaQ
#                   (+ mean only when d = 0 and D = 0).
s     <- 12
n     <- 300
sigma <- 1

sim_order    <- c(2, 1, 2)
sim_seasonal <- c(2, 1, 2)
fixed_coefs  <- c(ar1 = 0.5, ar2 = -0.2,
                  ma1 = 0.3, ma2 =  0.1,
                  sar1 = 0.4, sar2 = 0.2,
                  sma1 = -0.3, sma2 = -0.1)

y_ts <- simulate_sarima(order    = sim_order,
                        seasonal = sim_seasonal,
                        period   = s,
                        fixed    = fixed_coefs,
                        n        = n,
                        sigma    = sigma)

h     <- 12
train <- window(y_ts, end = time(y_ts)[n - h])
test  <- window(y_ts, start = time(y_ts)[n - h + 1])

print(
  autoplot(y_ts) +
    autolayer(test, series = "Test") +
    labs(title = "Simulated SARIMA", x = NULL, y = NULL, colour = NULL) +
    plot_theme
)

# Box-Jenkins identification ------------------------------------------------
dS_train  <- diff(train, lag = s)
ddS_train <- diff(dS_train)

print(ggAcf(train,     lag.max = 36) + labs(title = "ACF  y_t")            + plot_theme)
print(ggPacf(train,    lag.max = 36) + labs(title = "PACF y_t")            + plot_theme)
print(ggAcf(dS_train,  lag.max = 36) + labs(title = "ACF  Delta_12 y_t")   + plot_theme)
print(ggPacf(dS_train, lag.max = 36) + labs(title = "PACF Delta_12 y_t")   + plot_theme)
print(ggAcf(ddS_train, lag.max = 36) + labs(title = "ACF  Delta Delta_12 y_t")  + plot_theme)
print(ggPacf(ddS_train,lag.max = 36) + labs(title = "PACF Delta Delta_12 y_t")  + plot_theme)

# Unit-root tests -----------------------------------------------------------
adf_levels <- ur.df(train,     type = "drift", lags = 12, selectlags = "AIC")
adf_dS     <- ur.df(dS_train,  type = "drift", lags = 12, selectlags = "AIC")
adf_ddS    <- ur.df(ddS_train, type = "drift", lags = 12, selectlags = "AIC")

adf_tbl <- cbind(series = c("y_t", "Delta_12 y_t", "Delta Delta_12 y_t"),
                 rbind(adf_brief(adf_levels), adf_brief(adf_dS), adf_brief(adf_ddS)))
cat("\nADF (drift, AIC-selected lags):\n");                print(adf_tbl, row.names = FALSE)
cat("\nBG on ADF regression residuals -- y_t:\n");         print(adf_bg_check(adf_levels), row.names = FALSE)
cat("\nBG on ADF regression residuals -- Delta_12:\n");    print(adf_bg_check(adf_dS),     row.names = FALSE)
cat("\nBG on ADF regression residuals -- Delta Delta_12:\n"); print(adf_bg_check(adf_ddS), row.names = FALSE)

# SARIMA specs to compare ---------------------------------------------------
specs <- list(
  general  = list(order = c(2, 1, 2), seasonal = c(1, 1, 1)),
  drop_ar2 = list(order = c(1, 1, 2), seasonal = c(1, 1, 1)),
  drop_ma2 = list(order = c(1, 1, 1), seasonal = c(1, 1, 1)),
  drop_sar = list(order = c(1, 1, 1), seasonal = c(0, 1, 1)),
  drop_sma = list(order = c(1, 1, 1), seasonal = c(1, 1, 0))
)

fits   <- lapply(specs, fit_spec, y = train)
labels <- vapply(specs, label_spec, character(1))

cat("\nFitted SARIMA specs:\n")
print(do.call(rbind, Map(sarima_row, labels, fits)), row.names = FALSE)

# LR tests between adjacent (nested) specs ----------------------------------
lr_pairs <- list(
  c("general",  "drop_ar2"),
  c("drop_ar2", "drop_ma2"),
  c("drop_ma2", "drop_sar"),
  c("drop_ma2", "drop_sma")
)
lr_tbl <- do.call(rbind, lapply(lr_pairs, function(p) {
  cbind(general = label_spec(specs[[p[1]]]),
        specific = label_spec(specs[[p[2]]]),
        lr_test(fits[[p[1]]], fits[[p[2]]]))
}))
cat("\nLR tests (general vs nested specific):\n")
print(lr_tbl, row.names = FALSE)

# Information-criteria search -----------------------------------------------
model_aic <- auto.arima(train, max.d = 1, max.D = 1,
                        max.p = 3, max.q = 3, max.P = 2, max.Q = 2,
                        ic = "aic", stepwise = FALSE, approximation = FALSE)
model_bic <- auto.arima(train, max.d = 1, max.D = 1,
                        max.p = 3, max.q = 3, max.P = 2, max.Q = 2,
                        ic = "bic", stepwise = FALSE, approximation = FALSE)
cat("\nauto.arima (AIC):", arimaorder(model_aic),
    "\nauto.arima (BIC):", arimaorder(model_bic), "\n")

# Forecasts -----------------------------------------------------------------
fc_models <- c(setNames(fits, labels),
               list("auto AIC" = model_aic, "auto BIC" = model_bic))

for (nm in names(fc_models)) {
  print(
    autoplot(forecast(fc_models[[nm]], h = h)) +
      autolayer(y_ts, series = "Actual") +
      labs(title = paste("Forecast:", nm), x = NULL, y = NULL, colour = NULL) +
      plot_theme
  )
}

acc_tbl <- do.call(rbind, lapply(names(fc_models), function(nm) {
  a <- accuracy(forecast(fc_models[[nm]], h = h), test)["Test set", ]
  data.frame(model = nm, RMSE = round(a["RMSE"], 3), MAE = round(a["MAE"], 3),
             row.names = NULL)
}))
cat("\nOut-of-sample accuracy:\n")
print(acc_tbl, row.names = FALSE)
