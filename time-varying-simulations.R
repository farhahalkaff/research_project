# Distributional regression is a statistical framework that models multiple parameters of
# a response distribution as functions of covariates. It extends classical regression by
# allowing predictors to influence not only the mean but also the variance, skewness and
# kurtosis. Simulation studies are essential for understanding how these models behave under
# controlled conditions. They help clarify the assumptions embedded in the model structure and
# reveal how different components contribute to the observed data patterns.

# Here I provide two functions that allow us to simulate a time series with a time-varying
# standard deviation. A Gaussian Process (GP) is a flexible, nonparametric approach for
# modeling smooth, continuous functions over time or space. It is particularly useful for
# capturing latent temporal dynamics in distributional parameters such as variance, skewness
# or kurtosis, without requiring a fixed functional form. The GP framework allows uncertainty
# to be quantified at every point in time, which is essential when modeling complex, evolving
# systems. Its ability to encode assumptions about smoothness and correlation structure makes
# it well suited for simulating and inferring time-varying properties in distributional
# regression models.

# Load required packages
library(ggplot2)
library(tibble)
library(mvnfast)
library(mgcv)

#### Our functions ####
# roxygen2 code tags are used to document R functions in a structured and readable format that
# integrates seamlessly with package development workflows. They allow developers to specify
# function purpose, arguments, return values, examples, and authorship directly above the
# function definition. This improves code transparency, facilitates automated help file
# generation, and ensures reproducibility. Using roxygen2 also supports consistent
# documentation standards across collaborative projects and teaching materials

#' Simulate a Gaussian Process Realization
#'
#' This function generates a single realization from a zero-mean Gaussian
#' Process (GP) using a squared exponential covariance function. It is useful
#' for simulating smooth latent processes such as time-varying effects in
#' distributional regression models.
#'
#' @param N Integer. Number of time points or locations to simulate.
#' @param alpha Numeric. Marginal standard deviation of the GP, controlling
#'   the overall scale of variation.
#' @param rho Numeric. Length scale of the GP, controlling the smoothness
#'   of the process. Larger values produce smoother trajectories.
#'
#' @return A numeric vector of length N containing a single GP realization.
#'
#' @examples
#' # Simulate a smooth latent process over 100 time points
gp_draw <- sim_gp(N = 100, alpha = 1, rho = 20)
plot(gp_draw, type = 'l')
#'
#' @author Nicholas J Clark
sim_gp <- function(N, alpha, rho) {

  if (!requireNamespace("mvnfast", quietly = TRUE)) {
    stop("Package 'mvnfast' is required but not installed.")
  }

  Sigma <- alpha^2 *
    exp(-0.5 * ((outer(1:N, 1:N, "-") / rho)^2)) +
    diag(1e-9, N)

  mvnfast::rmvn(
    1,
    mu = rep(0, N),
    sigma = Sigma
  )[1, ]
}

#' Simulate Gaussian Time Series with GP-Based Time-Varying Variance and
#' Covariate-Driven Mean
#'
#' This function generates a univariate Gaussian time series with a mean that
#' responds to user-specified covariates via a linear model, and a time-varying
#' variance governed by a Gaussian Process (GP). The GP simulates smooth changes
#' in log-variance over time. Users can control the sampling interval and
#' specify covariates for the mean structure.
#'
#' @param covariates A data frame of covariates with n rows. These will be used
#'   to compute the mean via a linear predictor.
#' @param coefs A named numeric vector of regression coefficients. Names must
#'   match column names in `covariates`.
#' @param alpha Numeric. Marginal standard deviation of the GP controlling
#'   the scale of log-variance fluctuations.
#' @param rho Numeric. Length scale of the GP controlling the smoothness of
#'   log-variance fluctuations.
#' @param dt Numeric. Sampling interval (e.g., 1 for daily, 0.1 for sub-daily).
#'
#' @return A tibble with columns:
#'   - time: Time index
#'   - mean: Mean at each time point (from covariates)
#'   - log_sd: Simulated log standard deviation from the GP
#'   - sd: Standard deviation at each time point
#'   - value: Simulated Gaussian observation
#'
#' @examples
#' # Simulate with two covariates affecting the mean
n <- 200
covs <- tibble::tibble(
   x1 = rnorm(n),
   x2 = sin(seq(0, 2 * pi, length.out = n))
)
coefs <- c(x1 = 0.5, x2 = -1)
sim <- sim_time_varying_var(covariates = covs,
 coefs = coefs, alpha = 0.4, rho = 20, dt = 1)
plot(sim$time, sim$value, type = 'l')
#'
#' @author Nicholas J Clark
sim_time_varying_var <- function(
    covariates,
    coefs,
    alpha = 0.5,
    rho = 30,
    dt = 1
) {

  n <- nrow(covariates)

  # Time vector based on sampling interval
  time <- seq(0, by = dt, length.out = n)

  # Compute linear predictor for the mean
  if (!all(names(coefs) %in% colnames(covariates))) {
    stop("All names in 'coefs' must match columns in 'covariates'.")
  }
  mean <- as.numeric(as.matrix(covariates[, names(coefs)]) %*% coefs)

  # Simulate a Gaussian process draw for the time-varying log(SD);
  # we use the log scale here because variance and SD must be non-negative, but
  # we want to allow the underlying function to fluctuate unconstrained. Using
  # exp() converts this back to a non-negative scale
  log_sd <- sim_gp(n, alpha = alpha, rho = rho)
  sd <- exp(log_sd)

  # Simulate Gaussian observations to give the final time series
  value <- rnorm(n, mean = mean, sd = sd)

  # Return results as a tibble
  tibble::tibble(
    time = time,
    mean = mean,
    log_sd = log_sd,
    sd = sd,
    value = value
  )
}

#### Simulate one time series ####

# Simulate with two covariates affecting the mean
n <- 200
covs <- tibble::tibble(
  x1 = rep(1, n)
)
coefs <- c(x1 = 10)
sim <- sim_time_varying_var(
  covariates = covs,
  coefs = coefs,
  alpha = 1,
  rho = 20,
  dt = 1
)

# Plot the series with ggplot2
ggplot(sim, aes(x = time)) +
  geom_line(
    aes(y = value),
    color = "darkred",
    linewidth = 0.8
  ) +
  labs(
    x = "Time",
    y = "Value"
  ) +
  theme_minimal(base_size = 14)

# Use the mgcv package to model how the variance changes over time
?mgcv::gam
?mgcv::gaulss

mod <- gam(
  # For distributional models, we supply a list of formulae
  formula = list(

    # Formula for the mean (empty for now, apart from the intercept)
    value ~ 1,

    # Formula for the SD
    ~ s(time, k = 20)
  ),
  data = sim,
  family = gaulss()
)
summary(mod)

# Predict from the model for each timepoint to get a probabilistic
# prediction of what the SD is, plus any uncertainties
preds <- predict(mod, se.fit = TRUE)

# Extract the predicted SD
predicted_sd <- preds$fit[,2]

# Extract the uncertainty around the predicted SD
se_of_predicted_sd <- preds$se.fit[,2]

# Create mean and 95% upper / lower intervals of SD predictions.
# We need to convert back from the link scale, which in this case is
# the "logb" link (log(SD - b), where b = 0.01 by default)
sim$pred_sd <- exp(predicted_sd - 0.01)
sim$pred_sd_upper <- exp(predicted_sd + 1.96 * se_of_predicted_sd - 0.01)
sim$pred_sd_lower <- exp(predicted_sd - 1.96 * se_of_predicted_sd - 0.01)

# Plot the predicted and simulated SD with ggplot2
ggplot(sim, aes(x = time)) +
  geom_line(
    aes(y = sd),
    color = "darkred",
    linewidth = 0.8
  ) +

  # Add a ribbon for the prediction uncertainty
  geom_ribbon(
    aes(ymin = pred_sd_lower,
        ymax = pred_sd_upper
    ),
    fill = "grey70",
    alpha = 0.4) +
  labs(
    title = "Red = simulated SD; grey = estimated SD",
    x = "Time",
    y = "Value"
  ) +
  theme_minimal(base_size = 14)


### git issues ###
install.packages("usethis")
usethis::create_github_token()

install.packages("gitcreds")
gitcreds::gitcreds_set()


## Andrew's codes 

library(gamlss)
library(gamlss.dist)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

# Define time points
time <- 1:100

# Define time-varying parameters
mu <- 50 + 10 * sin(2 * pi * time / 50)      # location (mean)
sigma <- 5 + 2 * cos(2 * pi * time / 30)     # scale (standard deviation)
nu <- 5 * sin(2 * pi * time / 40)            # shape (skewness)
nu <- rep(0, length(time)) # no skewness 

# Simulate data using the SN1 (Skew Normal type 1) distribution
y <- numeric(length(time))
for (t in time) {
  y[t] <- rSN1(1, mu[t], sigma[t], nu[t])
}

# Create a data frame
df <- data.frame(Time = time, Value = y, Mu = mu, Sigma = sigma, Nu = nu)

# Plot the time series
ggplot(df, aes(x = Time, y = Value)) +
  geom_line(color = "steelblue") +
  labs(title = "Simulated Time Series with GAMLSS-like Parameters",
       x = "Time", y = "Value") +
  theme_minimal()


# Fit GAMLSS model: mean and variance can depend on Time
mod2 <- gamlss(Value ~ cs(Time), 
              sigma.formula = ~ cs(Time), 
              nu.formula = ~1,
              data = df, family = SN1)

# Predict mean and variance for each timepoint
df$PredictedMu <- fitted(mod2, "mu")
df$PredictedSigma <- fitted(mod2, "sigma")
df$PredictedVariance <- df$PredictedSigma^2

# Probabilistic prediction: 95% prediction intervals
df$Lower <- qSN1(0.025, mu = df$PredictedMu, sigma = df$PredictedSigma, nu = coef(mod2, "nu"))
df$Upper <- qSN1(0.975, mu = df$PredictedMu, sigma = df$PredictedSigma, nu = coef(mod2, "nu"))

# Plot observed, fitted mean, and prediction intervals
ggplot(df, aes(x = Time)) +
  geom_line(aes(y = Value), color = "steelblue", alpha = 0.6) +
  geom_line(aes(y = PredictedMu), color = "darkred", size = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "pink", alpha = 0.3) +
  labs(title = "GAMLSS Fit: Mean and 95% Prediction Interval",
       x = "Time", y = "Value") +
  theme_minimal()



## adding seasonaility 
set.seed(42)
time <- 1:100
# Add seasonality: period of 20 units
season_period <- 20
mu <- 50 + 10 * sin(2 * pi * time / 50) + 5 * sin(2 * pi * time / season_period)
sigma <- 5 + 2 * cos(2 * pi * time / 30)
nu <- rep(0, length(time)) # No skewness

y <- numeric(length(time))
for (t in time) {
  y[t] <- rSN1(1, mu[t], sigma[t], nu[t])
}
df <- data.frame(Time = time, Value = y)

# Fit GAMLSS model with seasonality (add sine/cosine terms)
df$sin_season <- sin(2 * pi * df$Time / season_period)
df$cos_season <- cos(2 * pi * df$Time / season_period)

mod4 <- gamlss(Value ~ cs(Time) + sin_season + cos_season,
              sigma.formula = ~ cs(Time),
              nu.start = 0, nu.fix = TRUE,
              data = df, family = SN1)

# Predict mean and variance for each time point
df$PredictedMu <- fitted(mod4, "mu")
df$PredictedSigma <- fitted(mod4, "sigma")
df$Lower <- qSN1(0.025, mu = df$PredictedMu, sigma = df$PredictedSigma, nu = 0)
df$Upper <- qSN1(0.975, mu = df$PredictedMu, sigma = df$PredictedSigma, nu = 0)

# Plot observed, fitted mean, and prediction intervals
ggplot(df, aes(x = Time)) +
  geom_line(aes(y = Value), color = "steelblue", alpha = 0.7) +
  geom_line(aes(y = PredictedMu), color = "firebrick", size = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "pink", alpha = 0.3) +
  labs(title = "Time Series with Seasonality in Mean (No Skewness)",
       x = "Time", y = "Value") +
  theme_minimal()



HR86 <- read.csv("datasets/HelgolandRoads_1986.csv")





####################################
### PHASE 1: Static distribution ###
####################################

# Phase 1: Static distribution → can we recover known parameters?
# ==============================================================
# Compares:
#  A) Mean-only model (gaussian()) — variance inferred from residuals
#  B) Location–scale model (gaulss()) — estimates mu and sigma jointly
#
# Metrics:
#  - Bias & RMSE for mu and sigma
#  - 95% CI coverage for mu
#  - Tail risk calibration: P(Y > c) true vs fitted

library(mgcv)
library(dplyr)
library(purrr)
library(tibble)

set.seed(1)

# -----------------------
# 1) Simulation settings
# -----------------------
n      <- 200           # sample size per dataset
nsim   <- 500           # number of simulated datasets
mu0    <- 0
sigma0 <- 2
thr_c  <- 3             # tail threshold (absolute scale, same units as Y)

# helper: true tail probability under Normal(mu0, sigma0)
true_tail_prob <- 1 - pnorm(thr_c, mean = mu0, sd = sigma0)

# -----------------------
# 2) One simulation step
# -----------------------
fit_once <- function() {
  y <- rnorm(n, mean = mu0, sd = sigma0)
  df <- data.frame(y = y)
  
  # A) Mean-only (gaussian) — constant variance, estimate mu; sigma via residual SD
  mA <- gam(y ~ 1, family = gaussian(), data = df)
  mu_hat_A    <- as.numeric(coef(mA)[1])
  se_mu_A     <- summary(mA)$se[1]
  sigma_hat_A <- sqrt(sum(residuals(mA, type = "response")^2) / (n - 1))
  
  # 95% CI coverage for mu
  cover_mu_A <- (mu0 >= mu_hat_A - 1.96*se_mu_A) & (mu0 <= mu_hat_A + 1.96*se_mu_A)
  
  # Tail prob using model A (Normal with mu_hat_A, sigma_hat_A)
  tail_A <- 1 - pnorm(thr_c, mean = mu_hat_A, sd = sigma_hat_A)
  
  # B) Location–scale (gaulss): two linear predictors, both intercept-only
  mB <- gam(list(y ~ 1, ~ 1), family = gaulss(), data = df)
  
  # Fitted on response scale:
  # predict(..., type="response") returns n x 2 matrix: col1=mu, col2=sigma
  predB <- as.data.frame(predict(mB, type = "response"))
  mu_hat_B    <- mean(predB$V1)      # constant by design, but average to be safe
  sigma_hat_B <- mean(predB$V2)
  
  # SE for mu intercept (on link scale equals identity here)
  se_mu_B <- summary(mB)$se[1]
  cover_mu_B <- (mu0 >= mu_hat_B - 1.96*se_mu_B) & (mu0 <= mu_hat_B + 1.96*se_mu_B)
  
  # Tail prob using model B
  tail_B <- 1 - pnorm(thr_c, mean = mu_hat_B, sd = sigma_hat_B)
  
  tibble(
    mu_hat_A, sigma_hat_A, cover_mu_A = as.integer(cover_mu_A), tail_A,
    mu_hat_B, sigma_hat_B, cover_mu_B = as.integer(cover_mu_B), tail_B
  )
}

# -----------------------
# 3) Run many simulations
# -----------------------
res <- map_dfr(1:nsim, ~fit_once())

# -----------------------
# 4) Summaries
# -----------------------
summ <- tibble(
  model   = c("gaussian(mean-only)", "gaulss(location-scale)"),
  mu_bias = c(mean(res$mu_hat_A - mu0), mean(res$mu_hat_B - mu0)),
  mu_rmse = c(sqrt(mean((res$mu_hat_A - mu0)^2)),
              sqrt(mean((res$mu_hat_B - mu0)^2))),
  mu_cov95 = c(mean(res$cover_mu_A), mean(res$cover_mu_B)),
  sd_bias  = c(mean(res$sigma_hat_A - sigma0), mean(res$sigma_hat_B - sigma0)),
  sd_rmse  = c(sqrt(mean((res$sigma_hat_A - sigma0)^2)),
               sqrt(mean((res$sigma_hat_B - sigma0)^2))),
  tail_true = rep(true_tail_prob, 2),
  tail_pred = c(mean(res$tail_A), mean(res$tail_B)),
  tail_abs_err = abs(c(mean(res$tail_A) - true_tail_prob,
                       mean(res$tail_B) - true_tail_prob))
)
print(summ, n=Inf)

# Optional: quick glance at the sampling skewness of SD estimators
skew_sd_A <- mean( ((res$sigma_hat_A - mean(res$sigma_hat_A)) / sd(res$sigma_hat_A))^3 )
skew_sd_B <- mean( ((res$sigma_hat_B - mean(res$sigma_hat_B)) / sd(res$sigma_hat_B))^3 )
cat("\nSampling skewness of SD estimates:\n",
    "  gaussian():", round(skew_sd_A, 3), "\n",
    "  gaulss()  :", round(skew_sd_B, 3), "\n")


######################################
### PHASE 2: Time-varying variance ###
######################################

# --- Libraries
library(mgcv)
library(gratia)
library(dplyr)
library(ggplot2)
set.seed(42)

# --- 1) Simulate data with constant mean and time‑varying SD
n     <- 800
time  <- seq(0, 10, length.out = n)    # continuous time
mu0   <- 0
# Choose any smooth σ(t); here: seasonal + linear drift, bounded away from 0
sigma_fun <- function(t) { 0.6 + 0.3*sin(2*pi*t) + 0.05*t } 
sigma_t   <- pmax(0.05, sigma_fun(time))

# Optional AR(1) structure (comment out if not needed)
phi <- 0.5
eps <- rnorm(n, 0, 1)
for (i in 2:n) eps[i] <- phi*eps[i-1] + sqrt(1-phi^2)*eps[i]
y <- mu0 + sigma_t * eps

df <- data.frame(y, time)

# --- 2) Fit location–scale GAM (constant mean, smooth SD)
m <- gam(list(y ~ 1, ~ s(time, k = 20)),
         family = gaulss(), data = df, method = "REML")

# If strong autocorrelation remains, consider bam(..., discrete=TRUE) or
# add AR via gamm(); but start simple to validate recovery.

# --- 3) Inspect σ(t): fitted smooth, derivatives (is it increasing?)
# Predict on response scale: columns are mu_hat and sigma_hat
pred <- as.data.frame(predict(m, type = "response"))
df$mu_hat    <- pred$V1
df$sigma_hat <- pred$V2

# Derivatives of sigma smooth
smooths(m)
der_sig <- derivatives(m, term = "s.1(time)", parm = "sigma", type = "central")
der_sig <- der_sig %>% 
  rename(data) %>% 
  mutate(sig_increasing = (lower > 0), sig_decreasing = (upper < 0))
## ERROR ##

library(gratia)
library(dplyr)
library(rlang)

# identify the sigma smooth label
sig_smooth <- grep(":sigma$|\\.2$|sigma$", smooths(m), value = TRUE)[1]

der_sig <- derivatives(
  m,
  term = "s.1(time)",
  type = "central",
  parameter = "sigma"
)

# find whichever column holds the x/time locations
time_col <- intersect(c("data","x","covariate","eval",".eval"), names(der_sig))[1]
if (is.na(time_col)) stop("Couldn't find time column in derivatives(). Names: ", paste(names(der_sig), collapse = ", "))

der_sig <- der_sig %>%
  rename(time = !!sym(time_col)) %>%
  mutate(sig_increasing = (lower > 0),
         sig_decreasing = (upper < 0))
# --- Plots: σ(t) and where it’s significantly rising/falling
ggplot(df, aes(time, sigma_hat)) +
  geom_line() +
  labs(title = "Estimated σ(t)",
       y = "σ̂(t)", x = "time")

ggplot(der_sig, aes(time, derivative)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  labs(title = "Derivative of σ(t) (with 95% CI)",
       y = "dσ̂/dt", x = "time")

# --- 4) “Extremeness” as exceedance probabilities (tail risk)
# Example: probability Y exceeds a fixed threshold c on the right tail.
c_thr <- 1.5  # choose a meaningful threshold in your units
df$p_exceed_right <- 1 - pnorm(c_thr, mean = df$mu_hat, sd = df$sigma_hat)

ggplot(df, aes(time, p_exceed_right)) +
  geom_line() +
  labs(title = "Tail exceedance probability over time",
       y = paste0("P(Y >", c_thr, ")"), x = "time")

# Two-sided exceedance around mean (|Y-μ|>c):
c_abs <- 2
df$p_exceed_two_sided <- 2 * (1 - pnorm(c_abs, mean = 0, sd = df$sigma_hat))
ggplot(df, aes(time, p_exceed_two_sided)) +
  geom_line() +
  labs(title = "Two-sided exceedance probability",
       y = paste0("P(|Y-μ| >", c_abs, ")"), x = "time")

# --- 5) Simple summaries answering “more or less extreme?”
# Trend test via average derivative sign (rough heuristic)
mean_deriv <- der_sig %>% summarize(prop_pos = mean(sig_increasing),
                                    prop_neg = mean(sig_decreasing))
print(mean_deriv)

# Compare early vs late periods (effect size on tail risk)
early <- df %>% slice(1:floor(n*0.3))
late  <- df %>% slice((ceiling(n*0.7)):n)
summary_tail <- tibble(
  period   = c("early","late"),
  mean_sigma = c(mean(early$sigma_hat), mean(late$sigma_hat)),
  mean_p_exceed = c(mean(early$p_exceed_right), mean(late$p_exceed_right))
)
print(summary_tail)


