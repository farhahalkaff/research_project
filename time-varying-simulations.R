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


install.packages("usethis")
usethis::create_github_token()

install.packages("gitcreds")
gitcreds::gitcreds_set()
