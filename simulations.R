# =========================
# Sim 1: time-varying variance (robust, no purrr)
# =========================

# ---- Packages ----
need <- c("gamlss", "gamlss.dist", "ggplot2")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
lapply(need, library, character.only = TRUE)

set.seed(123)

# ---- 1) SIMULATE ----
n          <- 800
t_index    <- seq_len(n)
t_scaled   <- (t_index - 1)/(n - 1)   # 0..1
mu0        <- 0
sigma0     <- 2
amp        <- 0.7                     # strength of variance modulation

sigma_t <- sigma0 * (1 + amp * sin(2*pi*t_scaled) + 0.2 * t_scaled)
sigma_t <- pmax(sigma_t, 0.2)         # keep positive
y <- rnorm(n, mean = mu0, sd = sigma_t)

dat <- data.frame(
  t = t_index,
  t_scaled = t_scaled,
  y = y,
  sd_true = sigma_t,
  mu_true = mu0
)

# ---- 2) FIT FULL-SERIES MODELS ----
# M0: constant variance
m0 <- gamlss(
  y ~ 1,
  sigma.fo = ~ 1,   # constant variance 
  family = NO(),
  data = dat,
  trace = FALSE
)

# M1: sigma varies smoothly with time
m1 <- gamlss(
  y ~ 1,                     # Sim 1: mean is constant
  sigma.fo = ~ pb(t_scaled),   # sigma(t)
  family = NO(),
  data = dat,
  trace = FALSE
)

aic_tbl <- data.frame(
  model = c("M0: const sigma", "M1: sigma(t)"),
  AIC   = c(GAIC(m0, k = 2), GAIC(m1, k = 2))
)
print(aic_tbl)

# ---- 3) ROLLING ONE-STEP-AHEAD EVALUATION ----
# Helper: CRPS for Normal (no extra pkg needed)
crps_norm_manual <- function(y, mean, sd) {
  # Gneiting & Raftery (2007) closed form
  if (sd <= 0) return(NA_real_)
  z <- (y - mean)/sd
  sd * (z*(2*pnorm(z)-1) + 2*dnorm(z) - 1/sqrt(pi))
}

k_start <- 200      # start after we have enough training data
k_end   <- n - 1

# Preallocate
K <- k_end - k_start + 1
logscore_M0 <- numeric(K); logscore_M1 <- numeric(K)
crps_M0     <- numeric(K); crps_M1     <- numeric(K)
cov95_M0    <- logical(K); cov95_M1    <- logical(K)

i <- 0
for (k in k_start:k_end) {
  i <- i + 1
  
  # Train/test split
  train <- dat[1:k, , drop = FALSE]
  test1 <- dat[k + 1L, , drop = FALSE]
  
  # Fit M0 safely
  m0_k <- try(
    gamlss(y ~ 1, sigma.fo = ~ 1, family = NO(), data = train, trace = FALSE),
    silent = TRUE
  )
  # Fit M1 safely
  m1_k <- try(
    gamlss(y ~ 1, sigma.fo = ~ pb(t_scaled), family = NO(), data = train, trace = FALSE),
    silent = TRUE
  )
  
  # If either fails (should be rare), carry NA and continue
  if (inherits(m0_k, "try-error") || inherits(m1_k, "try-error")) {
    logscore_M0[i] <- NA_real_; logscore_M1[i] <- NA_real_
    crps_M0[i]     <- NA_real_; crps_M1[i]     <- NA_real_
    cov95_M0[i]    <- NA;       cov95_M1[i]    <- NA
    next
  }
  
  # Predict params for the next point
  p0 <- predictAll(m0_k, newdata = test1, type = "response")
  p1 <- predictAll(m1_k, newdata = test1, type = "response")
  
  mu0_hat <- as.numeric(p0$mu);  sd0_hat <- as.numeric(p0$sigma)
  mu1_hat <- as.numeric(p1$mu);  sd1_hat <- as.numeric(p1$sigma)
  y_obs   <- test1$y
  
  # Log score (log predictive density under Normal)
  logscore_M0[i] <- dNO(y_obs, mu = mu0_hat, sigma = sd0_hat, log = TRUE)
  logscore_M1[i] <- dNO(y_obs, mu = mu1_hat, sigma = sd1_hat, log = TRUE)
  
  # CRPS
  crps_M0[i] <- crps_norm_manual(y_obs, mu0_hat, sd0_hat)
  crps_M1[i] <- crps_norm_manual(y_obs, mu1_hat, sd1_hat)
  
  # 95% coverage
  q0_lo <- qNO(0.025, mu = mu0_hat, sigma = sd0_hat)
  q0_hi <- qNO(0.975, mu = mu0_hat, sigma = sd0_hat)
  q1_lo <- qNO(0.025, mu = mu1_hat, sigma = sd1_hat)
  q1_hi <- qNO(0.975, mu = mu1_hat, sigma = sd1_hat)
  
  cov95_M0[i] <- (y_obs >= q0_lo) && (y_obs <= q0_hi)
  cov95_M1[i] <- (y_obs >= q1_lo) && (y_obs <= q1_hi)
}

summary_tbl <- data.frame(
  model          = c("M0: const sigma", "M1: sigma(t)"),
  mean_log_score = c(mean(logscore_M0, na.rm = TRUE), mean(logscore_M1, na.rm = TRUE)),
  mean_crps      = c(mean(crps_M0,     na.rm = TRUE), mean(crps_M1,     na.rm = TRUE)),
  cover95_rate   = c(mean(cov95_M0,    na.rm = TRUE), mean(cov95_M1,    na.rm = TRUE))
)
print(summary_tbl)

# ---- 4) VISUAL CHECKS ----
# True vs fitted sigma (full-series M1)
sig_hat <- as.numeric(fitted(m1, "sigma"))

print(ggplot(dat, aes(t)) +
        geom_line(aes(y = sd_true)) +
        geom_line(aes(y = sig_hat), colour = "red") +
        labs(title = "Sim 1: SD true (black) vs fitted sigma(t) (red)",
             x = "time", y = "SD"))

# PIT histogram for M1 (quick calibration feel)
pa_full <- predictAll(m1, type = "response")
u <- pNO(dat$y, mu = as.numeric(pa_full$mu), sigma = as.numeric(pa_full$sigma))
print(ggplot(data.frame(u = u), aes(u)) +
        geom_histogram(bins = 20) +
        labs(title = "PIT histogram (M1)", x = "u", y = "count"))
