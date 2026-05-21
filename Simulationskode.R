rm(list = ls())
par(mfrow = c(1, 1))
library(MASS)
library(corpcor)

# --------------------------------------------------
# 1. Parametre
# --------------------------------------------------
set.seed(123)
n_assets     <- 10
assets       <- paste0("Asset", 1:n_assets)
mu_annual    <- seq(0.01, 0.19, length.out = n_assets)
mu_true      <- log(1 + mu_annual) / 12
sigma_annual <- 0.10
sigma_asset  <- sigma_annual / sqrt(12)
rho          <- 0.20

Sigma_true <- matrix(rho * sigma_asset^2, n_assets, n_assets)
diag(Sigma_true) <- sigma_asset^2

roll_window <- 60
lookback    <- 12
tau         <- 0.05
confidence  <- 0.50   # 0.50 svarer til He & Litterman
delta_fixed <- 2.5

# --------------------------------------------------
# 2. Hjælpefunktioner
# --------------------------------------------------

build_mom_view <- function(R) {
  R  <- as.matrix(R)
  Tn <- nrow(R)
  
  make_p <- function(u) {
    scores <- colSums(R[(u - lookback + 1):u, , drop = FALSE])
    p      <- scores - mean(scores)
    p      <- pmax(p, 0)
    if (sum(p) > 0) p <- p / sum(p)
    p
  }
  
  ls <- sapply(lookback:(Tn - 1), function(u) sum(make_p(u) * R[u + 1, ]))
  ls <- ls[is.finite(ls)]
  
  list(P_row = make_p(Tn), q = mean(ls))
}

bl_weights <- function(R_insample) {
  w_mkt <- rep(1 / n_assets, n_assets)
  names(w_mkt) <- colnames(R_insample)
  
  Sigma <- cov(R_insample)
  r_m   <- as.numeric(R_insample %*% w_mkt)
  delta <- mean(r_m) / var(r_m)
  if (!is.finite(delta) || delta <= 0) return(w_mkt)
  
  pi   <- as.numeric(delta * Sigma %*% w_mkt)
  view <- tryCatch(build_mom_view(R_insample), error = function(e) NULL)
  if (is.null(view)) return(w_mkt)
  
  P           <- matrix(view$P_row, nrow = 1)
  colnames(P) <- colnames(R_insample)
  eps         <- 1e-8
  Sigma_r     <- Sigma + diag(eps, n_assets)
  tauSigma    <- tau * Sigma_r
  
  view_var <- as.numeric(P %*% tauSigma %*% t(P))
  Omega    <- diag((1 - confidence) / confidence * view_var + eps, 1)
  
  middle <- tryCatch(solve(P %*% tauSigma %*% t(P) + Omega), error = function(e) NULL)
  if (is.null(middle)) return(w_mkt)
  
  mu_bl <- pi + tauSigma %*% t(P) %*% middle %*% (view$q - P %*% pi)
  w     <- pmax(as.numeric((1 / delta) * solve(Sigma_r, mu_bl)), 0)
  names(w) <- colnames(R_insample)
  if (sum(w) > 0) w / sum(w) else w_mkt
}

# --------------------------------------------------
# 3. Simulation: 60 måneder IS → 1 måneds OOS
# --------------------------------------------------

n_sim <- 100000
cat(sprintf("Kører %d simulationer...\n", n_sim))

ret_bl  <- rep(NA_real_, n_sim)
ret_mkt <- rep(NA_real_, n_sim)
w_mat   <- matrix(NA_real_, n_sim, n_assets, dimnames = list(NULL, assets))
P_mat     <- matrix(NA_real_, n_sim, n_assets)
q_vec     <- rep(NA_real_, n_sim)

for (sim in 1:n_sim) {
  log_returns           <- mvrnorm(n = roll_window + 1, mu = mu_true, Sigma = Sigma_true)
  colnames(log_returns) <- assets
  
  R_insample <- log_returns[1:roll_window, ]
  r_oos      <- log_returns[roll_window + 1, ]
  w_bl       <- bl_weights(R_insample)
  
  ret_bl[sim]  <- sum(w_bl * r_oos)
  ret_mkt[sim] <- sum(rep(1 / n_assets, n_assets) * r_oos)
  w_mat[sim, ] <- w_bl
  
  view_s <- tryCatch(build_mom_view(R_insample), error = function(e) NULL)
  if (!is.null(view_s)) {
    P_mat[sim, ] <- view_s$P_row
    q_vec[sim]   <- view_s$q
  }
  
  if (sim %% 10000 == 0) cat(sprintf("  %d / %d\n", sim, n_sim))
}

# --------------------------------------------------
# 4. Resultater
# --------------------------------------------------

excess <- ret_bl - ret_mkt

res <- data.frame(
  Metrik          = c("Sharpe Ratio", "Ann. Afkast (%)",
                      "Ann. Risiko (%)", "Information Ratio"),
  Black_Litterman = c(mean(ret_bl) / sd(ret_bl) * sqrt(12),
                      mean(ret_bl) * 12 * 100,
                      sd(ret_bl)   * sqrt(12) * 100,
                      mean(excess) / sd(excess) * sqrt(12)),
  Equal_Weight    = c(mean(ret_mkt) / sd(ret_mkt) * sqrt(12),
                      mean(ret_mkt) * 12 * 100,
                      sd(ret_mkt)   * sqrt(12) * 100,
                      NA),
  row.names = NULL
)
res$Difference <- res$Black_Litterman - res$Equal_Weight
print(res, digits = 4, na.print = "-")

# --------------------------------------------------
# 5. Diagram
# --------------------------------------------------

deviations <- (colMeans(w_mat) - 0.10) * 100

bp <- barplot(deviations,
              col    = ifelse(deviations >= 0, "steelblue", "tomato"),
              main   = "Basisscenarie: Gennemsnitlig vægtafvigelse fra benchmark (10%)",
              ylab   = "Afvigelse (procentpoint)",
              xlab   = "Aktiv",
              border = NA, las = 2,
              ylim   = c(-8, 12))

abline(h = seq(-10, 10, by = 1),
       col = "grey90", lty = 1, lwd = 0.8)
abline(h = 0, col = "black", lwd = 1.2)

text(bp,
     deviations + ifelse(deviations >= 0, 0.45, -0.45),
     labels = sprintf("%.2f", deviations),
     cex = 0.8, col = "grey30")

# --------------------------------------------------
# 6. Boxplot: Vægtfordeling per aktiv (5%-95%)
# --------------------------------------------------

stats_mat <- apply(w_mat, 2, function(x) {
  c(quantile(x, 0.05),
    quantile(x, 0.25),
    quantile(x, 0.50),
    quantile(x, 0.75),
    quantile(x, 0.95))
})

bp <- list(
  stats = stats_mat,
  n     = rep(nrow(w_mat), n_assets),
  names = assets
)

bxp(bp,
    col     = "steelblue",
    border  = "black",
    main    = "Vægtfordeling per aktiv (5%–95% interval)",
    ylab    = "Vægt",
    xlab    = "Aktiv",
    las     = 2,
    outline = FALSE)

abline(h   = 0.10,
       lty = 2,
       col = "tomato",
       lwd = 1.5)

legend("topleft",
       legend = "Benchmark (10%)",
       col    = "tomato",
       lty    = 2,
       lwd    = 1.5,
       bty    = "n")

# --------------------------------------------------
# 7. View-plot: P og q
# --------------------------------------------------

valid <- which(!is.na(q_vec))
avg_P <- colMeans(P_mat[valid, ], na.rm = TRUE)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# P: Gns. momentum-view vægte
barplot(avg_P * 100,
        names.arg = assets, las = 2,
        col    = "darkorange", border = NA,
        main   = "Gns. momentum-view vægte (P)",
        ylab   = "Vægt i view (%)")

# q: Fordeling på tværs af simulationer
hist(q_vec[valid] * 12 * 100,
     breaks = 60, col = "steelblue", border = "white",
     main   = "Fordeling af momentum-view q",
     xlab   = "q (% p.a.)", ylab = "Antal simulationer")
abline(v   = mean(q_vec[valid], na.rm = TRUE) * 12 * 100,
       col = "tomato", lwd = 2, lty = 2)
legend("topright",
       legend = sprintf("Gns. q = %.2f%%", mean(q_vec[valid], na.rm = TRUE) * 12 * 100),
       col = "tomato", lty = 2, lwd = 2, bty = "n")

par(mfrow = c(1, 1))

cat(sprintf("Std. afvigelse på q: %.2f%%\n", sd(q_vec, na.rm = TRUE) * 12 * 100))

# --------------------------------------------------
# 8. IR som funktion af spredningen i forventet afkast
# --------------------------------------------------

spreads <- seq(0.01, 0.60, by = 0.01)

ir_vec       <- rep(NA_real_, length(spreads))
z_vec        <- rep(NA_real_, length(spreads))
absz_vec     <- rep(NA_real_, length(spreads))
num_vec      <- rep(NA_real_, length(spreads))
denom_vec    <- rep(NA_real_, length(spreads))
fallback_vec <- rep(NA_real_, length(spreads))
neff_vec     <- rep(NA_real_, length(spreads))
active_vec   <- rep(NA_real_, length(spreads))

n_sim_s <- 10000

cat("\nKører IR-analyse...\n")

for (s_idx in seq_along(spreads)) {
  
  spread <- spreads[s_idx]
  
  mu_s      <- seq(-spread, spread, length.out = n_assets) + 0.05
  mu_true_s <- log(1 + mu_s) / 12
  
  ret_bl_s  <- rep(NA_real_, n_sim_s)
  ret_mkt_s <- rep(NA_real_, n_sim_s)
  
  z_s        <- rep(NA_real_, n_sim_s)
  absz_s     <- rep(NA_real_, n_sim_s)
  num_s      <- rep(NA_real_, n_sim_s)
  denom_s    <- rep(NA_real_, n_sim_s)
  fallback_s <- rep(NA_real_, n_sim_s)
  neff_s     <- rep(NA_real_, n_sim_s)
  active_s   <- rep(NA_real_, n_sim_s)
  
  w_mkt <- rep(1 / n_assets, n_assets)
  
  for (sim in 1:n_sim_s) {
    
    log_returns <- mvrnorm(
      n     = roll_window + 1,
      mu    = mu_true_s,
      Sigma = Sigma_true
    )
    
    colnames(log_returns) <- assets
    
    R_insample <- log_returns[1:roll_window, ]
    r_oos      <- log_returns[roll_window + 1, ]
    
    w_bl <- bl_weights(R_insample)
    
    ret_bl_s[sim]  <- sum(w_bl * r_oos)
    ret_mkt_s[sim] <- mean(r_oos)
    
    # Porteføljemål
    neff_s[sim]   <- 1 / sum(w_bl^2)
    active_s[sim] <- 0.5 * sum(abs(w_bl - w_mkt))
    
    # Fallback måles som om vægtene ender lig equal weight
    fallback_s[sim] <- as.numeric(isTRUE(all.equal(
      as.numeric(w_bl),
      w_mkt,
      tolerance = 1e-8
    )))
    
    # View-diagnostik: z-score, tæller og nævner
    view_s <- tryCatch(build_mom_view(R_insample), error = function(e) NULL)
    
    if (!is.null(view_s) && is.finite(view_s$q)) {
      
      Sigma <- cov(R_insample)
      P     <- matrix(view_s$P_row, nrow = 1)
      
      pi <- as.numeric(delta_fixed * Sigma %*% w_mkt)
      
      num_now   <- as.numeric(view_s$q - P %*% pi)
      denom_now <- sqrt(as.numeric(P %*% Sigma %*% t(P)))
      
      num_s[sim]   <- num_now
      denom_s[sim] <- denom_now
      
      if (is.finite(denom_now) && denom_now > 0) {
        z_s[sim]    <- num_now / denom_now
        absz_s[sim] <- abs(z_s[sim])
      }
    }
  }
  
  excess_s <- ret_bl_s - ret_mkt_s
  
  ir_vec[s_idx]       <- mean(excess_s, na.rm = TRUE) / sd(excess_s, na.rm = TRUE) * sqrt(12)
  z_vec[s_idx]        <- mean(z_s, na.rm = TRUE)
  absz_vec[s_idx]     <- mean(absz_s, na.rm = TRUE)
  num_vec[s_idx]      <- mean(num_s, na.rm = TRUE)
  denom_vec[s_idx]    <- mean(denom_s, na.rm = TRUE)
  fallback_vec[s_idx] <- mean(fallback_s, na.rm = TRUE)
  neff_vec[s_idx]     <- mean(neff_s, na.rm = TRUE)
  active_vec[s_idx]   <- mean(active_s, na.rm = TRUE)
  
  cat(sprintf(
    "  Spredning = %.0f%%: IR = %.4f | z = %.4f | fallback = %.1f%% | N_eff = %.2f | AS = %.2f%%\n",
    spread * 100,
    ir_vec[s_idx],
    z_vec[s_idx],
    fallback_vec[s_idx] * 100,
    neff_vec[s_idx],
    active_vec[s_idx] * 100
  ))
}

# --------------------------------------------------
# 8. Plots for IR-analyse
# --------------------------------------------------

# --- Plot 1: IR ---
plot(
  spreads * 100,
  ir_vec,
  type = "b",
  pch  = 16,
  col  = "steelblue",
  lwd  = 2,
  main = "Information Ratio som funktion af spredningen",
  xlab = "Spredning i forventet afkast (%)",
  ylab = "Information Ratio",
  ylim = range(c(ir_vec, 0), na.rm = TRUE)
)

abline(h = 0, lty = 2, col = "tomato", lwd = 1.5)

# --- Plot 2: Gennemsnitlig |z-score| ---
ylim_abs <- range(c(absz_vec, 0, 2), na.rm = TRUE)

plot(
  spreads * 100,
  absz_vec,
  type = "b",
  pch  = 16,
  col  = "purple",
  lwd  = 2,
  main = "Gennemsnitlig |z-score|",
  xlab = "Spredning i forventet afkast (%)",
  ylab = "Gennemsnitlig |z|",
  ylim = ylim_abs
)

abline(h = 2, lty = 2, col = "tomato", lwd = 1.2)


# --- Plot 3: Fallback-rate ---
plot(
  spreads * 100,
  fallback_vec * 100,
  type = "b",
  pch  = 16,
  col  = "tomato",
  lwd  = 2,
  main = "Fallback-rate",
  xlab = "Spredning i forventet afkast (%)",
  ylab = "Fallback-rate (%)",
  ylim = c(0, 100)
)

# --- Plot 4: N_eff ---
plot(
  spreads * 100,
  neff_vec,
  type = "b",
  pch  = 16,
  col  = "darkgreen",
  lwd  = 2,
  main = "Effektivt antal aktiver (N_eff)",
  xlab = "Spredning (%)",
  ylab = "N_eff",
  ylim = c(1, n_assets)
)

abline(h = n_assets, lty = 2, col = "grey40")
abline(h = 1,        lty = 2, col = "tomato")

# --- Plot 5: Active Share ---
plot(
  spreads * 100,
  active_vec,
  type = "b",
  pch  = 16,
  col  = "darkorange",
  lwd  = 2,
  main = "Active Share",
  xlab = "Spredning (%)",
  ylab = "Active Share",
  ylim = c(0, 1)
)


