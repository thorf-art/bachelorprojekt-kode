rm(list = ls())
par(mfrow = c(1, 1))

library(readxl)

# ------------------------------------------------------------
# 1. Indlæs data
# ------------------------------------------------------------

raw <- read_excel(
  "/Users/thor/OneDrive/Skrivebord/6.semester/bachelor/Data/Data10.xlsx",
  sheet = 3
)
dates    <- as.Date(raw[[1]])
prices   <- as.matrix(raw[, -1])
assets   <- colnames(raw)[-1]
n_assets <- length(assets)

log_ret <- log(prices[-1, ] / prices[-nrow(prices), ])
colnames(log_ret) <- assets
dates_r <- dates[-1]
T_total <- nrow(log_ret)

# Risikofri rente (FRED TB3MS)
rf_raw <- read.csv(
  "/Users/thor/OneDrive/Skrivebord/6.semester/bachelor/Data/TB3MS.csv"
)
rf_raw$DATE       <- as.Date(rf_raw$observation_date)
rf_raw$rf_monthly <- log(1 + as.numeric(rf_raw$TB3MS) / 100) / 12
rf_map     <- setNames(rf_raw$rf_monthly, format(rf_raw$DATE, "%Y-%m"))
rf_monthly <- as.numeric(rf_map[format(dates_r, "%Y-%m")])

if (any(is.na(rf_monthly))) stop("Manglende risikofri rente for nogle måneder.")

log_ret_excess <- sweep(log_ret, 1, rf_monthly, "-")

cat(sprintf("Aktiver: %d | Observationer: %d | Periode: %s - %s\n",
            n_assets, T_total,
            format(dates_r[1], "%Y-%m"),
            format(dates_r[T_total], "%Y-%m")))


# ------------------------------------------------------------
# 2. Parametre
# ------------------------------------------------------------

roll_window  <- 36       # In-sample-vindue (måneder)
lookback     <- 12       # Momentum-lookback
tau          <- 0.05     # Black-Litterman skaleringsparameter
confidence   <- 0.50     # Baseline confidence i momentum-viewet


# Datadrevet shrinkage:
w_grid       <- seq(0, 1, by = 0.01)
n_init       <- 60       # Burn-in (måneder)
w_init       <- 0.50     # Neutral initialværdi
ir_tolerance <- 0.05     # Tolerance om maksimal Information Ratio

# Indeksering af OOS-perioden
n_oos     <- T_total - roll_window
oos_dates <- dates_r[(roll_window + 1):T_total]
rf_oos    <- rf_monthly[(roll_window + 1):T_total]
ret_during_oos <- log_ret[(roll_window + 1):T_total, ]

w_eq <- rep(1 / n_assets, n_assets)
names(w_eq) <- assets


# ------------------------------------------------------------
# 3. Hjælpefunktioner
# ------------------------------------------------------------

# Konstruér momentum-view (long-only P-vektor og forventet q)
build_mom_view <- function(R) {
  R  <- as.matrix(R)
  Tn <- nrow(R)
  
  make_p <- function(u) {
    scores <- colSums(R[(u - lookback + 1):u, , drop = FALSE])
    p      <- pmax(scores - mean(scores), 0)
    if (sum(p) > 0) p / sum(p) else p
  }
  
  ls <- sapply(lookback:(Tn - 1), function(u) sum(make_p(u) * R[u + 1, ]))
  ls <- ls[is.finite(ls)]
  
  list(P_row = make_p(Tn), q = mean(ls))
}

# Black-Litterman-vægte med EPO-shrinkage og momentum-view
bl_weights <- function(R_insample, shrink_w, confidence = 0.50, view = NULL) {
  
  w_mkt <- rep(1 / n_assets, n_assets)
  names(w_mkt) <- colnames(R_insample)
  
  if (confidence <= 0) return(w_mkt)
  
  # EPO-shrinkage af korrelationsmatricen
  Sigma_sample <- cov(R_insample)
  vol  <- sqrt(diag(Sigma_sample))
  Corr <- cov2cor(Sigma_sample)
  Corr_shrunk <- (1 - shrink_w) * Corr + shrink_w * diag(n_assets)
  Sigma <- diag(vol) %*% Corr_shrunk %*% diag(vol)
  
  # Risikoaversion fra equal weight-porteføljen
  r_m   <- as.numeric(R_insample %*% w_mkt)
  delta <- mean(r_m) / var(r_m)
  if (!is.finite(delta) || delta <= 0) return(w_mkt)
  
  # Prior via reverse optimization
  pi_vec <- as.numeric(delta * Sigma %*% w_mkt)
  
  # Momentum-view (Q-filter: bortfald hvis q <= 0)
  if (is.null(view)) view <- tryCatch(build_mom_view(R_insample), error = function(e) NULL)
  if (is.null(view) || !is.finite(view$q) || view$q <= 0) return(w_mkt)
  
  P <- matrix(view$P_row, nrow = 1)
  colnames(P) <- colnames(R_insample)
  
  eps      <- 1e-8
  Sigma_r  <- Sigma + diag(eps, n_assets)
  tauSigma <- tau * Sigma_r
  
  c_clip <- max(min(confidence, 0.999), 0.001)
  Omega  <- ((1 - c_clip) / c_clip) * (P %*% tauSigma %*% t(P)) + diag(eps, 1)
  
  middle <- tryCatch(solve(P %*% tauSigma %*% t(P) + Omega), error = function(e) NULL)
  if (is.null(middle)) return(w_mkt)
  
  mu_bl <- pi_vec + tauSigma %*% t(P) %*% middle %*% as.numeric(view$q - P %*% pi_vec)
  
  # Long-only og normalisering
  w <- tryCatch(
    pmax(as.numeric((1 / delta) * solve(Sigma_r, mu_bl)), 0),
    error = function(e) w_mkt
  )
  names(w) <- colnames(R_insample)
  if (sum(w) > 0) w / sum(w) else w_mkt
}

# Benchmarkstrategier
inv_vol_weights <- function(Sigma_sample) {
  w <- 1 / sqrt(diag(Sigma_sample))
  as.numeric(w / sum(w))
}
min_var_weights <- function(Sigma_sample) {
  w <- tryCatch(solve(Sigma_sample + diag(1e-8, n_assets), rep(1, n_assets)),
                error = function(e) rep(1 / n_assets, n_assets))
  w <- pmax(as.numeric(w / sum(w)), 0)
  if (sum(w) > 0) w / sum(w) else rep(1 / n_assets, n_assets)
}

# Performance- og porteføljemål
ann_metrics <- function(ret, ret_bench = NULL) {
  ex <- ret - rf_oos
  list(
    ann_excess = mean(ex) * 12 * 100,
    ann_vol    = sd(ret) * sqrt(12) * 100,
    sharpe     = mean(ex) / sd(ex) * sqrt(12),
    ir         = if (is.null(ret_bench)) NA
    else mean(ret - ret_bench) / sd(ret - ret_bench) * sqrt(12)
  )
}
calc_drawdown <- function(ret) {
  wealth <- exp(cumsum(ret))
  mean(wealth / cummax(wealth) - 1, na.rm = TRUE)
}
calc_turnover <- function(w_mat, ret_holding) {
  # ret_holding: log-afkast i hver holding-periode (n_oos x n_assets)
  n_per <- nrow(w_mat)
  to <- numeric(n_per - 1)
  
  for (i in 2:n_per) {
    r_simple <- exp(ret_holding[i - 1, ]) - 1
    drift    <- w_mat[i - 1, ] * (1 + r_simple)
    if (sum(drift, na.rm = TRUE) > 0) drift <- drift / sum(drift)
    to[i - 1] <- 0.5 * sum(abs(w_mat[i, ] - drift))
  }
  
  mean(to, na.rm = TRUE) * 12   # annualiseret
}
calc_neff <- function(w_mat) mean(1 / rowSums(w_mat^2), na.rm = TRUE)


# ------------------------------------------------------------
# 4. Pass 1: BL-afkast for hver kandidat w_EPO
# ------------------------------------------------------------

cat(sprintf("\nPass 1: BL-afkast for %d kandidater af w_EPO...\n", length(w_grid)))

ret_bl_grid <- matrix(NA_real_, n_oos, length(w_grid))

for (i in seq_len(n_oos)) {
  R_is <- log_ret_excess[i:(i + roll_window - 1), ]
  colnames(R_is) <- assets
  r_oos <- log_ret[i + roll_window, ]
  
  for (j in seq_along(w_grid)) {
    w_j <- bl_weights(R_is, w_grid[j], confidence)
    ret_bl_grid[i, j] <- sum(w_j * r_oos)
  }
}


# ------------------------------------------------------------
# 5. Pass 2: Datadrevet w_EPO + BL + benchmarks
# ------------------------------------------------------------

cat("Pass 2: Datadrevet w_EPO baseret på IR vs EW + benchmarks...\n")

ret_bl   <- ret_ew   <- ret_invvol <- ret_minvar <- rep(NA_real_, n_oos)
w_bl_mat <- w_ew_mat <- w_invvol_mat <- w_minvar_mat <-
  matrix(NA_real_, n_oos, n_assets, dimnames = list(NULL, assets))

w_epo_chosen <- rep(NA_real_, n_oos)
P_mat <- matrix(NA_real_, n_oos, n_assets, dimnames = list(NULL, assets))
q_vec <- rep(NA_real_, n_oos)

for (i in seq_len(n_oos)) {
  R_is <- log_ret_excess[i:(i + roll_window - 1), ]
  colnames(R_is) <- assets
  r_oos <- log_ret[i + roll_window, ]
  
  # Benchmarks (regnes først, så ret_ew[1:(i-1)] er tilgængelig)
  Sigma_sample <- cov(R_is)
  w_iv  <- inv_vol_weights(Sigma_sample)
  w_mv  <- min_var_weights(Sigma_sample)
  
  ret_ew[i]     <- sum(w_eq * r_oos)
  ret_invvol[i] <- sum(w_iv * r_oos)
  ret_minvar[i] <- sum(w_mv * r_oos)
  
  w_ew_mat[i, ]     <- w_eq
  w_invvol_mat[i, ] <- w_iv
  w_minvar_mat[i, ] <- w_mv
  
  # Vælg w_EPO ud fra historisk IR vs EW
  if (i <= n_init) {
    w_epo_chosen[i] <- w_init
  } else {
    irs <- sapply(seq_along(w_grid), function(j) {
      active <- ret_bl_grid[1:(i - 1), j] - ret_ew[1:(i - 1)]
      active <- active[is.finite(active)]
      if (length(active) < 12 || sd(active) == 0) return(NA_real_)
      mean(active) / sd(active)
    })
    
    if (all(!is.finite(irs))) {
      w_epo_chosen[i] <- w_init
    } else {
      max_ir <- max(irs, na.rm = TRUE)
      top    <- which(irs >= max_ir - ir_tolerance)
      w_epo_chosen[i] <- mean(w_grid[top])
    }
  }
  
  # BL-portefølje + gem view til senere analyse
  view_i <- tryCatch(build_mom_view(R_is), error = function(e) NULL)
  if (!is.null(view_i)) {
    P_mat[i, ] <- view_i$P_row
    q_vec[i]   <- view_i$q
  }
  
  w_bl <- bl_weights(R_is, w_epo_chosen[i], confidence, view = view_i)
  ret_bl[i]     <- sum(w_bl * r_oos)
  w_bl_mat[i, ] <- w_bl
}


# ============================================================
# OUTPUT: følger analysens rækkefølge
# ============================================================

# ------------------------------------------------------------
# 3.1 Korrelationsstruktur og valg af shrinkage-niveau
# ------------------------------------------------------------

Sigma_full <- cov(log_ret)
Corr_full  <- cov2cor(Sigma_full)
off_corr   <- Corr_full[upper.tri(Corr_full)]
vol_full   <- sqrt(diag(Sigma_full))
div_ratio  <- sum(w_eq * vol_full) /
  sqrt(as.numeric(t(w_eq) %*% Sigma_full %*% w_eq))

corr_summary <- data.frame(
  Metrik = c("Gns. absolut korrelation",
             "Andel |corr| > 0.50 (%)",
             "Diversification ratio EW"),
  Værdi  = c(mean(abs(off_corr)),
             mean(abs(off_corr) > 0.5) * 100,
             div_ratio)
)
cat("\n=== Korrelationsstruktur (fuldt sample) ===\n")
print(corr_summary, digits = 4, row.names = FALSE)

# Histogram over parvise korrelationer
h <- hist(off_corr, breaks = 30, col = "steelblue", border = "white",
          main = "Fordeling af parvise korrelationer",
          xlab = "Parvis korrelation", ylab = "Antal aktivpar")
abline(v = mean(off_corr), col = "tomato", lwd = 2, lty = 2)
abline(v = 0, col = "black")
text(0, max(h$counts) * 0.9, "0-korrelation", pos = 4, cex = 0.9)

# Datadrevet w_EPO over tid
plot(oos_dates, w_epo_chosen, type = "l", col = "steelblue",
     ylim = c(0, 1), xlab = "Dato", ylab = expression(w[EPO]),
     main = "Datadrevet shrinkage-parameter over tid")
abline(v = oos_dates[n_init], col = "firebrick", lty = 3)
legend("bottomright",
       c("Valgt w_EPO", "Slut på burn-in"),
       col = c("steelblue", "firebrick"),
       lty = c(1, 3), bty = "n")


# ------------------------------------------------------------
# 3.2 Modelvalidering: simpel momentum + 11+1-test
# ------------------------------------------------------------
# Simpel long-only momentumstrategi (top 1/3)
ret_mom <- ret_ew_simple <- rep(NA_real_, n_oos)
for (i in seq_len(n_oos)) {
  R_is  <- log_ret_excess[i:(i + roll_window - 1), ]
  r_oos <- as.numeric(log_ret[i + roll_window, ])
  scores     <- colSums(tail(R_is, lookback))
  top_assets <- which(scores >= quantile(scores, 2/3))
  w_mom <- numeric(n_assets); w_mom[top_assets] <- 1 / length(top_assets)
  ret_mom[i]       <- sum(w_mom * r_oos)
  ret_ew_simple[i] <- sum(w_eq  * r_oos)
}

cum_mom <- exp(cumsum(ret_mom))
cum_ew  <- exp(cumsum(ret_ew_simple))

plot(oos_dates, cum_mom, type = "l", col = "steelblue", lwd = 2,
     ylim = range(c(cum_mom, cum_ew)),
     xlab = "", ylab = "Kumuleret afkast (kr. 1 investeret)",
     main = "Simpel momentumstrategi vs equal weight")
lines(oos_dates, cum_ew, col = "firebrick", lwd = 2, lty = 2)
legend("topleft",
       c("Momentum (top 1/3)", "Equal Weight"),
       col = c("steelblue", "firebrick"),
       lty = c(1, 2), lwd = 2, bty = "n")

# 11+1-test: BL med fremadskuende P (q forbliver IS-baseret)
ret_cheat   <- rep(NA_real_, n_oos)
w_cheat_mat <- matrix(NA_real_, n_oos, n_assets, dimnames = list(NULL, assets))

for (i in seq_len(n_oos)) {
  R_is <- log_ret_excess[i:(i + roll_window - 1), ]
  colnames(R_is) <- assets
  r_next_ex <- log_ret_excess[i + roll_window, ]
  r_oos     <- log_ret[i + roll_window, ]
  
  # P bruger udvidet data (11 IS + 1 OOS), q bruger kun IS
  R_ext   <- rbind(R_is, matrix(r_next_ex, nrow = 1,
                                dimnames = list(NULL, assets)))
  view_ch <- list(P_row = build_mom_view(R_ext)$P_row,
                  q     = build_mom_view(R_is)$q)
  
  w_ch <- bl_weights(R_is, w_epo_chosen[i], confidence, view = view_ch)
  ret_cheat[i]     <- sum(w_ch * r_oos)
  w_cheat_mat[i, ] <- w_ch
}

# Vægtafvigelse i 11+1-test
dev_cheat <- (colMeans(w_cheat_mat, na.rm = TRUE) - 1 / n_assets) * 100

par(mar = c(12, 4, 4, 2))

bp <- barplot(dev_cheat,
              names.arg = assets,
              las       = 2,
              col       = ifelse(dev_cheat >= 0, "steelblue", "tomato"),
              border    = NA,
              cex.names = 0.7,
              ylim      = c(min(dev_cheat) - 0.5, max(dev_cheat) + 0.5),
              ylab      = "Afvigelse fra EW (procentpoint)",
              main      = "Gennemsnitlig vægtafvigelse fra Equal Weight benchmark - 11+1-test")
# Horisontale gridlinjer ved hver 0,5
abline(h = seq(floor(min(dev_cheat)) - 0.5,
               ceiling(max(dev_cheat)) + 0.5,
               by = 0.5),
       col = "grey90", lty = "dotted")
abline(h = 0, col = "black", lwd = 1.2)
# Værdier over/under barerne
text(x      = bp,
     y      = dev_cheat + ifelse(dev_cheat >= 0, 0.08, -0.08),
     labels = sprintf("%.2f", dev_cheat),
     cex    = 0.6,
     pos    = ifelse(dev_cheat >= 0, 3, 1))
par(mar = c(5, 4, 4, 2))

# 11+1 performance-tabel
m_bl    <- ann_metrics(ret_bl,    ret_ew)
m_ew    <- ann_metrics(ret_ew)
m_cheat <- ann_metrics(ret_cheat, ret_ew)

res_cheat <- data.frame(
  Metrik = c("Ann. Excess Afkast (%)", "Ann. Volatilitet (%)",
             "Sharpe Ratio", "Information Ratio vs EW"),
  BL_11plus1   = c(m_cheat$ann_excess, m_cheat$ann_vol, m_cheat$sharpe, m_cheat$ir),
  Equal_Weight = c(m_ew$ann_excess,    m_ew$ann_vol,    m_ew$sharpe,    NA)
)
cat("\n=== 11+1-test ===\n")
print(res_cheat, digits = 3, na.print = "-", row.names = FALSE)


# ------------------------------------------------------------
# 3.3 Realistisk out-of-sample performance
# ------------------------------------------------------------

m_iv <- ann_metrics(ret_invvol, ret_ew)
m_mv <- ann_metrics(ret_minvar, ret_ew)

res_full <- data.frame(
  Metrik = c("Ann. Excess Afkast (%)", "Ann. Volatilitet (%)",
             "Sharpe Ratio", "Information Ratio vs EW",
             "Gns. Drawdown (%)", "Ann. Turnover (%)", "Gns. N_eff"),
  Black_Litterman = c(m_bl$ann_excess, m_bl$ann_vol, m_bl$sharpe, m_bl$ir,
                      calc_drawdown(ret_bl) * 100,
                      calc_turnover(w_bl_mat, ret_during_oos) * 100,
                      calc_neff(w_bl_mat)),
  Inv_Vol = c(m_iv$ann_excess, m_iv$ann_vol, m_iv$sharpe, m_iv$ir,
              calc_drawdown(ret_invvol) * 100,
              calc_turnover(w_invvol_mat, ret_during_oos) * 100,
              calc_neff(w_invvol_mat)),
  Min_Var = c(m_mv$ann_excess, m_mv$ann_vol, m_mv$sharpe, m_mv$ir,
              calc_drawdown(ret_minvar) * 100,
              calc_turnover(w_minvar_mat, ret_during_oos) * 100,
              calc_neff(w_minvar_mat)),
  Equal_Weight = c(m_ew$ann_excess, m_ew$ann_vol, m_ew$sharpe, NA,
                   calc_drawdown(ret_ew) * 100,
                   calc_turnover(w_ew_mat, ret_during_oos) * 100,
                   calc_neff(w_ew_mat))

)
cat("\n=== Performance-tabel (realistisk OOS) ===\n")
print(res_full, digits = 3, na.print = "-", row.names = FALSE)


# ------------------------------------------------------------
# 3.4 BL mod EW: vægtafvigelse, P-vektor og q-fordeling
# ------------------------------------------------------------

# Vægtafvigelse fra EW
dev_bl <- (colMeans(w_bl_mat, na.rm = TRUE) - 1 / n_assets) * 100

par(mar = c(12, 4, 4, 2))

bp_bl <- barplot(dev_bl,
                 names.arg = assets,
                 las = 2,
                 col = ifelse(dev_bl >= 0, "steelblue", "tomato"),
                 border = NA,
                 cex.names = 0.7,
                 ylim = c(min(dev_bl, na.rm = TRUE) - 0.5,
                          max(dev_bl, na.rm = TRUE) + 0.5),
                 ylab = "Afvigelse fra EW (procentpoint)",
                 main = "Vægtafvigelse fra EW - realistisk BL")

abline(h = seq(floor(min(dev_bl, na.rm = TRUE)) - 0.5,
               ceiling(max(dev_bl, na.rm = TRUE)) + 0.5,
               by = 0.5),
       col = "grey90", lty = "dotted")

abline(h = 0, col = "black", lwd = 1.2)

text(x = bp_bl,
     y = dev_bl + ifelse(dev_bl >= 0, 0.08, -0.08),
     labels = sprintf("%.2f", dev_bl),
     cex = 0.6,
     pos = ifelse(dev_bl >= 0, 3, 1))

par(mar = c(5, 4, 4, 2))

# Gns. P-vægte og fordeling af q
valid <- which(is.finite(q_vec) & rowSums(P_mat, na.rm = TRUE) > 0)
avg_P <- colMeans(P_mat[valid, , drop = FALSE], na.rm = TRUE)

par(mfrow = c(1, 1), mar = c(12, 4, 4, 2))
barplot(avg_P * 100, names.arg = assets, las = 2,
        col = "darkorange", border = NA, cex.names = 0.7,
        ylab = "Vægt i view (%)",
        main = "Gns. P-vægte (momentum-view)")

par(mar = c(5, 4, 4, 2))
hist(q_vec[valid] * 12 * 100, breaks = 40, col = "steelblue", border = "white",
     main = "Fordeling af momentum-view q",
     xlab = "q (% p.a.)", ylab = "Antal observationer")
abline(v = mean(q_vec[valid], na.rm = TRUE) * 12 * 100,
       col = "tomato", lwd = 2, lty = 2)
legend("topright",
       sprintf("Gns. q = %.2f%%",
               mean(q_vec[valid], na.rm = TRUE) * 12 * 100),
       col = "tomato", lty = 2, lwd = 2, bty = "n")
par(mfrow = c(1, 1))


# ------------------------------------------------------------
# 3.5 Sensitivitetsanalyse: confidence
# ------------------------------------------------------------

conf_levels <- c(0.25, 0.50, 0.75)

run_confidence <- function(conf) {
  ret_c <- rep(NA_real_, n_oos)
  W_c   <- matrix(NA_real_, n_oos, n_assets, dimnames = list(NULL, assets))
  for (i in seq_len(n_oos)) {
    R_is  <- log_ret_excess[i:(i + roll_window - 1), ]
    colnames(R_is) <- assets
    r_oos <- log_ret[i + roll_window, ]
    w_c   <- bl_weights(R_is, w_epo_chosen[i], conf)
    ret_c[i]  <- sum(w_c * r_oos)
    W_c[i, ]  <- w_c
  }
  list(ret = ret_c, W = W_c)
}

conf_results <- lapply(conf_levels, run_confidence)

res_conf <- data.frame(
  Metrik = c("Ann. Excess Afkast (%)", "Ann. Volatilitet (%)",
             "Sharpe Ratio", "Information Ratio")
)
for (k in seq_along(conf_levels)) {
  m <- ann_metrics(conf_results[[k]]$ret, ret_ew)
  res_conf[[sprintf("BL_%d%%", conf_levels[k] * 100)]] <-
    c(m$ann_excess, m$ann_vol, m$sharpe, m$ir)
}
cat("\n=== Confidence-sensitivitet ===\n")
print(res_conf, digits = 3, row.names = FALSE)

# Gns. aktivvægte per confidence-niveau
w_means <- sapply(conf_results, function(x) colMeans(x$W, na.rm = TRUE))
colnames(w_means) <- sprintf("%d%%", conf_levels * 100)
rownames(w_means) <- assets
ord <- order(w_means[, ncol(w_means)], decreasing = TRUE)

par(mar = c(10, 4, 4, 1))
barplot(t(w_means[ord, ]) * 100, beside = TRUE,
        col = c("steelblue", "darkorange", "seagreen3"),
        names.arg = rownames(w_means)[ord], las = 2, cex.names = 0.65,
        ylab = "Gennemsnitlig vægt (%)",
        main = "Gennemsnitlig aktivvægt per confidence-niveau",
        ylim = c(0, max(w_means) * 100 + 1))
abline(h = 100 / n_assets, lty = 2, col = "firebrick", lwd = 1.5)
legend("topright", legend = colnames(w_means),
       fill = c("steelblue", "darkorange", "seagreen3"),
       bty = "n", title = "Confidence")
par(mar = c(5, 4, 4, 2))
