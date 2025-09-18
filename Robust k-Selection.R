# ============================
# Model-Selection on Speed Distributions (k = 2,3,4) — Robust Full Script
# ============================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(mclust)        # fast GMM & BIC (higher is better)
  library(mixtools)      # fallback EM (classic BIC lower is better)
})

# ----------------------------
# 0) INPUT & PARAMETERS
# ----------------------------
# Expect a data.frame/data.table 'black' with columns: MMSI, CALLSIGN, SPEED
# str(black)

# (Optional) known pelagic MMSIs — fill if pelagic-only mode is desired
pelagic_mmsi <- c(271073225,271072915,271072286,271062080,271073657,271072245,271073138,271072507,271072953,271073067,271073301,271072180,271072342,271072987,271072150,271072554,271073135,271073129,271069031,271073662,271072313,271069043,271073631,271073477,271073047,271072268,271072302,271072996,271069045,271073096,271072633,271072352,271072327,271073179,271073277,271073316,271073360,271072047,271069023,271062085,271073571,271072048,271073390,271073644,271072185,271073469,271072660,271072040,271073218,271073207,271073250,271072029,271072030,271056078,271072293,271073593,271072520,271073296,271073319)

use_pelagic_only <- TRUE   # FALSE => run on all TC-callsign vessels
min_n   <- 100             # minimum observations per vessel (50–200 reasonable)
max_n   <- 2000            # upper cap per vessel (for speed)
speed_min <- 0
speed_max <- 25
turkey_speed_cap <- 20     # upper speed threshold for initial filter (kn)

# (optional) signal requirement: at least 10 samples in 2–4 or 6–10 kn band
require_signal <- TRUE
signal_n_min   <- 10

# ----------------------------
# 1) DATA PREP
# ----------------------------
black$SPEED <- suppressWarnings(as.numeric(black$SPEED))
dt <- as.data.table(black)[
  grepl("^TC", CALLSIGN) & is.finite(SPEED) & SPEED >= speed_min & SPEED <= turkey_speed_cap,
  .(MMSI, SPEED)
]
if (use_pelagic_only) {
  dt <- dt[MMSI %in% pelagic_mmsi]
}

# Speed vectors per vessel
speeds_list <- split(dt$SPEED, dt$MMSI)
speeds_list <- lapply(speeds_list, function(x){
  x <- as.numeric(x); x <- x[is.finite(x)]
  x <- x[x >= speed_min & x <= speed_max]
  if (length(x) > max_n) x <- sample(x, max_n)
  x
})
# Minimum observation requirement
speeds_list <- speeds_list[sapply(speeds_list, length) >= min_n]

# Signal requirement (optional)
if (require_signal) {
  has_signal <- function(x){
    (sum(x >= 2 & x <= 4,  na.rm=TRUE) >= signal_n_min) ||
    (sum(x >= 6 & x <= 10, na.rm=TRUE) >= signal_n_min)
  }
  speeds_list <- speeds_list[sapply(speeds_list, has_signal)]
}

# ----------------------------
# 2) HELPERS
# ----------------------------
# Safe sd construction for scalar/vector 'sigmasq'
.make_sd <- function(sig2, k){
  if (is.null(sig2)) return(rep(NA_real_, k))
  s2 <- suppressWarnings(as.numeric(sig2))
  if (!length(s2) || any(!is.finite(s2))) return(rep(NA_real_, k))
  if (length(s2) == 1) rep(sqrt(s2), k) else sqrt(s2[seq_len(k)])
}

# Classic BIC (lower is better) — for mixtools
bic_from_loglik <- function(loglik, k, n){
  p <- 3*k - 1 # k means + k variances + (k-1) weights
  -2*loglik + p*log(n)
}

# Gaussian mixture density
mix_pdf <- function(x, mu, sd, pi){
  comp <- mapply(function(m, s, w) w * dnorm(x, mean=m, sd=s),
                 mu, sd, pi, SIMPLIFY = FALSE)
  rowSums(do.call(cbind, comp))
}

# Nearest towing component (prefer 2–4 kn; otherwise closest to 3 kn)
pick_towing_mu <- function(mu){
  mu <- as.numeric(mu)
  if (!length(mu)) return(NA_real_)
  in_band <- mu[mu >= 2 & mu <= 4]
  if (length(in_band)) in_band[which.min(abs(in_band - 3))]
  else mu[which.min(abs(mu - 3))]
}
pick_towing_idx <- function(mu){
  mu <- as.numeric(mu)
  if (!length(mu)) return(NA_integer_)
  in_band_idx <- which(mu >= 2 & mu <= 4)
  if (length(in_band_idx)) in_band_idx[which.min(abs(mu[in_band_idx] - 3))]
  else which.min(abs(mu - 3))
}

# Point label via posterior
label_obs <- function(x, pars){
  mu <- pars$mu; sd <- pars$sd; pi <- pars$pi
  if (any(is.null(list(mu,sd,pi))) || any(!lengths(list(mu,sd,pi)))) return(rep(NA_integer_, length(x)))
  K <- length(mu)
  dens <- sapply(seq_len(K), function(j) pi[j] * dnorm(x, mu[j], sd[j]))
  lab  <- max.col(dens, ties.method = "first")
  lab[!is.finite(rowSums(dens))] <- NA_integer_
  lab
}

# Pelagic candidate by component means (species-typical speed bands)
is_pelagic_by_means <- function(mu){
  mu <- as.numeric(mu)
  any(mu >= 2 & mu <= 4) && any(mu >= 6 & mu <= 10)
}

# ----------------------------
# 3) CORE FITTER: fit_k (mclust preferred, mixtools fallback)
# ----------------------------
fit_k <- function(x, k){
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (length(x) < 50 || sd(x) < 1e-6 || length(unique(x)) < 3){
    return(list(ok=FALSE, method=NA_character_, k=k,
                mu=NULL, sd=NULL, pi=NULL,
                score_unified=NA_real_, bic_native=NA_real_, bic_rule=NA_character_))
  }
  # 1) mclust (higher BIC = better) — reduce to single score
  fit_mc <- tryCatch(mclust::Mclust(x, G=k, modelNames=c("E","V")), error=function(e) NULL)
  if (!is.null(fit_mc)) {
    bic_val <- suppressWarnings(as.numeric(fit_mc$BIC))
    if (length(bic_val)) bic_val <- bic_val[1] else bic_val <- NA_real_
    if (is.finite(bic_val)) {
      mu  <- tryCatch(as.numeric(fit_mc$parameters$mean),  error=function(e) NULL)
      piw <- tryCatch(as.numeric(fit_mc$parameters$pro),   error=function(e) NULL)
      sig2<- tryCatch(fit_mc$parameters$variance$sigmasq,  error=function(e) NULL)
      sdv <- if (!is.null(mu)) .make_sd(sig2, length(mu)) else NULL
      return(list(
        ok = TRUE, method = "mclust", k = k,
        mu = mu, sd = sdv, pi = piw,
        score_unified = bic_val,        # unified: larger is better
        bic_native    = bic_val,
        bic_rule      = "higher_is_better"
      ))
    }
  }
  # 2) mixtools fallback (lower BIC = better) — unified score = -BIC
  fit_mx <- tryCatch(
    suppressWarnings(mixtools::normalmixEM(x, k=k, maxit=1000, verb=FALSE)),
    error=function(e) NULL
  )
  if (!is.null(fit_mx) && is.finite(fit_mx$loglik)) {
    bic <- bic_from_loglik(fit_mx$loglik, k, length(x))
    return(list(
      ok = TRUE, method = "mixtools", k = k,
      mu = as.numeric(fit_mx$mu),
      sd = as.numeric(fit_mx$sigma),
      pi = as.numeric(fit_mx$lambda),
      score_unified = -bic,             # unified scale
      bic_native    = bic,
      bic_rule      = "lower_is_better"
    ))
  }
  # 3) failure
  list(ok=FALSE, method=NA_character_, k=k,
       mu=NULL, sd=NULL, pi=NULL,
       score_unified=NA_real_, bic_native=NA_real_, bic_rule=NA_character_)
}

fit_both_k234 <- function(x){
  list(k2 = fit_k(x,2),
       k3 = fit_k(x,3),
       k4 = fit_k(x,4))
}

# Evidence categories (Kass & Raftery BIC-difference scale)
evidence_label <- function(delta){
  if (!is.finite(delta)) return(NA_character_)
  if (delta < 0) return(NA_character_)
  if (delta < 2)  "weak"
  else if (delta < 6)  "positive"
  else if (delta < 10) "strong"
  else                 "very strong"
}

# ----------------------------
# 4) FIT ALL VESSELS
# ----------------------------
mmsi_ids <- names(speeds_list)
rows <- vector("list", length(mmsi_ids))

for (i in seq_along(mmsi_ids)){
  id <- mmsi_ids[[i]]
  x  <- speeds_list[[id]]

  fits <- fit_both_k234(x)

  sc2 <- if (isTRUE(fits$k2$ok)) fits$k2$score_unified else NA_real_
  sc3 <- if (isTRUE(fits$k3$ok)) fits$k3$score_unified else NA_real_
  sc4 <- if (isTRUE(fits$k4$ok)) fits$k4$score_unified else NA_real_

  scores <- c(k2=sc2, k3=sc3, k4=sc4)
  if (!any(is.finite(scores))) {
    best_k <- NA_integer_; tow_mu <- NA_real_; delta_best2 <- NA_real_; best_label <- NA_character_
  } else {
    ord <- order(scores, decreasing = TRUE, na.last = NA)
    best_name <- names(scores)[ord[1]]
    best_k    <- switch(best_name, k2=2L, k3=3L, k4=4L)
    # best - second best (unified BIC difference)
    delta_best2 <- if (length(ord) >= 2) (scores[ord[1]] - scores[ord[2]]) else NA_real_
    best_label  <- paste(evidence_label(delta_best2), best_name)
    # towing mean (from winner's component means)
    best_mu <- switch(best_name,
                      k2 = fits$k2$mu,
                      k3 = fits$k3$mu,
                      k4 = fits$k4$mu)
    tow_mu <- if (length(best_mu)) pick_towing_mu(best_mu) else NA_real_
  }

  rows[[i]] <- data.frame(
    MMSI = as.numeric(id),
    score_k2 = sc2,
    score_k3 = sc3,
    score_k4 = sc4,
    best_k   = best_k,
    delta_best_vs_second = delta_best2,
    evidence_best = best_label,
    towing_mu_best = tow_mu
  )
}

res_dt <- data.table::rbindlist(rows, fill=TRUE)

# ----------------------------
# 5) SUMMARIES
# ----------------------------
cat("Vessel count:", nrow(res_dt), "\n")
print(table(res_dt$best_k, useNA="ifany"))

# Evidence distribution (by best model)
print(table(res_dt$evidence_best, useNA="ifany"))

# Towing mu summaries
summary_tbl <- res_dt[, .(
  vessels        = .N,
  n_k2           = sum(best_k==2, na.rm=TRUE),
  n_k3           = sum(best_k==3, na.rm=TRUE),
  n_k4           = sum(best_k==4, na.rm=TRUE),
  median_tow_mu  = median(towing_mu_best, na.rm=TRUE),
  q25_tow_mu     = quantile(towing_mu_best, 0.25, na.rm=TRUE),
  q75_tow_mu     = quantile(towing_mu_best, 0.75, na.rm=TRUE)
)]
print(summary_tbl)

# ----------------------------
# 6) REPRESENTATIVE FIT PLOTS (optional)
# ----------------------------
# Pick examples for three categories: very strong k=3, any k=2, very strong k=4 (if present)
pick_one <- function(res, key){ res[grep(key, res$evidence_best), ][1L]$MMSI }
id_k3_vs <- pick_one(res_dt, "^very strong k=3$")
id_k2_bt <- pick_one(res_dt, "^.*k=2$")
id_k4_vs <- pick_one(res_dt, "^very strong k=4$")

pick_ids <- unique(na.omit(c(id_k3_vs, id_k2_bt, id_k4_vs)))

if (length(pick_ids)) {
  op <- par(mfrow=c(length(pick_ids),1), mar=c(4,4,3,1))
  for (id in pick_ids) {
    x <- as.numeric(speeds_list[[as.character(id)]])
    x <- x[is.finite(x)]
    if (length(x) < 50) next
    f2 <- fit_k(x,2); f3 <- fit_k(x,3); f4 <- fit_k(x,4)
    hist(x, breaks=60, freq=FALSE, main=paste("MMSI", id, "- k=2/3/4 fits"), xlab="Speed (kn)")
    xs <- seq(0, max(12, quantile(x, 0.99, na.rm=TRUE)), length.out=400)
    if (isTRUE(f2$ok)) lines(xs, mix_pdf(xs, f2$mu, f2$sd, f2$pi), lwd=2, lty=2)
    if (isTRUE(f3$ok)) lines(xs, mix_pdf(xs, f3$mu, f3$sd, f3$pi), lwd=2, lty=1)
    if (isTRUE(f4$ok)) lines(xs, mix_pdf(xs, f4$mu, f4$sd, f4$pi), lwd=2, lty=3)
    legend("topright", c("k=2","k=3","k=4"), lty=c(2,1,3), lwd=2, bty="n")
  }
  par(op)
}

# ----------------------------
# 7) DOWNSTREAM: towing share and pelagic candidate counts (k2/k3/k4)
# ----------------------------
calc_downstream <- function(MMSI_vec){
  rows <- vector("list", length(MMSI_vec))
  for (i in seq_along(MMSI_vec)) {
    id <- MMSI_vec[i]
    x  <- as.numeric(speeds_list[[as.character(id)]])
    x  <- x[is.finite(x)]
    if (length(x) < 50) next
    f2 <- fit_k(x,2); f3 <- fit_k(x,3); f4 <- fit_k(x,4)

    # k=2
    tow_share_k2 <- pelagic_k2 <- NA
    if (isTRUE(f2$ok)) {
      cls2 <- label_obs(x, f2)
      tow_idx2 <- pick_towing_idx(f2$mu)
      tow_share_k2 <- mean(cls2 == tow_idx2, na.rm=TRUE)
      pelagic_k2   <- is_pelagic_by_means(f2$mu)
    }
    # k=3
    tow_share_k3 <- pelagic_k3 <- NA
    if (isTRUE(f3$ok)) {
      cls3 <- label_obs(x, f3)
      tow_idx3 <- pick_towing_idx(f3$mu)
      tow_share_k3 <- mean(cls3 == tow_idx3, na.rm=TRUE)
      pelagic_k3   <- is_pelagic_by_means(f3$mu)
    }
    # k=4
    tow_share_k4 <- pelagic_k4 <- NA
    if (isTRUE(f4$ok)) {
      cls4 <- label_obs(x, f4)
      tow_idx4 <- pick_towing_idx(f4$mu)
      tow_share_k4 <- mean(cls4 == tow_idx4, na.rm=TRUE)
      pelagic_k4   <- is_pelagic_by_means(f4$mu)
    }

    rows[[i]] <- data.frame(
      MMSI = as.numeric(id), n_obs = length(x),
      tow_share_k2 = tow_share_k2, tow_share_k3 = tow_share_k3, tow_share_k4 = tow_share_k4,
      pelagic_k2 = pelagic_k2, pelagic_k3 = pelagic_k3, pelagic_k4 = pelagic_k4
    )
  }
  df <- data.table::rbindlist(rows, fill=TRUE)
  list(
    per_vessel = df,
    summary = df[, .(
      vessels = .N,
      pelagic_k2_n = sum(pelagic_k2 %in% TRUE, na.rm=TRUE),
      pelagic_k3_n = sum(pelagic_k3 %in% TRUE, na.rm=TRUE),
      pelagic_k4_n = sum(pelagic_k4 %in% TRUE, na.rm=TRUE),
      tow_share_k2_med = median(tow_share_k2, na.rm=TRUE),
      tow_share_k3_med = median(tow_share_k3, na.rm=TRUE),
      tow_share_k4_med = median(tow_share_k4, na.rm=TRUE)
    )]
  )
}

down <- calc_downstream(as.numeric(names(speeds_list)))
print(down$summary)

# ----------------------------
# 8) (Optional) SAVE OUTPUTS
# ----------------------------
# data.table::fwrite(res_dt, "model_selection_k234_results.csv")
# data.table::fwrite(down$per_vessel, "downstream_per_vessel.csv")
# data.table::fwrite(down$summary, "downstream_summary.csv")

# how many vessels 'won' which k?
cat("\nCounts by best_k:\n")
print(table(res_dt$best_k, useNA="ifany"))

# evidence strength distribution (Kass & Raftery thresholds)
cat("\nEvidence (best vs second-best):\n")
print(sort(table(res_dt$evidence_best), decreasing = TRUE))

# towing_mu (mean 'towing speed' for the winning model) — expected band 2–4 kn
cat("\nTowing mu summary:\n")
print(summary(res_dt$towing_mu_best))
hist(res_dt$delta_best_vs_second,
     breaks = 40,
     main   = "ΔBIC (best − second best): larger = stronger evidence",
     xlab   = "ΔBIC (unified scale)")
abline(v = c(2,6,10), lty = c(2,3,3))  # weak / strong / very strong thresholds

boxplot(towing_mu_best ~ factor(best_k), data = res_dt,
        xlab = "Chosen k", ylab = "Towing mean (kn)",
        main = "Towing mean by chosen k")
abline(h = c(2,4), lty = 2)

cat("\nDownstream summary (k2/k3/k4):\n")
print(down$summary)

down_tbl <- with(down$summary, data.frame(
  model         = c("k=2","k=3","k=4"),
  pelagic_n     = c(pelagic_k2_n, pelagic_k3_n, pelagic_k4_n),
  tow_share_med = round(c(tow_share_k2_med, tow_share_k3_med, tow_share_k4_med), 3)
))
down_tbl
