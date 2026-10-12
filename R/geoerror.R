# =============================================================================
# geoerror.R  Correcting area-based index estimates for postal-code geocoding error
#
# Method (see data/derived/err/METHOD.md for the full description and references)
#   1. Calibration. From the NAR validation (scripts/err1_transitions.R), each address
#      has a true DA (reference location) and a PCCF+ DA. Merging any DA-level index
#      gives, for each stratum s, counts n_s(t, a) of addresses with true quintile t
#      and assigned quintile a. The classification matrix P_s(a | t) is estimated with
#      a Dirichlet posterior that shrinks sparse strata toward the national matrix for
#      the same SAC group (and delivery mode): alpha_s(t, .) = [n_s(t, .) + k * P_parent(. | t)] / deff_s(t),
#      where deff_s(t) is the design effect of clustering within postal codes, estimated by a
#      Bayesian bootstrap over postal codes (optional `clusters` input).
#   2. Correction (matrix method, estimated by EM). For a user's counts o_s(a) by
#      stratum and assigned quintile, the true-quintile counts c_s(t) maximize the
#      multinomial likelihood of o_s given P_s(a | t); EM keeps them non-negative. EM runs
#      to convergence (cap 5,000 iterations); convergence and boundary solutions are reported.
#      Events and population are corrected separately, so rates and rate ratios by
#      true quintile follow. Assumption: within a stratum, misclassification does not
#      depend on the outcome (nondifferential) and P_s(a | t) transports to the user's data.
#   3. Uncertainty. Monte Carlo over (i) Dirichlet draws of P_s(a | t) and (ii) Poisson
#      sampling of observed event counts (and population counts when the population is itself
#      a sample, sample_n = TRUE); percentile simulation intervals. P_s is informed only by the
#      register, so this is a modular ("cut") posterior combined with an ML fit per draw.
#   4. Continuous scores. Stratum-specific regression calibration E[true | assigned]
#      with the residual standard deviation as the per-record error term.
#
# Main functions
#   ge_build_calibration(tr, index, strata, k)   -> calibration object
#   ge_correct(obs, cal, draws, ...)             -> corrected distribution, rates, rate ratios
#   ge_calibrate_score(df, cal, score_col)       -> calibrated score and error SD per record
#
# Strata names used in the tables: prov (PCCF+ province code), sac (PCCF+ SAC group),
# dmt3 (PCCF+ delivery mode: "Urban LC", "Rural", "Other"). The user's data must carry
# the same columns, taken from their own PCCF+ output.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

GE_K <- 5L   # number of index categories (quintiles)

# ---- 1. Calibration ---------------------------------------------------------------

#' @param tr     DA transition table (err1): strata columns, ref_DAUID, pccf_DAUID, n
#' @param index  data.table(DAUID, q [1..5], score [optional]) for the index to calibrate
#' @param strata character vector of stratum columns, e.g. c("prov", "sac") or c("prov", "sac", "dmt3")
#' @param k      prior weight (pseudo-addresses per true-quintile row) toward the parent matrix
#' @param clusters optional postal-code-level transitions (err1 pc_transitions.rds, same rows
#'               as tr): strata columns, PCODE, ref_DAUID, pccf_DAUID, n. Enables the design effect.
#' @param bb_reps Bayesian bootstrap replicates for the design effect
ge_build_calibration <- function(tr, index, strata = c("prov", "sac"), k = 20,
                                 clusters = NULL, bb_reps = 200) {
  stopifnot(all(c("DAUID", "q") %in% names(index)))
  ix <- unique(as.data.table(index), by = "DAUID")
  has_score <- "score" %in% names(ix)
  d <- copy(tr)
  d[ix, `:=`(t = i.q, score_t = if (has_score) i.score else NA_real_), on = .(ref_DAUID = DAUID)]
  d[ix, `:=`(a = i.q, score_a = if (has_score) i.score else NA_real_), on = .(pccf_DAUID = DAUID)]
  d <- d[!is.na(t) & !is.na(a)]

  # Parent (national) matrices used for shrinkage and as fallbacks for strata absent from
  # the calibration: drop province, then SAC, then everything.
  parent <- setdiff(strata, "prov")
  levels_list <- unique(list(parent, setdiff(strata, c("prov", "sac")), character(0)))
  parent_matrix <- function(lv) {
    pc <- d[, .(n = sum(n)), by = c(lv, "t", "a")]
    g <- CJ(t = 1:GE_K, a = 1:GE_K)
    if (length(lv)) g <- unique(pc[, ..lv])[, as.list(g), by = lv]
    pc <- pc[g, on = c(lv, "t", "a")][is.na(n), n := 0]
    pc[, p_parent := (n + 0.5) / sum(n + 0.5), by = c(lv, "t")][]   # Jeffreys-smoothed
  }
  parents <- lapply(levels_list, parent_matrix)
  pcnt <- parents[[1]]
  cnt  <- d[, .(n = sum(n)), by = c(strata, "t", "a")]

  sgrid <- unique(cnt[, ..strata])[, as.list(CJ(t = 1:GE_K, a = 1:GE_K)), by = strata]
  cl <- cnt[sgrid, on = c(strata, "t", "a")][is.na(n), n := 0]
  setorderv(cl, c(strata, "t", "a"))
  cl[pcnt, p_parent := i.p_parent, on = c(parent, "t", "a")]
  cl[, n_t := sum(n), by = c(strata, "t")]

  # Design effect of clustering within postal codes. Addresses in one postal code share
  # their geocoding error, so the counts carry less information than independent draws.
  # A Bayesian bootstrap over postal codes (Rubin 1981) gives the variance of each row
  # P_s(a | t); deff_s(t) is its ratio to the variance under independent addresses. The
  # Dirichlet parameters are divided by deff, which keeps the posterior mean and inflates
  # its variance by the design effect.
  cl[, deff := 1]
  if (!is.null(clusters)) {
    cd <- as.data.table(clusters)[, c(strata, "PCODE", "ref_DAUID", "pccf_DAUID", "n"), with = FALSE]
    cd[ix, t := i.q, on = .(ref_DAUID = DAUID)]; cd[ix, a := i.q, on = .(pccf_DAUID = DAUID)]
    cd <- cd[!is.na(t) & !is.na(a), .(n = sum(n)), by = c(strata, "PCODE", "t", "a")]
    cl[, cell := .I]; cl[, row := .GRP, by = c(strata, "t")]
    cd[cl, cell := i.cell, on = c(strata, "t", "a")]
    pid <- as.integer(factor(cd$PCODE)); np <- max(pid)
    Bw <- matrix(0, nrow(cl), bb_reps)
    for (b in seq_len(bb_reps)) {
      v <- rowsum(cd$n * rexp(np)[pid], cd$cell)
      Bw[as.integer(rownames(v)), b] <- v
    }
    rt <- rowsum(Bw, cl$row)[cl$row, , drop = FALSE]
    Pb <- ifelse(rt > 0, Bw / rt, NA_real_)
    cl[, var_bb := apply(Pb, 1, var, na.rm = TRUE)]
    cl[, p_hat := fifelse(n_t > 0, n / n_t, 0)]
    cl[, deff := {
      num <- sum(var_bb, na.rm = TRUE); den <- sum(p_hat * (1 - p_hat)) / (n_t[1] + 1)
      if (den > 0 && is.finite(num)) min(max(1, num / den), n_t[1] + k) else 1 }, by = c(strata, "t")]
    cl[, c("cell", "row", "var_bb", "p_hat") := NULL]
  }
  cl[, alpha := (n + k * p_parent) / deff]          # mean-preserving: variance inflated by deff
  cl[, p_a_given_t := alpha / sum(alpha), by = c(strata, "t")]
  cl[, p_t_given_a := (n_t * p_a_given_t) / sum(n_t * p_a_given_t), by = c(strata, "a")]

  summ <- d[, .(N = sum(n),
                pct_misclassified = 100 * sum(n[t != a]) / sum(n),
                pct_off_2plus = 100 * sum(n[abs(t - a) >= 2]) / sum(n)), by = strata]
  kap <- cl[, .(kappa_linear = {
    m <- matrix(0, GE_K, GE_K); m[cbind(t, a)] <- n; N <- sum(m)
    w <- 1 - abs(outer(1:GE_K, 1:GE_K, "-")) / (GE_K - 1)
    if (N == 0) NA_real_ else {
      po <- sum(w * m) / N; pe <- sum(w * outer(rowSums(m), colSums(m))) / N^2
      (po - pe) / (1 - pe) } }), by = strata]
  summ <- summ[kap, on = strata]

  score <- NULL
  if (has_score) {
    s <- d[!is.na(score_t) & !is.na(score_a)]
    score <- s[, {
      w <- n; W <- sum(w)
      ma <- sum(w * score_a) / W; mt <- sum(w * score_t) / W
      va <- sum(w * (score_a - ma)^2) / W; vt <- sum(w * (score_t - mt)^2) / W
      cv <- sum(w * (score_a - ma) * (score_t - mt)) / W
      e  <- score_a - score_t
      slope <- if (va > 0) cv / va else NA_real_
      .(N = W, mean_assigned = ma, mean_true = mt, sd_assigned = sqrt(va), sd_true = sqrt(vt),
        correlation = if (va > 0 && vt > 0) cv / sqrt(va * vt) else NA_real_,
        error_bias = sum(w * e) / W, error_sd = sqrt(sum(w * (e - sum(w * e) / W)^2) / W),
        calib_slope = slope, calib_intercept = mt - slope * ma,
        calib_resid_sd = sqrt(pmax(vt - slope * cv, 0)))
    }, by = strata]
  }
  structure(list(class = cl, parent = pcnt, parents = parents, parent_levels = levels_list,
                 summary = summ, score = score, strata = strata, parent_strata = parent, k = k),
            class = "ge_calibration")
}

# ---- 2. EM for the matrix method, vectorized over strata x draws ---------------------

# o: matrix [R, K] observed counts by assigned category; P: array [R, K(t), K(a)] of P(a|t).
# Returns matrix [R, K] of estimated true counts (same total as o in each row), with
# attributes: converged (logical per row: max change in pi below tol before iter) and
# boundary (logical per row: some true category ends with share < bnd, i.e. the direct
# inversion of the matrix relation would give a negative count).
ge_em <- function(o, P, iter = 5000, tol = 1e-8, bnd = 1e-4) {
  R <- nrow(o); K <- ncol(o); N <- rowSums(o)
  em_step <- function(p, oo, PP, NN) {                    # one EM update for a set of rows
    q <- matrix(0, nrow(p), K)                            # q[r, a] = sum_t pi[r, t] P[r, t, a]
    for (a in 1:K) for (t in 1:K) q[, a] <- q[, a] + p[, t] * PP[, t, a]
    ratio <- ifelse(q > 0, oo / q, 0)
    new <- matrix(0, nrow(p), K)
    for (t in 1:K) for (a in 1:K) new[, t] <- new[, t] + ratio[, a] * PP[, t, a]
    p * new / NN
  }
  # SQUAREM acceleration (Varadhan and Roland 2008), vectorized over rows; converged rows
  # (largest change of a plain EM step below tol) leave the active set. An extrapolated
  # point outside the simplex falls back to two plain EM steps. `iter` counts EM steps.
  pi <- matrix(1 / K, R, K); conv <- N == 0; act <- which(!conv); steps <- 0
  while (length(act) && steps < iter) {
    p0 <- pi[act, , drop = FALSE]; oo <- o[act, , drop = FALSE]; PP <- P[act, , , drop = FALSE]; NN <- N[act]
    p1 <- em_step(p0, oo, PP, NN); r <- p1 - p0
    done <- apply(abs(r), 1, max) < tol
    pi[act[done], ] <- p1[done, , drop = FALSE]; conv[act[done]] <- TRUE
    go <- !done
    if (any(go)) {
      p0 <- p0[go, , drop = FALSE]; p1 <- p1[go, , drop = FALSE]; r <- r[go, , drop = FALSE]
      oo <- oo[go, , drop = FALSE]; PP <- PP[go, , , drop = FALSE]; NN <- NN[go]
      p2 <- em_step(p1, oo, PP, NN); v <- p2 - p1 - r
      al <- pmin(-sqrt(rowSums(r^2)) / pmax(sqrt(rowSums(v^2)), 1e-300), -1)
      pn <- p0 - 2 * al * r + al^2 * v
      bad <- !is.finite(rowSums(pn)) | apply(pn, 1, min) < 0
      pn[bad, ] <- p2[bad, ]
      pi[act[go], ] <- em_step(pn / rowSums(pn), oo, PP, NN)
    }
    act <- act[go]; steps <- steps + 3
  }
  out <- pi * N
  attr(out, "converged") <- conv
  attr(out, "boundary") <- N > 0 & apply(pi, 1, min) < bnd
  out
}

# Dirichlet draws of P(a|t): returns array [S*D, K, K] (stratum-major), from alpha [S, K, K]
ge_draw_P <- function(alpha, draws) {
  S <- dim(alpha)[1]; K <- dim(alpha)[2]
  out <- array(0, c(S * draws, K, K))
  idx <- rep(seq_len(S), each = draws)
  for (t in 1:K) {
    al <- alpha[idx, t, , drop = FALSE][, 1, ]
    g <- matrix(rgamma(S * draws * K, shape = as.vector(t(al))), ncol = K, byrow = TRUE)
    z <- rowSums(g) == 0                           # all-zero draw (only possible for tiny shapes)
    g[z, ] <- al[z, , drop = FALSE]
    out[, t, ] <- g / rowSums(g)
  }
  out
}

# Alpha array [S, K, K] for the user's strata; strata missing from the calibration fall
# back to the parent (national) matrix with the prior weight k.
ge_alpha <- function(cal, keys) {
  st <- cal$strata; K <- GE_K
  A <- array(NA_real_, c(nrow(keys), K, K))
  for (i in seq_len(nrow(keys))) {
    sub <- cal$class[keys[i], on = st, nomatch = NULL]
    if (nrow(sub) == K * K) {
      A[i, , ] <- matrix(sub$alpha, K, K, byrow = TRUE)           # rows t, columns a
    } else {                                                        # hierarchical fallback
      for (j in seq_along(cal$parent_levels)) {
        lv <- cal$parent_levels[[j]]
        ps <- if (length(lv)) cal$parents[[j]][keys[i, lv, with = FALSE], on = lv, nomatch = NULL] else cal$parents[[j]]
        if (nrow(ps) == K * K) break
      }
      if (nrow(ps) != K * K) stop("No calibration for stratum: ", paste(unlist(keys[i]), collapse = " / "))
      A[i, , ] <- matrix(cal$k * ps$p_parent, K, K, byrow = TRUE)
    }
  }
  A
}

# ---- 3. Correction with confidence intervals ----------------------------------------

#' @param obs  data.table with the calibration's stratum columns, `a` (assigned quintile 1..5),
#'             `n` (population or person-time) and optionally `events`
#' @param cal  ge_calibration object
#' @param draws Monte Carlo draws; @param conf confidence level
#' @param sample_events Poisson-resample event counts (sampling error of rates)
#' @param sample_n      Poisson-resample population counts (usually FALSE for census-type denominators)
#' @return list(distribution, rates, rate_ratios, by_stratum) with naive and corrected estimates
ge_correct <- function(obs, cal, draws = 1000, conf = 0.95, sample_events = TRUE,
                       sample_n = FALSE, ref_q = 1L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  st <- cal$strata; K <- GE_K
  obs <- as.data.table(obs)
  has_ev <- "events" %in% names(obs)
  agg <- obs[, .(n = sum(n), events = if (has_ev) sum(events) else 0), by = c(st, "a")]
  keys <- unique(agg[, ..st]); S <- nrow(keys)
  keys[, sid := .I]
  agg[keys, sid := i.sid, on = st]
  O_n <- matrix(0, S, K); O_e <- matrix(0, S, K)
  O_n[cbind(agg$sid, agg$a)] <- agg$n; O_e[cbind(agg$sid, agg$a)] <- agg$events

  A <- ge_alpha(cal, keys[, ..st])
  Pmean <- A / array(apply(A, c(1, 2), sum), dim(A))           # posterior mean P(a|t)

  # point estimates
  C_n <- ge_em(O_n, Pmean); C_e <- if (has_ev) ge_em(O_e, Pmean) else NULL

  # Monte Carlo
  Pd <- ge_draw_P(A, draws)
  idx <- rep(seq_len(S), each = draws)
  On_d <- O_n[idx, , drop = FALSE]; Oe_d <- O_e[idx, , drop = FALSE]
  if (sample_n) On_d[] <- rpois(length(On_d), On_d)
  if (has_ev && sample_events) Oe_d[] <- rpois(length(Oe_d), Oe_d)
  Cn_d <- ge_em(On_d, Pd); Ce_d <- if (has_ev) ge_em(Oe_d, Pd) else NULL
  # diagnostics over stratum x draw rows: EM convergence, and boundary solutions (a check
  # on whether the prior for P is compatible with the data; MacLehose and Gustafson 2012)
  wn <- rowSums(On_d)
  diag <- data.table(
    em_converged = mean(c(attr(Cn_d, "converged"), if (has_ev) attr(Ce_d, "converged"))),
    boundary_n = mean(attr(Cn_d, "boundary")[wn > 0]),
    boundary_n_pop = sum(wn[attr(Cn_d, "boundary")]) / sum(wn),
    boundary_events = if (has_ev) mean(attr(Ce_d, "boundary")[rowSums(Oe_d) > 0]) else NA_real_)
  dr <- rep(seq_len(draws), times = S)                             # draw id for each row
  tot_n <- rowsum(Cn_d, dr); tot_e <- if (has_ev) rowsum(Ce_d, dr) else NULL   # [draws, K]

  lo <- (1 - conf) / 2; hi <- 1 - lo
  qs <- function(m) t(apply(m, 2, quantile, c(lo, hi), na.rm = TRUE))

  # distribution
  naive_n <- colSums(O_n); corr_n <- colSums(C_n)
  share_d <- 100 * tot_n / rowSums(tot_n)
  dist <- data.table(q = 1:K, naive_n = naive_n, naive_pct = 100 * naive_n / sum(naive_n),
                     corrected_n = corr_n, corrected_pct = 100 * corr_n / sum(corr_n),
                     se_pct = apply(share_d, 2, sd))
  dist[, c("lo_pct", "hi_pct") := as.data.table(qs(share_d))]

  rates <- rr <- NULL
  if (has_ev) {
    naive_e <- colSums(O_e); corr_e <- colSums(C_e)
    rate_d <- tot_e / tot_n
    rates <- data.table(q = 1:K, naive_events = naive_e, naive_rate = naive_e / naive_n,
                        corrected_events = corr_e, corrected_rate = corr_e / corr_n,
                        se_rate = apply(rate_d, 2, sd))
    rates[, c("lo_rate", "hi_rate") := as.data.table(qs(rate_d))]
    # naive Poisson CI for the rate ratio vs reference quintile
    z <- qnorm(hi)
    lrr_d <- log(rate_d / rate_d[, ref_q])
    rr <- data.table(q = 1:K,
                     naive_rr = (naive_e / naive_n) / (naive_e[ref_q] / naive_n[ref_q]),
                     corrected_rr = (corr_e / corr_n) / (corr_e[ref_q] / corr_n[ref_q]),
                     se_log_rr = apply(lrr_d, 2, sd))
    se_naive <- sqrt(1 / naive_e + 1 / naive_e[ref_q]); se_naive[ref_q] <- 0
    rr[, `:=`(naive_lo = naive_rr * exp(-z * se_naive), naive_hi = naive_rr * exp(z * se_naive))]
    rr[, c("lo_rr", "hi_rr") := as.data.table(exp(qs(lrr_d)))]
  }
  by_stratum <- cbind(keys[rep(seq_len(S), each = K), ..st], q = rep(1:K, S),
                      observed_n = as.vector(t(O_n)), corrected_n = as.vector(t(C_n)))
  list(distribution = dist[], rates = rates, rate_ratios = rr, by_stratum = by_stratum,
       diagnostics = diag, draws = draws, conf = conf)
}

# ---- 4. Continuous scores: regression calibration -------------------------------------

#' Adds calibrated_score (E[true | assigned, stratum]) and error_sd (residual SD) to df.
ge_calibrate_score <- function(df, cal, score_col) {
  stopifnot(!is.null(cal$score))
  df <- as.data.table(copy(df))
  df[cal$score, `:=`(slope = i.calib_slope, intercept = i.calib_intercept, error_sd = i.calib_resid_sd),
     on = cal$strata]
  df[, calibrated_score := intercept + slope * get(score_col)]
  df[, c("slope", "intercept") := NULL][]
}

# =============================================================================
# 5. Postal-code-level method (record-level data with postal codes)
#
# The NAR is a census of residential addresses, so for every postal code p the
# distribution of true index quintiles among its addresses, pi_p(t), is known. A person
# whose only location is p has true quintile t with probability pi_p(t) (assuming residents
# are spread over the code's addresses as the addresses are). This uses the actual
# geography of each postal code rather than a stratum average.
#   Distribution: expected counts sum_i pi_i(t); intervals by multiple imputation of t.
#   Rates by quintile: binomial model P(Y = 1 | t) = lambda_t with t latent and prior pi_i(t),
#   fitted by Gibbs sampling (data augmentation; Greenland 2009). Cases are allocated
#   toward quintiles with higher lambda_t, which removes the attenuation that allocation
#   by pi alone would leave. Nondifferential assumption: Y depends on location only
#   through t.
# Postal codes absent from the table fall back to the user's assigned quintile.
# =============================================================================

#' @param addr  data.table(PCODE, ref_DAUID, n): residential address counts by postal code and true DA
#' @param index data.table(DAUID, q, score [optional])
#' @return data.table(PCODE, n_addr, p1..p5, mean_score, sd_score)
ge_pc_table <- function(addr, index) {
  ix <- unique(as.data.table(index), by = "DAUID")
  d <- copy(as.data.table(addr))
  d[ix, `:=`(q = i.q, score = if ("score" %in% names(ix)) i.score else NA_real_), on = .(ref_DAUID = DAUID)]
  d <- d[!is.na(q)]
  w <- dcast(d[, .(n = sum(n)), by = .(PCODE, q)], PCODE ~ q, value.var = "n", fill = 0)
  for (k in as.character(1:GE_K)) if (!k %in% names(w)) w[, (k) := 0]
  w[, n_addr := rowSums(.SD), .SDcols = as.character(1:GE_K)]
  for (k in 1:GE_K) w[, (paste0("p", k)) := get(as.character(k)) / n_addr]
  w[, (as.character(1:GE_K)) := NULL]
  sc <- d[!is.na(score), .(mean_score = sum(n * score) / sum(n),
                           sd_score = sqrt(pmax(sum(n * score^2) / sum(n) - (sum(n * score) / sum(n))^2, 0))),
          by = PCODE]
  sc[w, on = "PCODE"]
}

# Prior matrix [N, K] for records with postal code `pcode` and assigned quintile `a`
ge_pc_prior <- function(pcode, a, pc_table) {
  P <- matrix(0, length(pcode), GE_K)
  m <- match(pcode, pc_table$PCODE)
  hit <- !is.na(m)
  P[hit, ] <- as.matrix(pc_table[m[hit], paste0("p", 1:GE_K), with = FALSE])
  P[cbind(which(!hit), a[!hit])] <- 1                          # fallback: assigned quintile
  attr(P, "fallback_share") <- mean(!hit)
  P
}

# Vectorized multinomial allocation of counts m (length G) with probabilities W [G, K]
ge_rmultinom_rows <- function(m, W) {
  G <- length(m); K <- ncol(W); out <- matrix(0L, G, K)
  rem <- m; prem <- rep(1, G)
  for (k in 1:(K - 1)) {
    pk <- ifelse(prem > 0, pmin(pmax(W[, k] / prem, 0), 1), 0)
    out[, k] <- rbinom(G, rem, pk)
    rem <- rem - out[, k]; prem <- prem - W[, k]
  }
  out[, K] <- rem
  out
}

#' Distribution, rates and rate ratios by true quintile from record-level postal codes.
#' @param pcode, a  postal code and assigned quintile per record
#' @param y        0/1 outcome per record (NULL for distribution only)
#' @param pc_table output of ge_pc_table
#' @param iter, burn Gibbs iterations; @param prior Beta(a, b) prior for lambda_t
#' @param return_draws also return the retained draws of lambda (for convergence diagnostics)
#' @param lam_init starting rates by quintile (default: overall rate); dispersed starts for multiple chains
ge_pc_correct <- function(pcode, a, pc_table, y = NULL, iter = 1500, burn = 500, conf = 0.95,
                          prior = c(0.5, 0.5), ref_q = 1L, seed = NULL, return_draws = FALSE,
                          lam_init = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- GE_K
  P <- ge_pc_prior(pcode, a, pc_table)
  fb <- attr(P, "fallback_share")
  lo <- (1 - conf) / 2; hi <- 1 - lo
  # group identical (prior, outcome) patterns: records in the same postal code share pi
  grp <- if (is.null(y)) paste(pcode, a) else paste(pcode, a, y)
  g <- data.table(grp = grp, i = seq_along(grp))[, .(m = .N, i = i[1]), by = grp]
  W0 <- P[g$i, , drop = FALSE]; m <- g$m
  yy <- if (is.null(y)) rep(0L, nrow(g)) else y[g$i]

  expected <- colSums(P)
  if (is.null(y)) {                                          # distribution only: impute t
    sims <- t(replicate(iter - burn, colSums(ge_rmultinom_rows(m, W0))))
    sh <- 100 * sims / rowSums(sims)
    return(list(distribution = data.table(q = 1:K, naive_pct = 100 * tabulate(a, K) / length(a),
                                          corrected_pct = 100 * expected / sum(expected),
                                          se_pct = apply(sh, 2, sd),
                                          lo_pct = apply(sh, 2, quantile, lo), hi_pct = apply(sh, 2, quantile, hi)),
                fallback_share = fb))
  }
  lam <- if (is.null(lam_init)) rep(mean(y), K) else lam_init
  keep <- matrix(NA_real_, iter - burn, K); shares <- keep
  for (it in seq_len(iter)) {
    W <- W0 * outer(ifelse(yy == 1, 1, 0), lam) + W0 * outer(ifelse(yy == 1, 0, 1), 1 - lam)
    W <- W / rowSums(W)
    alloc <- ge_rmultinom_rows(m, W)
    cases <- colSums(alloc * yy); tot <- colSums(alloc)
    lam <- rbeta(K, prior[1] + cases, prior[2] + tot - cases)
    if (it > burn) { keep[it - burn, ] <- lam; shares[it - burn, ] <- 100 * tot / sum(tot) }
  }
  lrr <- log(keep / keep[, ref_q])
  naive <- data.table(q = 1:K)[, `:=`(n = tabulate(a, K), e = tabulate(a[y == 1], K))]
  se_naive <- sqrt(1 / naive$e + 1 / naive$e[ref_q]); se_naive[ref_q] <- 0
  naive_rr <- (naive$e / naive$n) / (naive$e[ref_q] / naive$n[ref_q])
  z <- qnorm(hi)
  list(distribution = data.table(q = 1:K, naive_pct = 100 * naive$n / sum(naive$n),
                                 corrected_pct = colMeans(shares), se_pct = apply(shares, 2, sd),
                                 lo_pct = apply(shares, 2, quantile, lo), hi_pct = apply(shares, 2, quantile, hi)),
       rates = data.table(q = 1:K, naive_rate = naive$e / naive$n, corrected_rate = colMeans(keep),
                          lo_rate = apply(keep, 2, quantile, lo), hi_rate = apply(keep, 2, quantile, hi)),
       rate_ratios = data.table(q = 1:K, naive_rr = naive_rr, naive_lo = naive_rr * exp(-z * se_naive),
                                naive_hi = naive_rr * exp(z * se_naive),
                                corrected_rr = exp(colMeans(lrr)), se_log_rr = apply(lrr, 2, sd),
                                lo_rr = exp(apply(lrr, 2, quantile, lo)), hi_rr = exp(apply(lrr, 2, quantile, hi))),
       fallback_share = fb, draws = if (return_draws) keep else NULL)
}
