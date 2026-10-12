# =============================================================================
# ERR2  Validation of the postal-code error correction (geoerror.R)
#
# Calibration is built from fold 1 (half of postal codes) and applied to simulated
# cohorts drawn from fold 2, whose true and PCCF+-assigned MSDI quintiles are known.
# Each cohort carries a synthetic outcome with a known gradient:
#   P(Y = 1 | true quintile t) = P0 * RR^((t - 1) / 4), so the Q5:Q1 rate ratio is RR = 2.
#
# Scenarios (REPS replicates each)
#   S1 national   simple random sample, n ~ 100,000
#   S2 rural      simple random sample from MIZ strata (assigned SAC), n ~ 30,000
#   S3 clustered  all addresses in 2,500 random postal codes (n ~ 20,000-40,000)
#   S4 skewed     inclusion probability proportional to true quintile (deprived cohort), n ~ 60,000
#   S5 SK rural   Saskatchewan MIZ strata, n ~ 10,000
# Methods
#   naive    estimates by assigned quintile (Poisson CI for the rate ratio)
#   A        correction stratified by province x SAC
#   B        correction stratified by province x SAC x delivery mode (PCCF+ DMT)
#   C        correction stratified by province x delivery mode
# Calibrations include the design effect of clustering within postal codes (pc_transitions).
# Monte Carlo draws: ERR_DRAWS (default 200).
#   oracle   estimates by true quintile (Poisson CI), the best achievable
# Metrics: bias and 95% CI coverage for the share in Q1 and Q5 (target: the cohort's true
# share) and for the Q5:Q1 rate ratio (target: RR = 2).
#
# Outputs: data/derived/err/validation_replicates.csv, data/derived/err/validation_summary.csv,
#          data/derived/err/calibration_fold1_*.rds
# Run: Rscript scripts/err2_validation.R   (env vars: ERR_REPS, ERR_DRAWS; ERR_SCEN="S1"
#      runs one scenario into validation_replicates_S1.csv so scenarios can run in parallel;
#      ERR_SUMMARY=1 only summarizes existing replicate files)
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
source(file.path(ROOT, "R", "geoerror.R"))
OUT <- file.path(DERIVED, "err")
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

REPS  <- as.integer(Sys.getenv("ERR_REPS", "100"))
SCEN  <- strsplit(Sys.getenv("ERR_SCEN", "S1,S2,S3,S4,S5"), ",")[[1]]
SUMMARY_ONLY <- nzchar(Sys.getenv("ERR_SUMMARY"))
METHODS <- strsplit(Sys.getenv("ERR_METHODS", "A,B,C"), ",")[[1]]   # A prov x SAC; B + DMT3; C prov x DMT3
TAG <- if (identical(METHODS, c("A", "B", "C"))) "" else paste0(paste(METHODS, collapse = ""), "_")
DRAWS <- as.integer(Sys.getenv("ERR_DRAWS", "200"))
P0 <- 0.05; RR <- 2
set.seed(20261008)

summarize_validation <- function() {
  v <- rbindlist(lapply(list.files(OUT, "^validation_replicates_S\\d\\.csv$", full.names = TRUE), fread))
  summ <- v[, .(reps = .N, n_mean = round(mean(n)),
                q5_true = mean(true_q5), q5_est = mean(q5), q5_bias = mean(q5 - true_q5),
                q5_cover = mean(q5_lo <= true_q5 & true_q5 <= q5_hi),
                q1_true = mean(true_q1), q1_est = mean(q1), q1_bias = mean(q1 - true_q1),
                q1_cover = mean(q1_lo <= true_q1 & true_q1 <= q1_hi),
                rr_mean = exp(mean(log(rr))), rr_bias_pct = 100 * (exp(mean(log(rr))) / RR - 1),
                rr_rmse_log = sqrt(mean((log(rr) - log(RR))^2)),
                rr_cover = mean(rr_lo <= RR & RR <= rr_hi),
                rr_ci_width = exp(mean(log(rr_hi / rr_lo))),
                em_converged = mean(em_converged), boundary_n_pop = mean(boundary_n_pop),
                boundary_n = mean(boundary_n), boundary_events = mean(boundary_events)), by = .(scenario, method)]
  num <- names(summ)[vapply(summ, is.numeric, TRUE)]
  summ[, (num) := lapply(.SD, round, 3), .SDcols = num]
  fwrite(summ, file.path(OUT, "validation_summary.csv"))
  print(summ)
}
if (SUMMARY_ONLY) { summarize_validation(); quit(save = "no") }

# ---- Index and calibration (fold 1) -------------------------------------------------

msdi <- as.data.table(read_excel(file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en",
                                           "1. EquivalenceTableCanada2021_en.xlsx"), sheet = "Data"))
idx_mat <- msdi[, .(DAUID = as.character(DA), q = as.integer(QuintMat), score = SCOREMAT)][!is.na(q)]

tr <- readRDS(file.path(OUT, "da_transitions.rds"))
ptr <- readRDS(file.path(OUT, "pc_transitions.rds"))[fold == 1]     # design effect (clustering)
strata_def <- list(A = c("prov", "sac"), B = c("prov", "sac", "dmt3"), C = c("prov", "dmt3"))
cal <- lapply(strata_def[METHODS], function(st) ge_build_calibration(tr[fold == 1], idx_mat, st, clusters = ptr))
for (m in names(cal)) saveRDS(cal[[m]], file.path(OUT, paste0("calibration_fold1_", m, ".rds")))
msg("Calibration strata:", paste(names(cal), vapply(cal, function(c) nrow(c$summary), 1L), collapse = " | "))

# ---- Test population (fold 2) ---------------------------------------------------------

te <- readRDS(file.path(OUT, "test_addresses.rds"))
te[idx_mat, t := i.q, on = .(ref_DAUID = DAUID)]
te[idx_mat, a := i.q, on = .(pccf_DAUID = DAUID)]
te <- te[!is.na(t) & !is.na(a)]
miz <- c("Strong MIZ", "Moderate MIZ", "Weak MIZ", "No MIZ")
msg(format(nrow(te), big.mark = ","), "fold-2 addresses with both quintiles")

draw_cohort <- function(s) {
  switch(s,
    S1 = te[runif(.N) < 100000 / .N],
    S2 = { r <- te[sac %in% miz]; r[runif(.N) < 30000 / .N] },
    S3 = { pcs <- sample(unique(te$PCODE), 2500); te[PCODE %in% pcs] },
    S4 = { p <- te$t / mean(te$t) * 60000 / nrow(te); te[runif(.N) < p] },
    S5 = { r <- te[prov == "47" & sac %in% miz]; r[runif(.N) < 10000 / .N] })
}

# ---- Simulation -----------------------------------------------------------------------

one_rep <- function(s, r) {
  co <- draw_cohort(s)
  co[, y := rbinom(.N, 1, P0 * RR^((t - 1) / 4))]
  tr_share <- co[, .N, keyby = t][, 100 * N / sum(N)]
  out <- list()
  # oracle and naive
  for (m in c("oracle", "naive")) {
    g <- if (m == "oracle") co[, .(n = .N, e = sum(y)), keyby = .(q = t)] else co[, .(n = .N, e = sum(y)), keyby = .(q = a)]
    share <- 100 * g$n / sum(g$n)
    rr <- (g$e[5] / g$n[5]) / (g$e[1] / g$n[1]); se <- sqrt(1 / g$e[5] + 1 / g$e[1])
    out[[m]] <- data.table(method = m, q1 = share[1], q5 = share[5], q1_lo = NA, q1_hi = NA, q5_lo = NA, q5_hi = NA,
                           rr = rr, rr_lo = rr * exp(-1.96 * se), rr_hi = rr * exp(1.96 * se),
                           em_converged = NA_real_, boundary_n = NA_real_, boundary_n_pop = NA_real_,
                           boundary_events = NA_real_)
  }
  for (m in METHODS) {
    st <- cal[[m]]$strata
    ob <- co[, .(n = .N, events = sum(y)), by = c(st, "a")]
    res <- ge_correct(ob, cal[[m]], draws = DRAWS, sample_n = TRUE)
    d <- res$distribution; rrt <- res$rate_ratios
    out[[m]] <- data.table(method = m, q1 = d$corrected_pct[1], q5 = d$corrected_pct[5],
                           q1_lo = d$lo_pct[1], q1_hi = d$hi_pct[1], q5_lo = d$lo_pct[5], q5_hi = d$hi_pct[5],
                           rr = rrt$corrected_rr[5], rr_lo = rrt$lo_rr[5], rr_hi = rrt$hi_rr[5],
                           res$diagnostics)
  }
  cbind(scenario = s, rep = r, n = nrow(co), true_q1 = tr_share[1], true_q5 = tr_share[5], rbindlist(out))
}

for (s in SCEN) {
  t0 <- Sys.time()
  v <- rbindlist(lapply(seq_len(REPS), function(r) one_rep(s, r)))
  fwrite(v, file.path(OUT, paste0("validation_replicates_", TAG, s, ".csv")))
  msg(s, "done in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
}
if (length(SCEN) == 5) summarize_validation()
msg("Done")
