# =============================================================================
# ERR8  Convergence diagnostics for the postal-code correction (Gibbs data augmentation)
#
# For cohorts drawn as in the validation (err5 scenarios S1, S2, S5; err7 cohort K2 rural
# delivery in CMA/CA), runs 4 chains with dispersed starting rates and the validation
# settings (800 iterations, first 300 discarded) and reports, for the log Q5:Q1 rate ratio,
# the rank-normalized split R-hat and bulk effective sample size (Vehtari et al. 2021).
# Output: data/derived/err/pc_mcmc_diagnostics.csv, pc_mcmc_chains.rds (retained draws)
# Run: Rscript scripts/err8_pc_diagnostics.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
source(file.path(ROOT, "R", "geoerror.R"))
OUT <- file.path(DERIVED, "err")
set.seed(20261009)
COHORTS <- 10; CHAINS <- 4; ITER <- 800; BURN <- 300; RR <- 2

# rank-normalized split R-hat and bulk ESS; x is a matrix [draws, chains]
split_chains <- function(x) { h <- floor(nrow(x) / 2); cbind(x[1:h, , drop = FALSE], x[(h + 1):(2 * h), , drop = FALSE]) }
rank_norm <- function(x) { r <- rank(x, ties.method = "average"); matrix(qnorm((r - 3 / 8) / (length(x) + 1 / 4)), nrow(x)) }
rhat_basic <- function(x) {
  n <- nrow(x); W <- mean(apply(x, 2, var)); B <- n * var(colMeans(x))
  sqrt(((n - 1) / n * W + B / n) / W)
}
rhat <- function(x) {
  s <- split_chains(x); f <- abs(s - median(s))
  max(rhat_basic(rank_norm(s)), rhat_basic(rank_norm(f)))
}
ess_bulk <- function(x) {
  s <- rank_norm(split_chains(x)); n <- nrow(s); m <- ncol(s)
  acor <- apply(s, 2, function(v) { v <- v - mean(v); f <- fft(c(v, rep(0, n)))   # autocorrelation by FFT
    ac <- Re(fft(f * Conj(f), inverse = TRUE))[1:n] / (2 * n); ac / ac[1] })
  cv <- apply(s, 2, var); W <- mean(cv); vplus <- (n - 1) / n * W + var(colMeans(s))
  rho <- 1 - (W - rowMeans(sweep(acor, 2, cv, "*"))) / vplus; rho[1] <- 1
  tau <- -1; k <- 1                          # Geyer initial positive sequence
  while (k + 1 <= n) { p <- rho[k] + rho[k + 1]; if (p < 0) break; tau <- tau + 2 * p; k <- k + 2 }
  m * n / max(tau, 1 / log10(m * n))
}

msdi <- as.data.table(read_excel(file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en",
                                           "1. EquivalenceTableCanada2021_en.xlsx"), sheet = "Data"))
ix <- msdi[, .(DAUID = as.character(DA), q = as.integer(QuintMat))][!is.na(q)]
x <- readRDS(file.path(OUT, "pc_addresses.rds"))
pt1 <- ge_pc_table(x[afold == 1, .(n = .N), by = .(PCODE, ref_DAUID)], ix)
te <- x[afold == 2]; rm(x); invisible(gc())
te[ix, t := i.q, on = .(ref_DAUID = DAUID)]; te[ix, a := i.q, on = .(pccf_DAUID = DAUID)]
te <- te[!is.na(t) & !is.na(a)]
miz <- c("Strong MIZ", "Moderate MIZ", "Weak MIZ", "No MIZ")
coh <- list(
  S1_national        = function() te[runif(.N) < 100000 / .N],
  S2_rural           = function() { r <- te[sac %in% miz]; r[runif(.N) < 30000 / .N] },
  S5_SK_rural        = function() { r <- te[prov == "47" & sac %in% miz]; r[runif(.N) < 10000 / .N] },
  K2_rural_in_CMA_CA = function() { r <- te[sac %in% c("CMA", "CA") & dmt3 == "Rural"]; r[runif(.N) < 30000 / .N] })

chains <- list()
res <- rbindlist(lapply(names(coh), function(k) rbindlist(lapply(seq_len(COHORTS), function(r) {
  co <- coh[[k]](); co[, y := rbinom(.N, 1, 0.05 * RR^((t - 1) / 4))]
  ch <- sapply(seq_len(CHAINS), function(c) {
    d <- ge_pc_correct(co$PCODE, co$a, pt1, y = co$y, iter = ITER, burn = BURN, return_draws = TRUE,
                       lam_init = mean(co$y) * exp(runif(5, log(0.5), log(2))))$draws
    log(d[, 5] / d[, 1])
  })
  chains[[paste(k, r)]] <<- ch
  data.table(cohort = k, rep = r, n = nrow(co), rhat = rhat(ch), ess_bulk = ess_bulk(ch),
             ess_per_chain = ess_bulk(ch) / CHAINS)
}))))
saveRDS(chains, file.path(OUT, "pc_mcmc_chains.rds"))   # log Q5:Q1 draws, [iteration, chain] per cohort
fwrite(res, file.path(OUT, "pc_mcmc_diagnostics.csv"))
print(res[, .(cohorts = .N, rhat_max = round(max(rhat), 3), rhat_median = round(median(rhat), 3),
              ess_min = round(min(ess_bulk)), ess_median = round(median(ess_bulk))), by = cohort])
