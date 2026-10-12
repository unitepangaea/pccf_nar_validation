# =============================================================================
# ERR7  Postal-code correction for cohorts with atypical delivery-mode mixes
#
# Complements err6 (stratum methods): the same cohort definitions are drawn from the
# address-split validation half (pc_addresses.rds, afold 2) and corrected with the
# postal-code method using fold-1 postal-code tables. True Q5:Q1 rate ratio 2.0.
#   K2 town dwellers in MIZ: urban letter-carrier codes in MIZ strata (n ~ 30,000)
#   K3 rural codes in CMA/CA: rural delivery codes in CMA/CA strata (n ~ 30,000)
#   K4 MIZ, typical mix (n ~ 30,000)
# Output: data/derived/err/dmt_eval/pc_atypical_cohorts.csv
# Run: Rscript scripts/err7_pc_atypical.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
source(file.path(ROOT, "R", "geoerror.R"))
OUT <- file.path(DERIVED, "err")
set.seed(20261009)
REPS <- 40; RR <- 2

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
  K2_town_in_MIZ     = function() { r <- te[sac %in% miz & dmt3 == "Urban LC"]; r[runif(.N) < 30000 / .N] },
  K3_rural_in_CMA_CA = function() { r <- te[sac %in% c("CMA", "CA") & dmt3 == "Rural"]; r[runif(.N) < 30000 / .N] },
  K4_MIZ_typical     = function() { r <- te[sac %in% miz]; r[runif(.N) < 30000 / .N] })

res <- rbindlist(lapply(names(coh), function(k) rbindlist(lapply(seq_len(REPS), function(r) {
  co <- coh[[k]](); co[, y := rbinom(.N, 1, 0.05 * RR^((t - 1) / 4))]
  out <- list()
  for (m in c("oracle", "naive")) {
    g <- if (m == "oracle") co[, .(n = .N, e = sum(y)), keyby = .(q = t)] else co[, .(n = .N, e = sum(y)), keyby = .(q = a)]
    rr <- (g$e[5] / g$n[5]) / (g$e[1] / g$n[1]); se <- sqrt(1 / g$e[5] + 1 / g$e[1])
    out[[m]] <- data.table(method = m, rr = rr, lo = rr * exp(-1.96 * se), hi = rr * exp(1.96 * se))
  }
  p <- ge_pc_correct(co$PCODE, co$a, pt1, y = co$y, iter = 800, burn = 300)$rate_ratios
  out$PC <- data.table(method = "postal code", rr = p$corrected_rr[5], lo = p$lo_rr[5], hi = p$hi_rr[5])
  cbind(cohort = k, rep = r, n = nrow(co), rbindlist(out))
}))))
summ <- res[, .(n = round(mean(n)), rr_mean = round(exp(mean(log(rr))), 3),
                bias_pct = round(100 * (exp(mean(log(rr))) / RR - 1), 1),
                coverage = mean(lo <= RR & RR <= hi), ci_width = round(exp(mean(log(hi / lo))), 2)),
            by = .(cohort, method)]
fwrite(summ, file.path(OUT, "dmt_eval", "pc_atypical_cohorts.csv")); print(summ)
