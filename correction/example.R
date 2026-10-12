# =============================================================================
# Example: correcting MSDI material-quintile estimates for postal-code geocoding error
#
# Uses only the published tables in correction/tables and a SYNTHETIC cohort of
# 20,000 people with rural Ontario postal codes. In real use, the assigned quintile `a`
# and the strata (prov, sac, dmt3) come from your own PCCF+ output: the PCCF+ DA joined
# to the MSDI, the PCCF+ province (PR), SAC type (SACtype) and delivery mode (DMT),
# grouped as described in correction/README.md.
#
# Run from the repository root: Rscript correction/example.R   (about 1 minute)
# =============================================================================

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "R", "geoerror.R"))
TAB <- file.path(ROOT, "correction", "tables")
set.seed(2026)

cal <- readRDS(file.path(TAB, "calibration_material_B.rds"))           # province x SAC x DMT3
pct <- as.data.table(read.csv(gzfile(file.path(TAB, "pc_quintile_probs_material.csv.gz")),
                             colClasses = c(PCODE = "character")))   # postal code -> true quintile shares

# ---- Synthetic cohort ------------------------------------------------------------------
# Rural Ontario postal codes (second character 0), sampled in proportion to addresses.
# True quintile t is drawn from each postal code's distribution; the "PCCF+" quintile a
# is simulated from the published classification probabilities for one stratum. Because
# these synthetic assigned quintiles ignore the postal code, the example shows how to call
# the functions, not their relative accuracy (see Table 5 of the paper for that).
pc <- pct[grepl("^[KLMNP]0", PCODE)]
co <- pc[data.table(PCODE = sample(pc$PCODE, 20000, replace = TRUE, prob = pc$n_addr)), on = "PCODE"]
pm <- as.matrix(co[, .(p1, p2, p3, p4, p5)])
co <- co[, .(PCODE, t = apply(pm, 1, function(p) sample(1:5, 1, prob = p)))]
co[, `:=`(prov = "35", sac = "Moderate MIZ", dmt3 = "Rural")]
P <- cal$class[prov == "35" & sac == "Moderate MIZ" & dmt3 == "Rural"]
co[, a := vapply(t, function(tt) sample(1:5, 1, prob = P[t == tt, p_a_given_t]), 1L)]
co[, y := rbinom(.N, 1, 0.05 * 2^((t - 1) / 4))]                       # true Q5:Q1 rate ratio 2.0

# ---- 1. Stratum correction (aggregate counts by stratum and assigned quintile) ---------
obs <- co[, .(n = .N, events = sum(y)), by = .(prov, sac, dmt3, a)]
st <- ge_correct(obs, cal, draws = 1000, sample_n = TRUE, seed = 1)
print(st$rate_ratios[, .(q, naive_rr, corrected_rr, lo_rr, hi_rr)])   # 95% simulation interval
print(st$diagnostics)                                                    # EM convergence, boundary share

# ---- 2. Postal-code correction (record-level postal codes) -----------------------------
pcr <- ge_pc_correct(co$PCODE, co$a, pct, y = co$y, iter = 1500, burn = 500, seed = 1)
print(pcr$rate_ratios[, .(q, naive_rr, corrected_rr, lo_rr, hi_rr)])  # 95% credible interval
cat("Share of records whose postal code is not in the table:", pcr$fallback_share, "\n")

# ---- Truth in the synthetic cohort ------------------------------------------------------
tr <- co[, .(rate = mean(y)), keyby = t]
cat("True rate ratios:", round(tr$rate / tr$rate[1], 2), "\n")
