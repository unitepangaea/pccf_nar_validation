# =============================================================================
# ERR6  Should geocoding error measures be stratified by delivery mode (DMT)?
#
# Evaluation only (does not change the paper or published tables). Outputs in
# data/derived/err/dmt_eval/.
#   1. Out-of-sample explanatory power of candidate stratifications for three error
#      measures (material-quintile misclassification, log10 positional error, same DA):
#      cell rates/means estimated on fold 1 (half of postal codes), scored on fold 2.
#   2. Heterogeneity: misclassification and positional error by DMT x province and
#      DMT x SAC, with the spread across provinces within each DMT.
#   3. Transportability: predicted vs actual misclassification for cohorts with typical
#      and atypical delivery-mode mixes, using fold-1 rates from each stratification.
#   4. Correction: stratum correction of the Q5:Q1 rate ratio (true 2.0) for atypical
#      cohorts with province x SAC, province x SAC x DMT and province x DMT strata.
# Run: Rscript scripts/err6_dmt_eval.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
source(file.path(ROOT, "R", "geoerror.R"))
OUT <- file.path(DERIVED, "err", "dmt_eval")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
set.seed(20261008)

sac_group <- function(s) fcase(s == "1", "CMA", s %chin% c("2", "3"), "CA", s == "4", "Strong MIZ",
                               s == "5", "Moderate MIZ", s == "6", "Weak MIZ", s == "7", "No MIZ",
                               s == "8", "Territories", default = "Unknown")
dmt_eff <- function(d, dd) fifelse(d == "Z" & !is.na(dd) & dd != "", dd, d)
dmt7_of <- function(d) fcase(d == "A", "A house", d == "B", "B apartment", d %chin% c("E", "G"), "E/G business",
                             d == "W", "W rural", d %chin% c("H", "T"), "H/T rural route",
                             d %chin% c("J", "K", "M", "R", "X"), "J/K/M PO box/GD", default = "Z/unknown")
dmt3_of <- function(d) fcase(d %chin% c("A", "B", "E", "G"), "Urban LC", d %chin% c("W", "H", "T"), "Rural",
                             default = "Other")

# ---- Address-level frame ------------------------------------------------------------

f_frame <- file.path(OUT, "frame.rds")
if (!file.exists(f_frame)) {
  msg("Building address frame")
  lf <- list.files(file.path(DERIVED, "linked"), "^linked_\\d+\\.rds$", full.names = TRUE)
  x <- rbindlist(lapply(lf, function(f) readRDS(f)[pccf_coded == TRUE,
         .(ADDR_GUID, PCODE, prov = pccf_PR, sac = sac_group(pccf_SACtype),
           dmt = dmt_eff(pccf_DMT, pccf_DMTDIFF), ref_source, ref_DAUID, pccf_DAUID)]))
  x[, `:=`(dmt7 = dmt7_of(dmt), dmt3 = dmt3_of(dmt))]
  x <- x[!(ref_source == "blockface" & dmt3 == "Urban LC")][, ref_source := NULL]
  d2 <- readRDS(file.path(DERIVED, "p2", "p2_analysis.rds"))[pccf_coded == TRUE, .(ADDR_GUID, dist_m)]
  x <- d2[x, on = "ADDR_GUID"]; rm(d2)
  msdi <- as.data.table(read_excel(file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en",
                                             "1. EquivalenceTableCanada2021_en.xlsx"), sheet = "Data"))
  ix <- msdi[, .(DAUID = as.character(DA), q = as.integer(QuintMat))][!is.na(q)]
  x[ix, t := i.q, on = .(ref_DAUID = DAUID)]
  x[ix, a := i.q, on = .(pccf_DAUID = DAUID)]
  x[, `:=`(miscl = as.integer(t != a), ldist = log10(dist_m + 1), same_da = as.integer(ref_DAUID == pccf_DAUID))]
  pc <- unique(x$PCODE)
  fold_of <- setNames(sample(rep(1:2, length.out = length(pc))), pc)
  x[, fold := fold_of[PCODE]]
  x[, c("ADDR_GUID", "dist_m") := NULL]
  saveRDS(x, f_frame)
}
# ERR6_ONLY_CORRECTION=1 reruns only section 4 (after changes to the correction method)
ONLY4 <- nzchar(Sys.getenv("ERR6_ONLY_CORRECTION"))
if (!ONLY4) {
x <- readRDS(f_frame)
msg(format(nrow(x), big.mark = ","), "addresses")

# ---- 1. Out-of-sample explanatory power ------------------------------------------------

groupings <- list(
  "None (national)" = character(0), "Province" = "prov", "SAC" = "sac", "DMT" = "dmt7",
  "Province x SAC" = c("prov", "sac"), "Province x DMT" = c("prov", "dmt7"), "SAC x DMT" = c("sac", "dmt7"),
  "Province x SAC x DMT3" = c("prov", "sac", "dmt3"), "Province x SAC x DMT" = c("prov", "sac", "dmt7"))

score <- function(g) {
  tr <- x[fold == 1]; te <- x[fold == 2]
  # misclassification (binary, MSDI available) --------------------------------
  tb <- tr[!is.na(miscl)]; eb <- te[!is.na(miscl)]
  p0 <- mean(tb$miscl)
  cb <- if (length(g)) tb[, .(n = .N, k = sum(miscl)), by = g] else data.table(n = nrow(tb), k = sum(tb$miscl))
  cb[, p := (k + 20 * p0) / (n + 20)]                       # light shrinkage toward national
  pe <- if (length(g)) cb[eb, on = g, p] else rep(cb$p, nrow(eb))
  pe[is.na(pe)] <- p0
  ll  <- sum(eb$miscl * log(pe) + (1 - eb$miscl) * log(1 - pe))
  ll0 <- sum(eb$miscl * log(p0) + (1 - eb$miscl) * log(1 - p0))
  brier <- mean((eb$miscl - pe)^2); brier0 <- mean((eb$miscl - p0)^2)
  # log distance -----------------------------------------------------------------
  m0 <- mean(tr$ldist)
  cd <- if (length(g)) tr[, .(n = .N, m = mean(ldist)), by = g] else data.table(n = nrow(tr), m = m0)
  md <- if (length(g)) cd[te, on = g, m] else rep(m0, nrow(te)); md[is.na(md)] <- m0
  r2d <- 1 - mean((te$ldist - md)^2) / mean((te$ldist - m0)^2)
  # same DA ------------------------------------------------------------------------
  s0 <- mean(tr$same_da)
  cs <- if (length(g)) tr[, .(n = .N, s = mean(same_da)), by = g] else data.table(n = nrow(tr), s = s0)
  ms <- if (length(g)) cs[te, on = g, s] else rep(s0, nrow(te)); ms[is.na(ms)] <- s0
  r2s <- 1 - mean((te$same_da - ms)^2) / mean((te$same_da - s0)^2)
  data.table(cells = nrow(cb), cells_lt100 = sum(cb$n < 100),
             misclass_pseudoR2 = 1 - ll / ll0, misclass_brier_skill = 1 - brier / brier0,
             logdist_R2 = r2d, sameDA_R2 = r2s)
}
msg("1. Explanatory power")
ex <- rbindlist(lapply(names(groupings), function(n) cbind(stratification = n, score(groupings[[n]]))))
num <- names(ex)[vapply(ex, is.double, TRUE)]; ex[, (num) := lapply(.SD, round, 4), .SDcols = num]
fwrite(ex, file.path(OUT, "explanatory_power.csv")); print(ex)

# ---- 2. Heterogeneity --------------------------------------------------------------------

msg("2. Heterogeneity")
ab <- c("10"="NL","11"="PE","12"="NS","13"="NB","24"="QC","35"="ON","46"="MB","47"="SK","48"="AB","59"="BC","60"="YT","61"="NT","62"="NU")
hp <- x[!is.na(miscl), .(N = .N, pct_misclassified = round(100 * mean(miscl), 1),
                         median_m = round(10^median(ldist) - 1), pct_same_da = round(100 * mean(same_da), 1)),
        by = .(dmt7, prov)][, province := ab[prov]]
fwrite(hp, file.path(OUT, "dmt_by_province.csv"))
spread <- hp[N >= 1000, .(provinces = .N, min_pct = min(pct_misclassified), max_pct = max(pct_misclassified),
                          range_pts = max(pct_misclassified) - min(pct_misclassified),
                          min_median_m = min(median_m), max_median_m = max(median_m)), by = dmt7]
nat <- x[!is.na(miscl), .(N = .N, national_pct = round(100 * mean(miscl), 1)), by = dmt7]
spread <- nat[spread, on = "dmt7"][order(-N)]
fwrite(spread, file.path(OUT, "dmt_spread_across_provinces.csv")); print(spread)
hs <- x[!is.na(miscl), .(N = .N, pct_misclassified = round(100 * mean(miscl), 1)), by = .(dmt7, sac)]
fwrite(hs, file.path(OUT, "dmt_by_sac.csv"))
print(dcast(hs[N >= 1000], dmt7 ~ sac, value.var = "pct_misclassified"))
mix <- x[, .N, by = .(sac, dmt3)][, pct := round(100 * N / sum(N), 1), by = sac]
fwrite(mix, file.path(OUT, "dmt_mix_within_sac.csv"))
print(dcast(mix, sac ~ dmt3, value.var = "pct"))

# ---- 3. Transportability of stratum error rates -------------------------------------------

msg("3. Transportability")
te <- x[fold == 2 & !is.na(miscl)]
miz <- c("Strong MIZ", "Moderate MIZ", "Weak MIZ", "No MIZ")
cohorts <- list(
  K1_national        = function() te[runif(.N) < 100000 / .N],
  K2_town_in_MIZ     = function() { r <- te[sac %in% miz & dmt3 == "Urban LC"]; r[runif(.N) < 30000 / .N] },
  K3_rural_in_CMA_CA = function() { r <- te[sac %in% c("CMA", "CA") & dmt3 == "Rural"]; r[runif(.N) < 30000 / .N] },
  K4_MIZ_typical     = function() { r <- te[sac %in% miz]; r[runif(.N) < 30000 / .N] },
  K5_ON_rural_codes  = function() { r <- te[prov == "35" & dmt7 == "W rural"]; r[runif(.N) < 20000 / .N] })
pred_sets <- groupings[c("Province x SAC", "Province x DMT", "Province x SAC x DMT3", "Province x SAC x DMT")]
tr1 <- x[fold == 1 & !is.na(miscl)]
rates <- lapply(pred_sets, function(g) tr1[, .(p = mean(miscl)), by = g])
tp <- rbindlist(lapply(names(cohorts), function(k) rbindlist(lapply(1:50, function(r) {
  co <- cohorts[[k]]()
  rbindlist(lapply(names(pred_sets), function(g) {
    pr <- rates[[g]][co, on = pred_sets[[g]], p]; pr[is.na(pr)] <- mean(tr1$miscl)
    data.table(cohort = k, rep = r, stratification = g, actual = 100 * mean(co$miscl), predicted = 100 * mean(pr))
  }))
}))))
tsum <- tp[, .(actual = round(mean(actual), 1), predicted = round(mean(predicted), 1),
               mean_abs_error_pts = round(mean(abs(predicted - actual)), 2)), by = .(cohort, stratification)]
fwrite(tsum, file.path(OUT, "transportability.csv")); print(tsum)
}

# ---- 4. Correction of rate ratios for atypical cohorts ---------------------------------------

msg("4. Correction")
tr_tab <- readRDS(file.path(DERIVED, "err", "da_transitions.rds"))
msdi <- as.data.table(read_excel(file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en",
                                           "1. EquivalenceTableCanada2021_en.xlsx"), sheet = "Data"))
ix <- msdi[, .(DAUID = as.character(DA), q = as.integer(QuintMat))][!is.na(q)]
# calibration from the err1 fold 1 (err1 used its own postal-code split; this frame's folds differ,
# so cohorts below are drawn from err1's fold-2 test file to keep calibration and test disjoint)
ptr <- readRDS(file.path(DERIVED, "err", "pc_transitions.rds"))[fold == 1]   # design effect
cal <- list(A = ge_build_calibration(tr_tab[fold == 1], ix, c("prov", "sac"), clusters = ptr),
            B = ge_build_calibration(tr_tab[fold == 1], ix, c("prov", "sac", "dmt3"), clusters = ptr),
            C = ge_build_calibration(tr_tab[fold == 1], ix, c("prov", "dmt3"), clusters = ptr))
miz <- c("Strong MIZ", "Moderate MIZ", "Weak MIZ", "No MIZ")
tst <- readRDS(file.path(DERIVED, "err", "test_addresses.rds"))
tst[ix, t := i.q, on = .(ref_DAUID = DAUID)]; tst[ix, a := i.q, on = .(pccf_DAUID = DAUID)]
tst <- tst[!is.na(t) & !is.na(a)]
ccoh <- list(
  K2_town_in_MIZ     = function() { r <- tst[sac %in% miz & dmt3 == "Urban LC"]; r[runif(.N) < 30000 / .N] },
  K3_rural_in_CMA_CA = function() { r <- tst[sac %in% c("CMA", "CA") & dmt3 == "Rural"]; r[runif(.N) < 30000 / .N] },
  K4_MIZ_typical     = function() { r <- tst[sac %in% miz]; r[runif(.N) < 30000 / .N] })
RR <- 2
cr <- rbindlist(lapply(names(ccoh), function(k) rbindlist(lapply(1:40, function(r) {
  co <- ccoh[[k]](); co[, y := rbinom(.N, 1, 0.05 * RR^((t - 1) / 4))]
  g <- co[, .(n = .N, e = sum(y)), keyby = .(q = a)]
  nrr <- (g$e[5] / g$n[5]) / (g$e[1] / g$n[1]); se <- sqrt(1 / g$e[5] + 1 / g$e[1])
  out <- list(data.table(method = "naive", rr = nrr, lo = nrr * exp(-1.96 * se), hi = nrr * exp(1.96 * se)))
  for (m in names(cal)) {
    ob <- co[, .(n = .N, events = sum(y)), by = c(cal[[m]]$strata, "a")]
    cc <- ge_correct(ob, cal[[m]], draws = 200, sample_n = TRUE); res <- cc$rate_ratios
    out[[m]] <- data.table(method = m, rr = res$corrected_rr[5], lo = res$lo_rr[5], hi = res$hi_rr[5],
                           cc$diagnostics)
  }
  cbind(cohort = k, rep = r, rbindlist(out, fill = TRUE))
}))))
csum <- cr[, .(rr_mean = round(exp(mean(log(rr))), 3), bias_pct = round(100 * (exp(mean(log(rr))) / RR - 1), 1),
               coverage = mean(lo <= RR & RR <= hi), ci_width = round(exp(mean(log(hi / lo))), 2),
               em_converged = mean(em_converged), boundary_n_pop = mean(boundary_n_pop)),
           by = .(cohort, method)]
csum[, method := fcase(method == "A", "Province x SAC", method == "B", "Province x SAC x DMT3",
                       method == "C", "Province x DMT3", default = method)]
fwrite(csum, file.path(OUT, "correction_atypical_cohorts.csv")); print(csum)
msg("Done")
