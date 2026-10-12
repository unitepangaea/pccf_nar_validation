# =============================================================================
# ERR5  Postal-code-level correction: tables and validation
#
# 1. Address-level 50/50 split (within postal codes). Fold-1 addresses build the
#    postal-code table pi_p(t) (share of each postal code's residential addresses in
#    each true MSDI quintile); fold-2 addresses form the test population, so a cohort's
#    own addresses never inform its correction.
# 2. Same five cohort scenarios and synthetic outcome as err2 (true Q5:Q1 RR = 2), with
#    methods: oracle, naive, PC (ge_pc_correct: Gibbs data augmentation).
# 3. Final postal-code tables from all addresses for material and social quintiles:
#    data/derived/err/pc_quintile_probs_<dim>.csv.gz
#
# Run: Rscript scripts/err5_pc_validation.R  (env: ERR_REPS, ERR_SCEN as in err2)
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
source(file.path(ROOT, "R", "geoerror.R"))
OUT <- file.path(DERIVED, "err")
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
REPS <- as.integer(Sys.getenv("ERR_REPS", "100"))
SCEN <- strsplit(Sys.getenv("ERR_SCEN", "S1,S2,S3,S4,S5"), ",")[[1]]
P0 <- 0.05; RR <- 2

msdi <- as.data.table(read_excel(file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en",
                                           "1. EquivalenceTableCanada2021_en.xlsx"), sheet = "Data"))
idx <- list(material = msdi[, .(DAUID = as.character(DA), q = as.integer(QuintMat), score = SCOREMAT)][!is.na(q)],
            social   = msdi[, .(DAUID = as.character(DA), q = as.integer(QuintSoc), score = SCORESOC)][!is.na(q)])

# ---- Address-level data (built once) ---------------------------------------------------

f_addr <- file.path(OUT, "pc_addresses.rds")
if (!file.exists(f_addr)) {
  msg("Building address-level file")
  sac_group <- function(s) fcase(s == "1", "CMA", s %chin% c("2", "3"), "CA", s == "4", "Strong MIZ",
                                 s == "5", "Moderate MIZ", s == "6", "Weak MIZ", s == "7", "No MIZ",
                                 s == "8", "Territories", default = "Unknown")
  lf <- list.files(file.path(DERIVED, "linked"), "^linked_\\d+\\.rds$", full.names = TRUE)
  x <- rbindlist(lapply(lf, function(f) readRDS(f)[pccf_coded == TRUE,
         .(PCODE, prov = pccf_PR, sac = sac_group(pccf_SACtype),
           dmt3 = fcase(fifelse(pccf_DMT == "Z" & !is.na(pccf_DMTDIFF) & pccf_DMTDIFF != "", pccf_DMTDIFF, pccf_DMT) %chin% c("A", "B", "E", "G"), "Urban LC",
                        fifelse(pccf_DMT == "Z" & !is.na(pccf_DMTDIFF) & pccf_DMTDIFF != "", pccf_DMTDIFF, pccf_DMT) %chin% c("W", "H", "T"), "Rural",
                        default = "Other"),
           ref_source, ref_DAUID, pccf_DAUID)]))
  x <- x[!(ref_source == "blockface" & dmt3 == "Urban LC")][, ref_source := NULL]
  set.seed(20261009)
  x[, afold := sample(rep(1:2, length.out = .N))]
  saveRDS(x, f_addr)
}
x <- readRDS(f_addr)
msg(format(nrow(x), big.mark = ","), "addresses; fold sizes", paste(x[, .N, keyby = afold]$N, collapse = " / "))

# ---- Final postal-code tables (all addresses) ------------------------------------------

for (dim in names(idx)) {
  f <- file.path(OUT, sprintf("pc_quintile_probs_%s.csv.gz", dim))
  if (!file.exists(f)) {
    pt <- ge_pc_table(x[, .(n = .N), by = .(PCODE, ref_DAUID)], idx[[dim]])
    fwrite(pt, f)
    msg("Postal-code table", dim, ":", format(nrow(pt), big.mark = ","), "postal codes")
  }
}

# ---- Validation ------------------------------------------------------------------------

pt1 <- ge_pc_table(x[afold == 1, .(n = .N), by = .(PCODE, ref_DAUID)], idx$material)
te <- x[afold == 2]
te[idx$material, t := i.q, on = .(ref_DAUID = DAUID)]
te[idx$material, a := i.q, on = .(pccf_DAUID = DAUID)]
te <- te[!is.na(t) & !is.na(a)]
rm(x); invisible(gc())
miz <- c("Strong MIZ", "Moderate MIZ", "Weak MIZ", "No MIZ")
set.seed(20261010)

draw_cohort <- function(s) switch(s,
  S1 = te[runif(.N) < 100000 / .N],
  S2 = { r <- te[sac %in% miz]; r[runif(.N) < 30000 / .N] },
  S3 = { pcs <- sample(unique(te$PCODE), 2500); te[PCODE %in% pcs] },
  S4 = { p <- te$t / mean(te$t) * 60000 / nrow(te); te[runif(.N) < p] },
  S5 = { r <- te[prov == "47" & sac %in% miz]; r[runif(.N) < 10000 / .N] })

one_rep <- function(s, r) {
  co <- draw_cohort(s)
  co[, y := rbinom(.N, 1, P0 * RR^((t - 1) / 4))]
  tr_share <- co[, .N, keyby = t][, 100 * N / sum(N)]
  out <- list()
  for (m in c("oracle", "naive")) {
    g <- if (m == "oracle") co[, .(n = .N, e = sum(y)), keyby = .(q = t)] else co[, .(n = .N, e = sum(y)), keyby = .(q = a)]
    share <- 100 * g$n / sum(g$n)
    rr <- (g$e[5] / g$n[5]) / (g$e[1] / g$n[1]); se <- sqrt(1 / g$e[5] + 1 / g$e[1])
    out[[m]] <- data.table(method = m, q1 = share[1], q5 = share[5], q1_lo = NA, q1_hi = NA, q5_lo = NA, q5_hi = NA,
                           rr = rr, rr_lo = rr * exp(-1.96 * se), rr_hi = rr * exp(1.96 * se))
  }
  res <- ge_pc_correct(co$PCODE, co$a, pt1, y = co$y, iter = 800, burn = 300)
  d <- res$distribution; rrt <- res$rate_ratios
  out$PC <- data.table(method = "PC", q1 = d$corrected_pct[1], q5 = d$corrected_pct[5],
                       q1_lo = d$lo_pct[1], q1_hi = d$hi_pct[1], q5_lo = d$lo_pct[5], q5_hi = d$hi_pct[5],
                       rr = rrt$corrected_rr[5], rr_lo = rrt$lo_rr[5], rr_hi = rrt$hi_rr[5])
  cbind(scenario = s, rep = r, n = nrow(co), true_q1 = tr_share[1], true_q5 = tr_share[5],
        fallback = res$fallback_share, rbindlist(out, fill = TRUE))
}

for (s in SCEN) {
  t0 <- Sys.time()
  v <- rbindlist(lapply(seq_len(REPS), function(r) one_rep(s, r)))
  fwrite(v, file.path(OUT, paste0("pc_validation_replicates_", s, ".csv")))
  msg(s, "done in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)),
      "| fallback share", round(mean(v$fallback), 4))
}
msg("Done")
