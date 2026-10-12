# =============================================================================
# ERR1  Geocoding-error transition tables for postal-code geocoding (PCCF+ 8D)
#
# Builds the inputs for the error-correction method (R/geoerror.R):
#   da_transitions.rds  counts of addresses by stratum x reference DA x PCCF+ DA, with
#                       a 50/50 split by postal code (fold 1 = calibration, fold 2 =
#                       validation). Any DA-level index (MSDI, CAN-Marg, CIMD...) can be
#                       merged onto the two DAUID columns to build its own calibration.
#   stratum_positional_error.csv  positional-error and DA-agreement summaries by stratum
#   pc_transitions.rds  the same counts by postal code (for the clustering design effect)
#   test_addresses.rds  address-level fold-2 records used by the validation simulations
#
# Strata are defined only from information a data user has after running PCCF+:
#   prov  province of the PCCF+ assignment (pccf_PR)
#   sac   SAC type of the PCCF+ assignment, grouped: CMA (1), CA (2-3), Strong MIZ (4),
#         Moderate MIZ (5), Weak MIZ (6), No MIZ (7), Territories (8), Unknown (other)
#   dmt3  delivery mode from the PCCF+ DMT: urban letter carrier (A, B, E, G), rural
#         (W, H, T), other (J, K, M, R, X, retired/unknown)
#
# Sample: P3 main sample (PCCF+ coded; urban letter-carrier addresses with a blockface
# reference excluded because that comparison is partly circular).
# Run: Rscript scripts/err1_transitions.R
# =============================================================================

suppressPackageStartupMessages(library(data.table))

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
OUT  <- file.path(DERIVED, "err")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
SEED <- 20261007

sac_group <- function(s) fcase(s == "1", "CMA", s %chin% c("2", "3"), "CA", s == "4", "Strong MIZ",
                               s == "5", "Moderate MIZ", s == "6", "Weak MIZ", s == "7", "No MIZ",
                               s == "8", "Territories", default = "Unknown")
dmt3 <- function(d, dd) {
  d <- fifelse(d == "Z" & !is.na(dd) & dd != "", dd, d)
  fcase(d %chin% c("A", "B", "E", "G"), "Urban LC", d %chin% c("W", "H", "T"), "Rural", default = "Other")
}

# ---- Assemble address-level records ------------------------------------------------

msg("Reading P1 linked files")
lf <- list.files(file.path(DERIVED, "linked"), "^linked_\\d+\\.rds$", full.names = TRUE)
x <- rbindlist(lapply(lf, function(f) readRDS(f)[pccf_coded == TRUE,
  .(ADDR_GUID, PCODE, prov = pccf_PR, sac_code = pccf_SACtype, pccf_DMT, pccf_DMTDIFF,
    ref_source, ref_DAUID, pccf_DAUID)]))
msg("Reading P2 distances")
d2 <- readRDS(file.path(DERIVED, "p2", "p2_analysis.rds"))[pccf_coded == TRUE, .(ADDR_GUID, dist_m)]
x <- d2[x, on = "ADDR_GUID"]; rm(d2); invisible(gc())

x[, `:=`(sac = sac_group(sac_code), dmt3 = dmt3(pccf_DMT, pccf_DMTDIFF))]
x <- x[!(ref_source == "blockface" & dmt3 == "Urban LC")]           # main sample
msg(format(nrow(x), big.mark = ","), "addresses in main sample")

# 50/50 split by postal code: all addresses sharing a postal code fall in the same fold
set.seed(SEED)
pc <- unique(x$PCODE)
fold_of <- setNames(sample(rep(1:2, length.out = length(pc))), pc)
x[, fold := fold_of[PCODE]]
msg("Fold sizes:", paste(x[, .N, keyby = fold]$N, collapse = " / "))

# ---- DA transition table -------------------------------------------------------------

tr <- x[, .(n = .N), by = .(fold, prov, sac, dmt3, ref_DAUID, pccf_DAUID)]
saveRDS(tr, file.path(OUT, "da_transitions.rds"))
msg("DA transitions:", format(nrow(tr), big.mark = ","), "rows")

# Postal-code-level transitions: used to estimate the design effect of clustering within
# postal codes (Bayesian bootstrap over postal codes in ge_build_calibration)
ptr <- x[, .(n = .N), by = .(fold, PCODE, prov, sac, dmt3, ref_DAUID, pccf_DAUID)]
saveRDS(ptr, file.path(OUT, "pc_transitions.rds"))
msg("Postal-code transitions:", format(nrow(ptr), big.mark = ","), "rows")

# ---- Positional error by stratum -----------------------------------------------------

pe <- function(d) d[, .(N = .N, median_m = round(median(dist_m)), p90_m = round(quantile(dist_m, .9)),
                        mean_log10_m = round(mean(log10(dist_m + 1)), 3),
                        pct_same_da = round(100 * mean(ref_DAUID == pccf_DAUID), 1))]
pos <- rbind(
  x[, pe(.SD), by = .(prov, sac)][, dmt3 := "All"],
  x[, pe(.SD), by = .(prov, sac, dmt3)], fill = TRUE)
setcolorder(pos, c("prov", "sac", "dmt3"))
setorder(pos, prov, sac, dmt3)
fwrite(pos, file.path(OUT, "stratum_positional_error.csv"))

# ---- Fold-2 address records for validation ---------------------------------------------

saveRDS(x[fold == 2, .(ADDR_GUID, PCODE, prov, sac, dmt3, ref_DAUID, pccf_DAUID)],
        file.path(OUT, "test_addresses.rds"))
msg("Done")
