# =============================================================================
# ERR3  Final calibration tables for correcting MSDI estimates from postal codes
#
# Uses all addresses (both folds). For each MSDI dimension (material, social; national
# quintiles) and each stratification:
#   A  province x SAC group               (stratum = what PCCF+ assigns)
#   B  province x SAC group x delivery mode
# writes
#   calibration_<dim>_<A|B>.rds          ge_calibration objects for ge_correct()
#   classification_<dim>_<A|B>.csv       long table: stratum, t, a, n, alpha,
#                                        P(assigned | true), P(true | assigned)
#   stratum_error_<dim>_<A|B>.csv        per stratum: N, % misclassified, % off >= 2
#                                        quintiles, weighted kappa, score bias, score
#                                        error SD, correlation, calibration slope/intercept,
#                                        residual SD (error term), positional error
# and an Excel workbook with all tables (data/derived/err/postal_code_error_tables.xlsx).
# Run: Rscript scripts/err3_tables.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl); library(writexl) })
ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
source(file.path(ROOT, "R", "geoerror.R"))
OUT <- file.path(DERIVED, "err")
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

msdi <- as.data.table(read_excel(file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en",
                                           "1. EquivalenceTableCanada2021_en.xlsx"), sheet = "Data"))
index <- list(
  material = msdi[, .(DAUID = as.character(DA), q = as.integer(QuintMat), score = SCOREMAT)][!is.na(q)],
  social   = msdi[, .(DAUID = as.character(DA), q = as.integer(QuintSoc), score = SCORESOC)][!is.na(q)])
strata_sets <- list(A = c("prov", "sac"), B = c("prov", "sac", "dmt3"), C = c("prov", "dmt3"))

tr  <- readRDS(file.path(OUT, "da_transitions.rds"))
ptr <- readRDS(file.path(OUT, "pc_transitions.rds"))      # postal-code level: clustering design effect
set.seed(20261009)
pos <- fread(file.path(OUT, "stratum_positional_error.csv"), colClasses = list(character = "prov"))
# positional error for province x DMT3 strata, from the DMT evaluation frame (err6)
fr <- readRDS(file.path(OUT, "dmt_eval", "frame.rds"))
pos_c <- fr[, .(median_m = round(10^median(ldist) - 1), p90_m = round(10^quantile(ldist, 0.9) - 1),
                pct_same_da = round(100 * mean(same_da), 1)), by = .(prov, dmt3)]
rm(fr); invisible(gc())
prov_names <- c("10" = "NL", "11" = "PE", "12" = "NS", "13" = "NB", "24" = "QC", "35" = "ON", "46" = "MB",
                "47" = "SK", "48" = "AB", "59" = "BC", "60" = "YT", "61" = "NT", "62" = "NU")

book <- list()
for (dim in names(index)) for (m in names(strata_sets)) {
  st <- strata_sets[[m]]
  cal <- ge_build_calibration(tr, index[[dim]], st, clusters = ptr)
  saveRDS(cal, file.path(OUT, sprintf("calibration_%s_%s.rds", dim, m)))
  cls <- cal$class[, c(st, "t", "a", "n", "deff", "alpha", "p_a_given_t", "p_t_given_a"), with = FALSE]
  fwrite(cls, file.path(OUT, sprintf("classification_%s_%s.csv", dim, m)))

  se <- cal$summary[cal$score[, !"N"], on = st]
  pj <- switch(m, A = pos[dmt3 == "All", !"dmt3"], B = pos[dmt3 != "All"], C = pos_c)
  se <- pj[, c(st, "median_m", "p90_m", "pct_same_da"), with = FALSE][se, on = st]
  se[, province := prov_names[prov]]
  setcolorder(se, c("province", st))
  num <- names(se)[vapply(se, is.double, TRUE)]
  se[, (num) := lapply(.SD, signif, 4), .SDcols = num]
  setorderv(se, st)
  fwrite(se, file.path(OUT, sprintf("stratum_error_%s_%s.csv", dim, m)))
  book[[sprintf("error_%s_%s", dim, m)]] <- se
  book[[sprintf("classif_%s_%s", dim, m)]] <- cls
  msg(dim, m, ":", nrow(se), "strata")
}

readme <- data.table(sheet = names(book), content = ifelse(grepl("^error", names(book)),
  "Per-stratum error summary: N addresses; % misclassified; % off by >=2 quintiles; linear weighted kappa; score bias (assigned - true), error SD, correlation; regression-calibration slope, intercept and residual SD (per-record error term); positional error median/p90 (m); % same DA.",
  "Classification table: stratum, true quintile t, assigned quintile a, address count n, design effect of clustering within postal codes (deff, by stratum and t), Dirichlet alpha ((n + k x national P) / deff), P(assigned | true), P(true | assigned)."))
readme <- rbind(data.table(sheet = "Strata", content = "A = province x SAC group of the PCCF+ assignment; B = A x PCCF+ delivery mode (Urban LC: DMT A,B,E,G; Rural: W,H,T; Other); C = province x PCCF+ delivery mode. B is recommended for the stratum correction. Built from 12,412,300 NAR residential addresses geocoded with PCCF+ 8D; MSDI 2021 national quintiles."), readme)
write_xlsx(c(list(README = readme), book), file.path(OUT, "postal_code_error_tables.xlsx"))
msg("Done")
