# =============================================================================
# P3 Analysis: positional error, EVI and MSDI misclassification
# PCCF+ geocoding validation paper (preprint target 2026-10-23)
#
# Input:  data/derived/p2/p2_analysis.rds
# Output: data/derived/p3/*.csv and data/derived/p3/p3_results.xlsx (one sheet per table)
#
# Samples (column `sample` in every table)
#   main       all reference sources, EXCLUDING urban letter-carrier addresses
#              (DMT A, B, E, G) whose reference is the NAR blockface point: PCCF+
#              places those postal codes on block-face points too, so the comparison
#              is partly circular (decision 2026-10-04)
#   building   NAR building point only
#   oda        ODA point only
#   blockface  NAR blockface point, excluding urban letter carrier (as in main)
#   urban_blockface_excluded   the records left out of main (reported, not modelled)
#
# Descriptives (by delivery mode, rurality, delivery mode x rurality, province,
# remoteness decile): positional error distribution, DA/CSD agreement; EVI bias,
# MAE, RMSE, correlation; MSDI misclassification, direction, weighted kappa.
#
# Models (main and building samples; delivery mode x rurality + province):
#   m1  log10(distance + 1)            OLS
#   m2  |EVI difference|, 3x3 window   OLS
#   m3  material quintile misclassified  logistic
#   m4  social quintile misclassified    logistic
# Predictors are all categorical, so each model is fitted exactly on cell-level
# sufficient statistics (counts for logistic; n, mean and within-cell SS for OLS).
# Verified on PEI against individual-level lm()/glm(): identical estimates and SEs.
# Delivery mode x rurality cells with < MIN_CELL records are not modelled (sparse
# cells give non-finite logistic estimates); they keep their descriptive results.
# Standard errors assume independent records (no DA clustering; see limitations).
#
# Run: Rscript scripts/p3_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(writexl)
})

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
IN   <- Sys.getenv("P3_IN",  file.path(DERIVED, "p2", "p2_analysis.rds"))  # env overrides for testing
OUT  <- Sys.getenv("P3_OUT", file.path(DERIVED, "p3"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
MIN_CELL <- 100   # cells with fewer records are reported but flagged

# ---- Data ---------------------------------------------------------------------

msg("Reading P2")
x <- readRDS(IN)
x <- x[pccf_coded == TRUE, c("ref_source", "delivery_mode", "rurality_sac", "PROV_CODE", "ref_IOR",
                             "dist_m", "same_da", "same_csd",
                             paste0("evi_", rep(c("ref_", "pccf_"), each = 3), c("pix", "b1", "b2"), "_mean"),
                             "evi_diff_b1", "ref_QuintMat", "pccf_QuintMat", "ref_QuintSoc", "pccf_QuintSoc"),
       with = FALSE]
invisible(gc())

dm_levels <- c("Urban letter carrier: house", "Urban letter carrier: apartment",
               "Business / large volume receiver", "Rural postal code",
               "Rural route / suburban service", "PO box / general delivery / other",
               "Unknown / retired")
x[, delivery_mode := factor(delivery_mode, dm_levels)]
x[, urban_lc := delivery_mode %in% dm_levels[1:3]]

rur_levels <- c("CMA", "CA", "Strong MIZ", "Moderate MIZ", "Weak MIZ", "No MIZ", "Territories")
x[, rurality := factor(fcase(rurality_sac %in% c("CA (tracted)", "CA (non-tracted)"), "CA",
                             rurality_sac == "Territories (outside CA)", "Territories",
                             default = as.character(rurality_sac)), rur_levels)]
x[, province := factor(PROV_CODE)]
# address-weighted deciles of the CSD-level index (ties can merge adjacent deciles)
ior_br <- unique(quantile(x$ref_IOR, 0:10 / 10, na.rm = TRUE))
x[, ior_decile := cut(ref_IOR, ior_br, include.lowest = TRUE, dig.lab = 3)]

samples <- list(
  main      = x[, !(ref_source == "blockface" & urban_lc)],
  building  = x[, ref_source == "building"],
  oda       = x[, ref_source == "oda"],
  blockface = x[, ref_source == "blockface" & !urban_lc],
  urban_blockface_excluded = x[, ref_source == "blockface" & urban_lc])
msg("Samples:", paste(names(samples), vapply(samples, sum, 1L), sep = " = ", collapse = " | "))

by_sample <- function(fun, by) rbindlist(lapply(names(samples), function(s)
  cbind(sample = s, fun(x[samples[[s]]], by))), fill = TRUE)

flag_small <- function(d) d[, small_cell := N < MIN_CELL][]

# ---- Descriptives ---------------------------------------------------------------

err_desc <- function(d, by) flag_small(d[, .(
  N = .N,
  median_m = round(median(dist_m)), q1_m = round(quantile(dist_m, .25)),
  q3_m = round(quantile(dist_m, .75)), p90_m = round(quantile(dist_m, .90)),
  p95_m = round(quantile(dist_m, .95)), mean_m = round(mean(dist_m)),
  pct_le_100m = round(100 * mean(dist_m <= 100), 1), pct_le_1km = round(100 * mean(dist_m <= 1000), 1),
  pct_gt_5km = round(100 * mean(dist_m > 5000), 1),
  pct_same_da = round(100 * mean(same_da), 1), pct_same_csd = round(100 * mean(same_csd), 1)),
  keyby = by])

evi_desc <- function(d, by) {
  out <- rbindlist(lapply(c("pix", "b1", "b2"), function(w) {
    r <- paste0("evi_ref_", w, "_mean"); p <- paste0("evi_pccf_", w, "_mean")
    d[!is.na(get(r)) & !is.na(get(p)), .(
      window = w, N = .N,
      mean_ref = round(mean(get(r)), 3), mean_pccf = round(mean(get(p)), 3),
      bias = round(mean(get(p) - get(r)), 4), mae = round(mean(abs(get(p) - get(r))), 4),
      rmse = round(sqrt(mean((get(p) - get(r))^2)), 4),
      pearson_r = round(if (.N > 2) cor(get(r), get(p)) else NA_real_, 3)), keyby = by]
  }))
  flag_small(out)
}

wkappa <- function(a, b, weight = c("linear", "quadratic")) {
  weight <- match.arg(weight)
  ok <- !is.na(a) & !is.na(b); if (sum(ok) < 2) return(NA_real_)
  t <- table(factor(a[ok], 1:5), factor(b[ok], 1:5)); n <- sum(t)
  d <- abs(outer(1:5, 1:5, "-")); w <- if (weight == "linear") 1 - d / 4 else 1 - d^2 / 16
  po <- sum(w * t) / n; pe <- sum(w * outer(rowSums(t), colSums(t))) / n^2
  if (pe == 1) return(NA_real_)
  (po - pe) / (1 - pe)
}

msdi_desc <- function(d, by) {
  out <- rbindlist(lapply(c(material = "Mat", social = "Soc"), function(dim) {
    r <- paste0("ref_Quint", dim); p <- paste0("pccf_Quint", dim)
    d[!is.na(get(r)) & !is.na(get(p)), .(
      dimension = if (dim == "Mat") "material" else "social", N = .N,
      pct_misclassified = round(100 * mean(get(p) != get(r)), 1),
      pct_pccf_less_deprived = round(100 * mean(get(p) < get(r)), 1),
      pct_pccf_more_deprived = round(100 * mean(get(p) > get(r)), 1),
      pct_off_2plus = round(100 * mean(abs(get(p) - get(r)) >= 2), 1),
      kappa_linear = round(wkappa(get(r), get(p), "linear"), 3),
      kappa_quadratic = round(wkappa(get(r), get(p), "quadratic"), 3)), keyby = by]
  }))
  flag_small(out)
}

groupings <- list(delivery_mode = "delivery_mode", rurality = "rurality",
                  dm_x_rurality = c("delivery_mode", "rurality"), province = "province",
                  ior_decile = "ior_decile", overall = character(0))

tables <- list()
tables$t1_sample <- by_sample(function(d, by) d[, .N, keyby = .(ref_source, delivery_mode, rurality)], NULL)
for (g in names(groupings)) {
  msg("Descriptives by", g)
  tables[[paste0("t2_error_", g)]] <- by_sample(err_desc,  groupings[[g]])
  tables[[paste0("t3_evi_", g)]]   <- by_sample(evi_desc,  groupings[[g]])
  tables[[paste0("t4_msdi_", g)]]  <- by_sample(msdi_desc, groupings[[g]])
}

# Confusion matrices (main sample): reference quintile (rows) x PCCF+ quintile (columns)
conf <- function(dim) {
  m <- x[samples$main, table(ref = get(paste0("ref_Quint", dim)), pccf = get(paste0("pccf_Quint", dim)))]
  d <- as.data.table(unclass(m), keep.rownames = "ref_quintile")
  setnames(d, -1, paste0("pccf_Q", colnames(m))); d
}
tables$t5_confusion_material <- conf("Mat")
tables$t5_confusion_social   <- conf("Soc")

# ---- Models on cell-level sufficient statistics --------------------------------

drop_sparse <- function(d) {
  d <- d[delivery_mode != "Unknown / retired" & !is.na(rurality)]
  big <- d[, .N, by = .(delivery_mode, rurality)][N >= MIN_CELL]
  d[big, on = .(delivery_mode, rurality), nomatch = NULL]
}

model_cells <- function(d, y) {
  d <- drop_sparse(d[!is.na(get(y))])
  d[, .(n = .N, ybar = mean(get(y)), ss_within = sum((get(y) - mean(get(y)))^2)),
    by = .(delivery_mode, rurality, province)
  ][, `:=`(delivery_mode = droplevels(delivery_mode), rurality = droplevels(rurality),
           province = droplevels(province))][]
}

# Exact OLS from cell means: estimates from weighted lm; residual SS = within + between
ols_cells <- function(cells, f) {
  fit <- lm(f, data = cells, weights = n)
  rss <- sum(cells$ss_within) + sum(cells$n * residuals(fit)^2)
  df  <- sum(cells$n) - length(coef(fit)[!is.na(coef(fit))])
  s2  <- rss / df
  V   <- s2 * summary(fit)$cov.unscaled
  list(fit = fit, rss = rss, df = df, s2 = s2, V = V)
}

coef_table <- function(est, V, model, sample, scale = "identity") {
  se <- sqrt(diag(V)); z <- est / se
  data.table(model = model, sample = sample, term = names(est), estimate = est, se = se,
             ci_low = est - 1.96 * se, ci_high = est + 1.96 * se, p = 2 * pnorm(-abs(z)),
             scale = scale)
}

# Adjusted predictions for each delivery mode x rurality cell, averaged over the
# province distribution of the modelled sample (marginal standardization)
adj_pred <- function(fit, cells, type = "response") {
  pw <- cells[, .(w = sum(n)), by = province][, w := w / sum(w)]
  grid <- CJ(delivery_mode = levels(cells$delivery_mode), rurality = levels(cells$rurality),
             province = pw$province)
  grid[, `:=`(delivery_mode = factor(delivery_mode, levels(cells$delivery_mode)),
              rurality = factor(rurality, levels(cells$rurality)),
              province = factor(province, levels(cells$province)))]
  grid[, pred := suppressWarnings(predict(fit, newdata = grid, type = type))]
  grid <- merge(grid, pw, by = "province")
  out <- grid[, .(adjusted_prediction = sum(pred * w)), by = .(delivery_mode, rurality)]
  nn <- cells[, .(N = sum(n)), by = .(delivery_mode, rurality)]
  out <- merge(out, nn, by = c("delivery_mode", "rurality"), all.x = TRUE)
  out[is.na(N), `:=`(N = 0L, adjusted_prediction = NA_real_)]   # not modelled
  out[, modelled := N > 0L][]
}

# delivery mode x rurality + province; terms with a single level in the data are dropped
model_formulas <- function(cells, y) {
  has <- function(v) nlevels(cells[[v]]) > 1L
  core <- c("delivery_mode", "rurality")[c(has("delivery_mode"), has("rurality"))]
  adj  <- if (has("province")) "province" else NULL
  inter <- if (length(core) == 2L) "delivery_mode:rurality" else NULL
  list(full = reformulate(c(core, inter, adj), y), main = reformulate(c(core, adj), y))
}

coefs <- list(); tests <- list(); preds <- list()
for (s in c("main", "building")) {
  d <- x[samples[[s]]]
  d[, `:=`(log10_dist = log10(dist_m + 1), abs_evi_b1 = abs(evi_diff_b1))]

  for (m in list(list(id = "m1_log10_distance", y = "log10_dist"),
                 list(id = "m2_abs_evi_diff_b1", y = "abs_evi_b1"))) {
    msg("Model", m$id, "|", s)
    cells <- model_cells(d, m$y)
    fm <- model_formulas(cells, "ybar")
    full <- ols_cells(cells, fm$full)
    red  <- ols_cells(cells, fm$main)
    est <- coef(full$fit); keep <- !is.na(est)
    coefs[[paste(m$id, s)]] <- coef_table(est[keep], full$V, m$id, s)
    df1 <- red$df - full$df
    Fst <- ((red$rss - full$rss) / df1) / full$s2
    tests[[paste(m$id, s)]] <- data.table(model = m$id, sample = s, test = "delivery mode x rurality",
      statistic = "F", value = Fst, df1 = df1, df2 = full$df, p = pf(Fst, df1, full$df, lower.tail = FALSE),
      N = sum(cells$n))
    p <- adj_pred(full$fit, cells)
    # back-transformed: geometric mean of (distance + 1), minus 1
    if (m$y == "log10_dist") p[, adjusted_geomean_m := round(10^adjusted_prediction - 1)]
    preds[[paste(m$id, s)]] <- cbind(model = m$id, sample = s, p)
  }

  for (m in list(list(id = "m3_material_misclassified", r = "ref_QuintMat", p = "pccf_QuintMat"),
                 list(id = "m4_social_misclassified",   r = "ref_QuintSoc", p = "pccf_QuintSoc"))) {
    msg("Model", m$id, "|", s)
    cells <- drop_sparse(d[!is.na(get(m$r)) & !is.na(get(m$p))])[,
               .(n = .N, k = sum(get(m$p) != get(m$r))), by = .(delivery_mode, rurality, province)]
    cells[, `:=`(delivery_mode = droplevels(delivery_mode), rurality = droplevels(rurality),
                 province = droplevels(province))]
    fm <- model_formulas(cells, "cbind(k, n - k)")
    full <- glm(fm$full, family = binomial, data = cells)
    red  <- glm(fm$main, family = binomial, data = cells)
    est <- coef(full); keep <- !is.na(est)
    ct <- coef_table(est[keep], vcov(full, complete = FALSE), m$id, s, scale = "log-odds")
    ct[, `:=`(odds_ratio = exp(estimate), or_low = exp(ci_low), or_high = exp(ci_high))]
    coefs[[paste(m$id, s)]] <- ct
    lr <- red$deviance - full$deviance; df1 <- red$df.residual - full$df.residual
    tests[[paste(m$id, s)]] <- data.table(model = m$id, sample = s, test = "delivery mode x rurality",
      statistic = "LR chi-square", value = lr, df1 = df1, df2 = NA_real_,
      p = pchisq(lr, df1, lower.tail = FALSE), N = sum(cells$n))
    p <- adj_pred(full, cells)
    p[, adjusted_pct := round(100 * adjusted_prediction, 1)]
    preds[[paste(m$id, s)]] <- cbind(model = m$id, sample = s, p)
  }
}
tables$m_coefficients  <- rbindlist(coefs, fill = TRUE)
tables$m_interaction_tests <- rbindlist(tests, fill = TRUE)
tables$m_adjusted_predictions <- rbindlist(preds, fill = TRUE)

# ---- Write ---------------------------------------------------------------------

for (t in names(tables)) fwrite(tables[[t]], file.path(OUT, paste0(t, ".csv")))
names(tables) <- substr(names(tables), 1, 31)                 # Excel sheet-name limit
write_xlsx(tables, file.path(OUT, "p3_results.xlsx"))
msg("Wrote", length(tables), "tables to", OUT)

msg("Key results (main sample):")
print(tables$t2_error_delivery_mode[sample == "main", .(delivery_mode, N, median_m, p90_m, pct_same_da)])
print(tables$t4_msdi_delivery_mode[sample == "main", .(dimension, delivery_mode, N, pct_misclassified, kappa_linear)])
print(tables$m_interaction_tests)
