# =============================================================================
# P4 Figures, manuscript tables and sensitivity tables
# PCCF+ geocoding validation paper (preprint target 2026-10-23)
#
# Inputs:  data/derived/p2/p2_analysis.rds, data/derived/p3/*.csv, data/derived/p1_reference_*.csv
# Outputs: data/derived/p4/figures/fig<N>_*.png (300 dpi) and .pdf
#          data/derived/p4/tables/manuscript_tables.xlsx (+ one CSV per table)
#
# Main sample and groupings follow scripts/p3_analysis.R: urban letter-carrier
# addresses (DMT A, B, E, G) with a blockface reference are excluded.
#
# Figures
#   1  Cumulative distribution of positional error by delivery mode
#   2  Delivery mode x rurality: median error (m) and material-quintile misclassification (%)
#   3  EVI (3x3, 2025-26) at the reference vs PCCF+ point, urban letter carrier vs rural postal code
#   4  Material-quintile agreement (row %), urban letter carrier vs rural postal code
# Tables
#   1  Sample by reference source, delivery mode and rurality
#   2  Positional error and area agreement by delivery mode and rurality
#   3  Exposure misclassification (EVI, MSDI) by delivery mode
#   S1 Reference-point accuracy by province (ODA and blockface vs NAR building point)
#   S2 Sensitivity: key results by reference-source sample
#   S3 Sensitivity: EVI window and year
#   S4 Adjusted predictions by delivery mode x rurality (P3 models)
#
# Colours: reference categorical palette (validated: adjacent CVD dE >= 9.1;
# aqua/yellow < 3:1 contrast, so series are direct-labelled and in a legend).
# Run: Rscript scripts/p4_figures_tables.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(writexl)
})

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
DER  <- DERIVED
FIG  <- file.path(DER, "p4", "figures")
TAB  <- file.path(DER, "p4", "tables")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# ---- Style ------------------------------------------------------------------------

ink <- list(primary = "#0b0b0b", secondary = "#52514e", muted = "#898781",
            grid = "#e1e0d9", axis = "#c3c2b7", surface = "#fcfcfb", na = "#f0efec")
series <- c("Urban letter carrier" = "#2a78d6", "Rural postal code" = "#eb6834",
            "Rural route / suburban" = "#1baf7a", "PO box / general delivery" = "#eda100")
seq_blue   <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
seq_orange <- c("#fde4d8", "#f5a47f", "#eb6834", "#b84a1c", "#7a2e0e")

theme_paper <- function(base = 9) {
  theme_minimal(base_size = base, base_family = "sans") +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          panel.grid.major = element_line(colour = ink$grid, linewidth = 0.3),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour = ink$axis, linewidth = 0.3),
          axis.text = element_text(colour = ink$secondary),
          axis.title = element_text(colour = ink$secondary),
          plot.title = element_text(colour = ink$primary, face = "bold", size = base + 1),
          plot.subtitle = element_text(colour = ink$secondary),
          legend.text = element_text(colour = ink$secondary),
          legend.title = element_text(colour = ink$secondary),
          strip.text = element_text(colour = ink$primary, face = "bold", hjust = 0))
}

# Each figure is also saved as <name>_ms.png without its title and subtitle, for the
# manuscript, where the caption carries that text.
save_fig <- function(p, name, w, h, p_ms = p + labs(title = NULL, subtitle = NULL), h_ms = h - 0.5) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  ggsave(file.path(FIG, paste0(name, "_ms.png")), p_ms, width = w, height = h_ms, dpi = 300, bg = "white")
  msg("Figure", name)
}
FIGURES_ONLY <- nzchar(Sys.getenv("P4_FIGURES_ONLY"))   # skip tables when re-rendering figures

# ---- Data (main sample, as in P3) -------------------------------------------------------

msg("Reading P2")
x <- readRDS(file.path(DER, "p2", "p2_analysis.rds"))
x <- x[pccf_coded == TRUE, .(ref_source, delivery_mode, rurality_sac, dist_m,
                             evi_ref_b1_mean, evi_pccf_b1_mean, ref_QuintMat, pccf_QuintMat)]
invisible(gc())
x[, urban_lc := delivery_mode %in% c("Urban letter carrier: house", "Urban letter carrier: apartment",
                                     "Business / large volume receiver")]
x <- x[!(ref_source == "blockface" & urban_lc)]
x[, dm4 := fcase(urban_lc, "Urban letter carrier",
                 delivery_mode == "Rural postal code", "Rural postal code",
                 delivery_mode == "Rural route / suburban service", "Rural route / suburban",
                 delivery_mode == "PO box / general delivery / other", "PO box / general delivery",
                 default = NA_character_)]
x[, dm4 := factor(dm4, names(series))]
msg(format(nrow(x), big.mark = ","), "addresses in main sample")

rd <- function(f) fread(file.path(DER, "p3", f))
dm_short <- c("Urban letter carrier: house" = "Urban LC: house",
              "Urban letter carrier: apartment" = "Urban LC: apartment",
              "Business / large volume receiver" = "Business / LVR",
              "Rural postal code" = "Rural postal code",
              "Rural route / suburban service" = "Rural route / suburban",
              "PO box / general delivery / other" = "PO box / general delivery",
              "Unknown / retired" = "Unknown / retired")
rur_levels <- c("CMA", "CA", "Strong MIZ", "Moderate MIZ", "Weak MIZ", "No MIZ", "Territories")

# ---- Figure 1: ECDF of positional error ---------------------------------------------------

probs <- seq(0.001, 0.999, by = 0.001)
ec <- x[!is.na(dm4), .(dist_m = quantile(pmax(dist_m, 1), probs), p = probs), by = dm4]
n_dm <- x[!is.na(dm4), .N, by = dm4]
lab <- ec[round(p, 3) == 0.5][n_dm, on = "dm4"]
fmt_d <- function(d) ifelse(d >= 1000, paste0(round(d / 1000, 1), " km"), paste0(round(d), " m"))
# medians sit on the same 50% line, so they are named in the legend rather than labelled in place
leg <- setNames(lab[, sprintf("%s (median %s, n = %s)", dm4, fmt_d(dist_m), comma(N))], as.character(lab$dm4))
ref_lines <- data.table(d = c(100, 1000, 5000))
f1 <- ggplot(ec, aes(dist_m, p, colour = dm4)) +
  geom_vline(data = ref_lines, aes(xintercept = d), colour = ink$grid, linewidth = 0.4) +
  geom_line(linewidth = 0.7, lineend = "round") +
  geom_point(data = lab, aes(fill = dm4), size = 2.4, shape = 21, colour = "white", stroke = 0.7) +
  scale_x_log10(labels = function(b) ifelse(b >= 1000, paste0(b / 1000, " km"), paste0(b, " m")),
                breaks = c(1, 10, 100, 1000, 10000, 100000), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(labels = percent, breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.01))) +
  scale_colour_manual(values = series, labels = leg[names(series)], name = NULL,
                      guide = guide_legend(ncol = 2, override.aes = list(linewidth = 1.2))) +
  scale_fill_manual(values = series, guide = "none") +
  coord_cartesian(xlim = c(1, 2e5)) +
  labs(title = "Positional error of PCCF+ by delivery mode",
       subtitle = "Cumulative share of addresses within a given distance of the reference location; dots mark medians",
       x = "Distance from reference location to PCCF+ point (log scale)", y = "Cumulative share of addresses") +
  theme_paper() + theme(legend.position = "top", legend.justification = "left")
save_fig(f1, "fig1_error_ecdf", 7, 4.6)

# ---- Figure 2: delivery mode x rurality heatmaps ------------------------------------------

dm_rows <- c("Urban letter carrier: house", "Urban letter carrier: apartment", "Rural postal code",
             "Rural route / suburban service", "PO box / general delivery / other")
heat_data <- function(d, value) {
  d <- d[sample == "main" & delivery_mode %in% dm_rows]
  d[, v := get(value)]
  d[small_cell == TRUE, v := NA_real_]
  # complete the grid so combinations with no addresses show as "none"
  g <- CJ(delivery_mode = dm_rows, rurality = rur_levels)
  d <- d[, .(delivery_mode, rurality, v, small_cell)][g, on = .(delivery_mode, rurality)]
  d[, status := fifelse(is.na(small_cell), "none", fifelse(small_cell, "small", "ok"))]
  d[, `:=`(dm = factor(dm_short[delivery_mode], rev(dm_short[dm_rows])),
           rurality = factor(rurality, rur_levels))]
  d
}
cell_note <- function(status) fifelse(status == "none", "none", "n<100")
lum_text <- function(hex) ifelse(colSums(col2rgb(hex) * c(0.299, 0.587, 0.114)) / 255 > 0.55, ink$primary, "white")

e2 <- heat_data(rd("t2_error_dm_x_rurality.csv"), "median_m")
e2[, label := fifelse(is.na(v), cell_note(status), fifelse(v >= 1000, sprintf("%.1f km", v / 1000), sprintf("%d m", as.integer(v))))]
pal_e <- gradient_n_pal(seq_blue)
e2[, fillhex := fifelse(is.na(v), ink$na, pal_e(squish(rescale(log10(pmax(v, 1)), from = log10(c(30, 10000))))))]
e2[, txt := fifelse(is.na(v), ink$muted, lum_text(fillhex))]
p2a <- ggplot(e2, aes(rurality, dm)) +
  geom_tile(aes(fill = fillhex), colour = "white", linewidth = 1) +
  geom_text(aes(label = label, colour = txt), size = 2.5) +
  scale_fill_identity() + scale_colour_identity() +
  labs(title = "A  Median positional error", x = NULL, y = NULL) +
  theme_paper() + theme(panel.grid.major = element_blank(), axis.line.x = element_blank())

m2 <- heat_data(rd("t4_msdi_dm_x_rurality.csv")[dimension == "material"], "pct_misclassified")
m2[, label := fifelse(is.na(v), cell_note(status), sprintf("%.0f%%", v))]
pal_m <- gradient_n_pal(seq_orange)
m2[, fillhex := fifelse(is.na(v), ink$na, pal_m(squish(rescale(v, from = c(0, 60)))))]
m2[, txt := fifelse(is.na(v), ink$muted, lum_text(fillhex))]
p2b <- ggplot(m2, aes(rurality, dm)) +
  geom_tile(aes(fill = fillhex), colour = "white", linewidth = 1) +
  geom_text(aes(label = label, colour = txt), size = 2.5) +
  scale_fill_identity() + scale_colour_identity() +
  labs(title = "B  Material deprivation quintile misclassified", x = "Rurality of reference location (SAC)", y = NULL) +
  theme_paper() + theme(panel.grid.major = element_blank(), axis.line.x = element_blank())

f2 <- (p2a / p2b) + plot_annotation(
  title = "Error depends more on delivery mode than on rurality",
  subtitle = "Main sample; n<100 = fewer than 100 addresses (not shown); none = no addresses",
  theme = theme(plot.title = element_text(face = "bold", colour = ink$primary, size = 10),
                plot.subtitle = element_text(colour = ink$secondary, size = 9)))
save_fig(f2, "fig2_dm_x_rurality", 7.5, 5.8, p_ms = p2a / p2b)

# ---- Figure 3: EVI agreement ------------------------------------------------------------------

ev <- x[dm4 %in% c("Urban letter carrier", "Rural postal code") &
          !is.na(evi_ref_b1_mean) & !is.na(evi_pccf_b1_mean)]
ev[, `:=`(bx = round(evi_ref_b1_mean, 2), by = round(evi_pccf_b1_mean, 2))]
eb <- ev[, .N, by = .(dm4, bx, by)][, share := N / sum(N), by = dm4]
eb <- eb[bx >= -0.05 & bx <= 0.8 & by >= -0.05 & by <= 0.8]
stats3 <- ev[, .(r = cor(evi_ref_b1_mean, evi_pccf_b1_mean),
                 mae = mean(abs(evi_pccf_b1_mean - evi_ref_b1_mean)),
                 bias = mean(evi_pccf_b1_mean - evi_ref_b1_mean), N = .N), by = dm4]
stats3[, label := sprintf("n = %s\nr = %.2f\nMAE = %.3f\nbias = %+.3f", comma(N), r, mae, bias)]
f3 <- ggplot(eb, aes(bx, by)) +
  geom_tile(aes(fill = share)) +
  geom_abline(slope = 1, intercept = 0, colour = ink$muted, linewidth = 0.4) +
  geom_label(data = stats3, aes(x = -0.03, y = 0.78, label = label), hjust = 0, vjust = 1,
             size = 2.5, colour = ink$primary, lineheight = 0.95, fill = "white",
             label.size = 0, label.padding = unit(0.25, "lines")) +
  facet_wrap(~dm4) +
  scale_fill_gradientn(colours = seq_blue, trans = "log10", name = "Share of\naddresses",
                       labels = function(b) paste0(format(signif(100 * b, 1), scientific = FALSE,
                                                                     drop0trailing = TRUE, trim = TRUE), "%")) +
  coord_equal(xlim = c(-0.05, 0.8), ylim = c(-0.05, 0.8), expand = FALSE) +
  labs(title = "Greenness at the PCCF+ point versus the reference location",
       subtitle = paste("MODIS EVI, growing-season mean 2025-26, 3 x 3 cell window (about 700 m); line shows agreement.",
                        "Horizontal bands: many addresses share one PCCF+ point", sep = "\n"),
       x = "EVI at reference location", y = "EVI at PCCF+ point") +
  theme_paper() + theme(panel.grid.major = element_blank(), panel.spacing = unit(1, "lines"))
save_fig(f3, "fig3_evi_agreement", 7, 4)

# ---- Figure 4: MSDI confusion matrices --------------------------------------------------------

cm <- x[dm4 %in% c("Urban letter carrier", "Rural postal code") & !is.na(ref_QuintMat) & !is.na(pccf_QuintMat),
        .N, by = .(dm4, ref = ref_QuintMat, pccf = pccf_QuintMat)]
cm[, pct := 100 * N / sum(N), by = .(dm4, ref)]
cm[, fillhex := gradient_n_pal(seq_blue)(rescale(pct, from = c(0, 100)))]
cm[, txt := lum_text(fillhex)]
f4 <- ggplot(cm, aes(factor(pccf), factor(ref, 5:1))) +
  geom_tile(aes(fill = fillhex), colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.0f", pct), colour = txt), size = 2.6) +
  scale_fill_identity() + scale_colour_identity() +
  facet_wrap(~dm4) + coord_equal() +
  labs(title = "Material deprivation quintile: reference location versus PCCF+ point",
       subtitle = "Row percentages; quintile 1 = least deprived. Diagonal = correctly classified",
       x = "Quintile assigned at PCCF+ point", y = "Quintile at reference location") +
  theme_paper() + theme(panel.grid.major = element_blank(), axis.line.x = element_blank())
save_fig(f4, "fig4_msdi_confusion", 7, 3.8)

# ---- Tables -------------------------------------------------------------------------------------

if (FIGURES_ONLY) quit(save = "no")

tabs <- list()

t1 <- rd("t1_sample.csv")[sample == "main"]
t1a <- t1[, .(N = sum(N)), by = .(ref_source)][, `:=`(group = "Reference source", level = ref_source)]
t1b <- t1[, .(N = sum(N)), by = .(delivery_mode)][, `:=`(group = "Delivery mode", level = delivery_mode)]
t1c <- t1[, .(N = sum(N)), by = .(rurality)][, `:=`(group = "Rurality (SAC)", level = rurality)]
tab1 <- rbind(t1a[, .(group, level, N)], t1b[, .(group, level, N)], t1c[, .(group, level, N)])
tab1[, pct := round(100 * N / sum(N), 1), by = group]
tabs$Table1_sample <- tab1

e_dm  <- rd("t2_error_delivery_mode.csv")[sample == "main"][, .(group = "Delivery mode", level = delivery_mode, N,
            median_m, IQR_m = paste0(q1_m, "-", q3_m), p90_m, pct_le_100m, pct_gt_5km, pct_same_da, pct_same_csd)]
e_rur <- rd("t2_error_rurality.csv")[sample == "main"][, .(group = "Rurality (SAC)", level = rurality, N,
            median_m, IQR_m = paste0(q1_m, "-", q3_m), p90_m, pct_le_100m, pct_gt_5km, pct_same_da, pct_same_csd)]
e_all <- rd("t2_error_overall.csv")[sample == "main"][, .(group = "All", level = "All addresses", N,
            median_m, IQR_m = paste0(q1_m, "-", q3_m), p90_m, pct_le_100m, pct_gt_5km, pct_same_da, pct_same_csd)]
tabs$Table2_error <- rbind(e_all, e_dm, e_rur)

v3 <- rd("t3_evi_delivery_mode.csv")[sample == "main" & window == "b1", .(level = delivery_mode, EVI_N = N,
            EVI_ref = mean_ref, EVI_pccf = mean_pccf, EVI_bias = bias, EVI_MAE = mae, EVI_r = pearson_r)]
m3 <- dcast(rd("t4_msdi_delivery_mode.csv")[sample == "main"], delivery_mode ~ dimension,
            value.var = c("pct_misclassified", "pct_off_2plus", "kappa_linear"))
setnames(m3, "delivery_mode", "level")
tabs$Table3_exposure <- merge(v3, m3, by = "level", sort = FALSE)

tabs$TableS1_reference_accuracy <- fread(file.path(DER, "p1_reference_validation.csv"))

s2e <- rd("t2_error_delivery_mode.csv")[delivery_mode %in% c("Urban letter carrier: house", "Rural postal code"),
                                        .(sample, delivery_mode, N, median_m, p90_m, pct_same_da)]
s2m <- rd("t4_msdi_delivery_mode.csv")[dimension == "material" & delivery_mode %in% c("Urban letter carrier: house", "Rural postal code"),
                                       .(sample, delivery_mode, mat_misclassified = pct_misclassified, kappa_linear)]
s2v <- rd("t3_evi_delivery_mode.csv")[window == "b1" & delivery_mode %in% c("Urban letter carrier: house", "Rural postal code"),
                                      .(sample, delivery_mode, EVI_MAE = mae, EVI_bias = bias)]
tabs$TableS2_by_reference_source <- Reduce(function(a, b) merge(a, b, by = c("sample", "delivery_mode"), all = TRUE),
                                           list(s2e, s2m, s2v))

# EVI window x year sensitivity (main sample)
msg("Sensitivity: EVI window and year")
rm(x, ev, eb, cm); invisible(gc())
y <- readRDS(file.path(DER, "p2", "p2_analysis.rds"))
y <- y[pccf_coded == TRUE]
y[, urban_lc := delivery_mode %in% c("Urban letter carrier: house", "Urban letter carrier: apartment",
                                     "Business / large volume receiver")]
y <- y[!(ref_source == "blockface" & urban_lc)]
y[, grp := fcase(urban_lc, "Urban letter carrier", delivery_mode == "Rural postal code", "Rural postal code",
                 default = NA_character_)]
s3 <- rbindlist(lapply(c("pix", "b1", "b2"), function(w) rbindlist(lapply(c("2025", "2026", "mean"), function(yr) {
  r <- paste0("evi_ref_", w, "_", yr); p <- paste0("evi_pccf_", w, "_", yr)
  y[!is.na(grp) & !is.na(get(r)) & !is.na(get(p)),
    .(window = w, year = yr, N = .N, bias = round(mean(get(p) - get(r)), 4),
      MAE = round(mean(abs(get(p) - get(r))), 4), r = round(cor(get(r), get(p)), 3)), by = grp]
}))))
rm(y); invisible(gc())
tabs$TableS3_evi_window_year <- s3[order(grp, window, year)]

tabs$TableS4_adjusted_predictions <- rd("m_adjusted_predictions.csv")[sample == "main"]

for (t in names(tabs)) fwrite(tabs[[t]], file.path(TAB, paste0(t, ".csv")))
write_xlsx(tabs, file.path(TAB, "manuscript_tables.xlsx"))
msg("Tables written:", paste(names(tabs), collapse = ", "))
