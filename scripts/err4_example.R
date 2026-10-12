# =============================================================================
# ERR4  Worked example: correcting a deprivation gradient estimated from postal codes
#
# A synthetic cohort of 50,000 addresses is drawn from Ontario's MIZ strata (rural
# Ontario) among the address-level validation half (afold 2), with an outcome whose true
# Q5:Q1 rate ratio is 2. The study is analysed as a user would: only postal code, the
# PCCF+-assigned quintile, province, SAC and delivery mode are used.
#   Postal-code method: postal-code table from the other address half (afold 1).
#   Stratum method: province x SAC x DMT3 (calibration_material_B; all addresses; illustration only).
#
# Outputs: data/derived/err/example_rural_ontario.csv and .png
# Run: Rscript scripts/err4_example.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl); library(ggplot2) })
ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
source(file.path(ROOT, "R", "geoerror.R"))
OUT <- file.path(DERIVED, "err")
set.seed(42)

msdi <- as.data.table(read_excel(file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en",
                                           "1. EquivalenceTableCanada2021_en.xlsx"), sheet = "Data"))
idx <- msdi[, .(DAUID = as.character(DA), q = as.integer(QuintMat), score = SCOREMAT)][!is.na(q)]
cal <- readRDS(file.path(OUT, "calibration_material_B.rds"))

x <- readRDS(file.path(OUT, "pc_addresses.rds"))
pt1 <- ge_pc_table(x[afold == 1, .(n = .N), by = .(PCODE, ref_DAUID)], idx)
co <- x[afold == 2 & prov == "35" & sac %in% c("Strong MIZ", "Moderate MIZ", "Weak MIZ", "No MIZ")]
rm(x); invisible(gc())
co[idx, t := i.q, on = .(ref_DAUID = DAUID)]
co[idx, a := i.q, on = .(pccf_DAUID = DAUID)]
co <- co[!is.na(t) & !is.na(a)][sample(.N, 50000)]
co[, case := rbinom(.N, 1, 0.05 * 2^((t - 1) / 4))]

# Stratum method (aggregate data) and postal-code method (record-level data)
obs <- co[, .(n = .N, events = sum(case)), by = .(prov, sac, dmt3, a)]
rs <- ge_correct(obs, cal, draws = 1000, sample_n = TRUE, seed = 1)
rp <- ge_pc_correct(co$PCODE, co$a, pt1, y = co$case, iter = 2000, burn = 500, seed = 1)

truth <- co[, .(true_n = .N, true_rate = mean(case)), keyby = .(q = t)]
truth[, `:=`(true_pct = 100 * true_n / sum(true_n), true_rr = true_rate / true_rate[1])]
fmt <- function(e, l, h, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), e, l, h)
out <- data.table(quintile = 1:5,
  true_pct = round(truth$true_pct, 1),
  naive_pct = round(rs$distribution$naive_pct, 1),
  stratum_pct = fmt(rs$distribution$corrected_pct, rs$distribution$lo_pct, rs$distribution$hi_pct, 1),
  postalcode_pct = fmt(rp$distribution$corrected_pct, rp$distribution$lo_pct, rp$distribution$hi_pct, 1),
  true_rr = round(truth$true_rr, 2),
  naive_rr = fmt(rs$rate_ratios$naive_rr, rs$rate_ratios$naive_lo, rs$rate_ratios$naive_hi),
  stratum_rr = fmt(rs$rate_ratios$corrected_rr, rs$rate_ratios$lo_rr, rs$rate_ratios$hi_rr),
  postalcode_rr = fmt(rp$rate_ratios$corrected_rr, rp$rate_ratios$lo_rr, rp$rate_ratios$hi_rr),
  postalcode_se_log_rr = round(rp$rate_ratios$se_log_rr, 3))
fwrite(out, file.path(OUT, "example_rural_ontario.csv"))
print(out)

pd <- rbind(rs$rate_ratios[, .(q, est = naive_rr, lo = naive_lo, hi = naive_hi, method = "Naive (assigned quintile)")],
            rs$rate_ratios[, .(q, est = corrected_rr, lo = lo_rr, hi = hi_rr, method = "Stratum correction (province x SAC x DMT)")],
            rp$rate_ratios[, .(q, est = corrected_rr, lo = lo_rr, hi = hi_rr, method = "Postal-code correction")])
pd[, method := factor(method, unique(method))]
g <- ggplot(pd, aes(q, est, colour = method)) +
  geom_hline(yintercept = 1, colour = "#c3c2b7", linewidth = 0.3) +
  geom_line(data = data.table(q = 1:5, est = 2^((0:4) / 4)), aes(q, est), inherit.aes = FALSE,
            colour = "#898781", linewidth = 0.5) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.45), size = 0.3, linewidth = 0.6) +
  scale_colour_manual(values = c("#eb6834", "#1baf7a", "#2a78d6"), name = NULL) +
  scale_y_log10(breaks = c(0.8, 1, 1.25, 1.5, 2, 2.5, 3, 4)) +
  scale_x_continuous(breaks = 1:5) +
  labs(title = "Rural Ontario: rate ratio by material deprivation quintile",
       subtitle = sprintf("Synthetic cohort of 50,000 addresses; grey line = population gradient (Q5:Q1 = 2); this cohort's true Q5:Q1 = %.2f; 95%% intervals",
                          truth$true_rr[5]),
       x = "Material deprivation quintile (1 = least deprived)", y = "Rate ratio vs quintile 1 (log scale)") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top", legend.justification = "left", panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
        plot.title = element_text(face = "bold"), plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(OUT, "example_rural_ontario.png"), g, width = 6.5, height = 4.2, dpi = 300, bg = "white")
# manuscript version (Figure 5): no title or subtitle, the caption carries them
ggsave(file.path(OUT, "example_rural_ontario_ms.png"), g + labs(title = NULL, subtitle = NULL),
       width = 6.5, height = 3.8, dpi = 300, bg = "white")
cat("Wrote example table and figure\n")
