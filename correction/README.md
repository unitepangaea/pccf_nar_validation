# Correcting deprivation estimates for postal-code geocoding error

This folder lets a study that geocoded residential postal codes with PCCF+ Version 8D
report the expected geocoding error of its MSDI 2021 quintiles and correct its estimates,
without running the validation pipeline. The methods are described in Section 2.8 of the
paper; the code is `R/geoerror.R` (requires the `data.table` package).

Run the example from the repository root: `Rscript correction/example.R`.

## Which method to use

| Your data | Method | Function | Table |
|---|---|---|---|
| Record-level postal codes | Postal-code correction (recommended) | `ge_pc_correct()` | `pc_quintile_probs_<dim>.csv.gz` |
| Counts by stratum and assigned quintile only | Stratum correction | `ge_correct()` | `calibration_<dim>_B.rds` |
| Any (to report the error) | Error measures | — | `stratum_error_<dim>_<A,B,C>.csv`, `postal_code_error_tables.xlsx` |

`<dim>` is `material` or `social`. Strata: `A` province × SAC type; `B` province × SAC type ×
delivery mode (recommended for the stratum correction); `C` province × delivery mode.

## Preparing your data

From your PCCF+ output, for each person:

- `a`: MSDI 2021 national quintile of the PCCF+ DA (1 = least deprived)
- `prov`: PCCF+ province code (`PR`, e.g. `"35"`)
- `sac`: PCCF+ SAC type grouped as `CMA` (1), `CA` (2, 3), `Strong MIZ` (4), `Moderate MIZ` (5),
  `Weak MIZ` (6), `No MIZ` (7), `Territories` (8)
- `dmt3`: PCCF+ delivery mode type grouped as `Urban LC` (A, B, E, G), `Rural` (W, H, T) or
  `Other` (J, K, M, R, X and retired codes); for DMT `Z`, use `DMTDIFF` when present
- `PCODE`: six-character postal code without a space (postal-code method)

## Stratum correction

```r
source("R/geoerror.R")
cal <- readRDS("correction/tables/calibration_material_B.rds")
obs <- my_data[, .(n = .N, events = sum(case)), by = .(prov, sac, dmt3, a)]
res <- ge_correct(obs, cal, draws = 1000, sample_n = FALSE)
res$rate_ratios     # corrected Q:Q1 rate ratios with 95% simulation intervals
res$distribution    # corrected share of the population by quintile
res$diagnostics     # EM convergence and share of the population with boundary solutions
```

Use `sample_n = TRUE` when the population counts are themselves a sample (for example a
cohort) and `FALSE` for census-based denominators. A large boundary share means that the
classification probabilities fit the observed counts poorly in those strata, and the
intervals should be interpreted with caution.

## Postal-code correction

```r
source("R/geoerror.R")
pct <- data.table::as.data.table(read.csv(gzfile("correction/tables/pc_quintile_probs_material.csv.gz"),
                                           colClasses = c(PCODE = "character")))
res <- ge_pc_correct(my_data$PCODE, my_data$a, pct, y = my_data$case, iter = 1500, burn = 500)
res$rate_ratios     # posterior means with 95% credible intervals
res$fallback_share  # share of records whose postal code is not in the table
```

Records whose postal code is not in the table keep their assigned quintile.

## Assumptions

- Misclassification is unrelated to the outcome except through the true quintile.
- Stratum method: selection into the study does not depend on the assigned quintile given
  the true quintile and the stratum.
- Postal-code method: study members are distributed across the residential addresses of
  each postal code as the addresses are, and selection does not depend on the true quintile
  within a postal code.
- Postal codes are current (December 2025 PCCF+) and correctly recorded.

## Tables (CC BY 4.0)

| File | Content |
|---|---|
| `calibration_<dim>_<A,B,C>.rds` | Calibration objects for `ge_correct()`: counts, design effects and Dirichlet parameters by stratum, true and assigned quintile |
| `classification_<dim>_<A,B,C>.csv` | The same classification probabilities as CSV |
| `stratum_error_<dim>_<A,B,C>.csv` | Error measures by stratum: % misclassified, weighted kappa, score error, calibration slope, positional error |
| `postal_code_error_tables.xlsx` | All of the above in one workbook |
| `pc_quintile_probs_<dim>.csv.gz` | For each postal code in the National Address Register (June 2026): number of residential addresses and the share of addresses in each true MSDI quintile |
