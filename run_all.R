# =============================================================================
# Runs the full pipeline in order. Each step is a separate R process, so a failure
# stops the run at that step and later steps can be rerun individually.
# Approximate times on a 2024 laptop (8 cores, 16 GB): P1 several hours (PCCF+ engine),
# P2 25 min, P3 30 min, P4 5 min, err1-err8 about 3 hours in total.
# Run from the repository root: Rscript run_all.R
# =============================================================================

steps <- c(
  "scripts/00_download_data.R",   # check inputs (downloads NAR and ODA if missing)
  "scripts/p1_data_assembly.R",   # NAR addresses, reference locations, PCCF+ 8D, linkage
  "scripts/p2_variables.R",       # distance, EVI, MSDI variables
  "scripts/p3_analysis.R",        # descriptive tables and models
  "scripts/p4_figures_tables.R",  # figures 1-4 and tables 1-3, S1-S4
  "scripts/err1_transitions.R",   # transition tables and folds for the correction methods
  "scripts/err2_validation.R",    # stratum correction validation (S1-S5)
  "scripts/err5_pc_validation.R", # postal-code correction validation and postal-code tables
  "scripts/err6_dmt_eval.R",      # choice of strata; atypical cohorts (stratum methods)
  "scripts/err7_pc_atypical.R",   # atypical cohorts (postal-code method)
  "scripts/err3_tables.R",        # published calibration and error tables
  "scripts/err4_example.R",       # Figure 5 example
  "scripts/err8_pc_diagnostics.R")# MCMC diagnostics for the postal-code method

for (s in steps) {
  message("\n==== ", s, " (", format(Sys.time(), "%H:%M"), ")")
  rc <- system2(file.path(R.home("bin"), "Rscript"), shQuote(s))
  if (rc != 0) stop("Step failed: ", s)
}
message("\nPipeline complete. Outputs are in data/derived; compare with outputs/ and correction/tables/.")
