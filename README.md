# Postal code geocoding error and exposure misclassification in Canada

Code and derived tables for:

> Peters PA, McGaughey T. Postal Code Geocoding Error and Exposure Misclassification in
> Canada: A National Validation of PCCF+ against National Address Registers. *Preprints*
> 2026. (DOI to be added)

The study geocoded 12.4 million residential addresses from Statistics Canada's National
Address Register (NAR) with the Postal Code Conversion File Plus (PCCF+) Version 8D,
measured positional error, census-area agreement and misclassification of greenness (MODIS
EVI) and area-based deprivation (MSDI 2021), and developed methods that correct
deprivation-based estimates for postal-code geocoding error.

## Contents

| Folder | Content |
|---|---|
| `correction/` | **Correction tool**: published error tables, postal-code tables and a worked example. Start here to correct your own estimates. |
| `R/` | Functions: `geoerror.R` (correction methods), `oda_match.R` (address matching), `pccf_engine_patch.R` (run-time fix for the PCCF+ R engine) |
| `scripts/` | Pipeline that reproduces the paper (see below) |
| `outputs/` | Aggregate results from the paper: tables, figures, validation summaries, supplementary tables |
| `data/` | Input data folders (not tracked, except the Index of Remoteness and README files) |
| `config.R` | All file paths, relative to the repository root |

## Correcting your own estimates

If your data were geocoded with PCCF+ Version 8D and you use MSDI 2021 quintiles, see
[`correction/README.md`](correction/README.md). Only R and the `data.table` package are
needed:

```sh
Rscript correction/example.R
```

## Reproducing the paper

1. **Install R packages** (R 4.6 or later): `Rscript install_packages.R`
2. **Get the public data**: `Rscript scripts/00_download_data.R`. It downloads the NAR and
   ODA, checks for the files that must be downloaded by hand (Geographic Attribute File,
   dissemination area boundaries, MSDI) and prints where to get them. See
   [`data/raw/README.md`](data/raw/README.md).
3. **Get MODIS EVI**: create a free NASA Earthdata login, add `EARTHDATA_USERNAME` and
   `EARTHDATA_PASSWORD` to `~/.Renviron`, then `Rscript scripts/00_modis_appeears.R`
   (about 10 GB; AppEEARS processing takes several hours).
4. **Add PCCF+**: PCCF+ is licensed by Statistics Canada and is not included. Place your
   own copy in `data/licensed/PCCF8D/` (see [`data/licensed/README.md`](data/licensed/README.md)).
5. **Run the pipeline**: `Rscript run_all.R`, or run the scripts one at a time in the order
   listed in `run_all.R`. All scripts are run from the repository root.

Requirements: about 30 GB of free disk space and 16 GB of memory. To keep data outside the
repository, copy `config_local_example.R` to `config_local.R` and set the paths there.

| Script | Step | Paper |
|---|---|---|
| `p1_data_assembly.R` | NAR residential addresses, reference locations (building, ODA, blockface), PCCF+ 8D, DA, SAC, remoteness, MSDI | 2.1–2.4 |
| `p2_variables.R` | Positional error, EVI, deprivation variables | 2.5–2.6 |
| `p3_analysis.R` | Descriptive tables and regression models | 2.7, 3.1–3.7 |
| `p4_figures_tables.R` | Figures 1–4, Tables 1–3, S1–S4 | 3 |
| `err1_transitions.R` | Transition counts and postal-code folds | 2.8 |
| `err2_validation.R` | Validation of the stratum correction (S1–S5) | 2.8.4, Table 5 |
| `err3_tables.R` | Published error and calibration tables | Tables S5, S6 |
| `err4_example.R` | Worked example | Figure 5 |
| `err5_pc_validation.R` | Postal-code tables and validation of the postal-code correction | 2.8.3, Table 5 |
| `err6_dmt_eval.R` | Choice of strata, error by delivery mode; atypical cohorts (stratum correction) | 2.8.5, Tables 4, S9 |
| `err7_pc_atypical.R` | Atypical cohorts (postal-code correction) | Table 5 |
| `err8_pc_diagnostics.R` | MCMC convergence diagnostics | 3.8 |

The pipeline writes to `data/derived/`. The versions used in the paper are in `outputs/`
and `correction/tables/`. Simulation results depend on random seeds set in each script.

## Data and licences

- Code: MIT licence (`LICENSE`).
- Derived tables in `outputs/` and `correction/tables/`: CC BY 4.0 (`LICENSE-DATA.md`).
- Input data remain under their own terms (`data/raw/README.md`). PCCF+ is used under
  licence from Statistics Canada; only stratum-level aggregates derived from it are
  published here.

## Citation

See `CITATION.cff`.
