# pccf_nar_validation

Code, tables and an R correction tool for a national validation of PCCF+ postal code
geocoding against the National Address Register (12.4 million addresses). Reports
positional error and deprivation misclassification by delivery mode, and corrects
MSDI-based rates and rate ratios for geocoding error.

> Peters PA, McGaughey T. Postal Code Geocoding Error and Exposure Misclassification in
> Canada: A National Validation of PCCF+ Against National Address Registers. *Preprints*
> 2026 (posted 13 October 2026). https://doi.org/10.20944/preprints202610.1204.v1

## Description

Canadian health studies commonly assign area-based exposures from residential postal codes
using Statistics Canada's Postal Code Conversion File Plus (PCCF+). The study geocoded 12.4
million residential addresses from the National Address Register (NAR) with PCCF+ Version 8D
and compared each PCCF+ point with the address's building location. Error depended mainly on
delivery mode. Median positional error was 65 m for urban houses served by letter carriers
and 3.1 km for rural postal codes, and the material deprivation quintile of the 2021 Material
and Social Deprivation Index (MSDI) was misclassified for 5.0% and 46.2% of these addresses,
respectively.

The `correction/` folder lets researchers report and correct geocoding error in their own
data without rerunning the validation. It provides error measures by province, Statistical
Area Classification (SAC) type and delivery mode, classification probabilities for a stratum
correction of aggregate counts, and the distribution of true deprivation quintiles within each
postal code for a record-level correction. The functions in `R/geoerror.R` return corrected
rates and rate ratios by true quintile with simulation or credible intervals, and
`correction/example.R` demonstrates both methods. The correction needs only R and the
`data.table` package, along with the province, SAC type, delivery mode type and postal code
that PCCF+ already reports for each record.

The `scripts/` folder and `run_all.R` reproduce every table and figure in the paper. The
`outputs/` folder holds the aggregate results reported in the paper, including
Supplementary Tables S1 to S9.

The public inputs (NAR, Open Database of Addresses, Geographic Attribute File, boundary
files, MSDI and MODIS EVI) can be downloaded with the included scripts. PCCF+ is licensed
by Statistics Canada and is not included, so reproducing the full pipeline requires a
licensed copy, whereas the correction tables can be used without one. The tables apply to
PCCF+ Version 8D (December 2025 postal codes) and the 2021 MSDI.

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

See `CITATION.cff`, or cite the preprint above (https://doi.org/10.20944/preprints202610.1204.v1).
