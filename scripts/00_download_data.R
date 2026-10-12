# =============================================================================
# 00  Download and extract the public input data into data/raw
#
# Automatic (Statistics Canada, www150):
#   NAR-202606/   National Address Register, June 2026 release (46-26-0002)
#   ODA/          Open Database of Addresses v1 (46-26-0001), 10 provinces/territories
# Included in the repository (Statistics Canada Open Licence):
#   IOR-2021/     Index of Remoteness 2021 (17-26-0001)
# Manual (the host sites block scripted downloads; the script checks for the files
# and prints instructions for any that are missing):
#   GAF-2021/          Geographic Attribute File 2021 (92-151-X): 2021_92-151_X.csv
#   lda_000b21a_e/     2021 dissemination area cartographic boundary file (Shapefile)
#   Canada2021_en/     INSPQ Material and Social Deprivation Index 2021, Canada
# Separate script (needs a free NASA Earthdata login):
#   MODIS/        scripts/00_modis_appeears.R
# Licensed (not distributed): PCCF+ Version 8D, see data/licensed/README.md
#
# Existing files are kept, so the script can be rerun after a manual download.
# Disk: about 8 GB for the files below, plus 10 GB for MODIS.
# Run: Rscript scripts/00_download_data.R
# =============================================================================

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))
options(timeout = max(7200, getOption("timeout")))

STATCAN <- "https://www150.statcan.gc.ca/n1/pub"
NAR_URL <- file.path(STATCAN, "46-26-0002", "2022001", "202606.zip")
ODA_URL <- file.path(STATCAN, "46-26-0001", "2021001", "ODA_%s_v1.zip")
ODA_PROV <- c("AB", "BC", "MB", "NB", "NS", "NT", "ON", "PE", "QC", "SK")

is_zip <- function(f) file.exists(f) && !inherits(try(utils::unzip(f, list = TRUE), silent = TRUE), "try-error")

fetch <- function(url, dest) {
  if (is_zip(dest)) return(TRUE)
  message("Downloading ", url)
  ok <- !inherits(try(utils::download.file(url, dest, mode = "wb", quiet = FALSE), silent = TRUE), "try-error")
  if (!ok || !is_zip(dest)) { unlink(dest); message("  failed: ", url); return(FALSE) }
  TRUE
}

# extract a zip to a temporary folder and move files matching `pattern` into `target`
extract_to <- function(zip, target, pattern) {
  tmp <- tempfile("unz"); dir.create(tmp)
  utils::unzip(zip, exdir = tmp)
  f <- list.files(tmp, pattern = pattern, recursive = TRUE, full.names = TRUE)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  file.copy(f, file.path(target, basename(f)), overwrite = TRUE)
  unlink(tmp, recursive = TRUE)
  length(f)
}

dl_dir <- file.path(RAW, "_downloads"); dir.create(dl_dir, recursive = TRUE, showWarnings = FALSE)

# ---- NAR ----------------------------------------------------------------------------
nar_ok <- length(list.files(file.path(RAW, "NAR-202606", "Addresses"), "^Address_\\d+.*\\.csv$")) > 0 &&
          length(list.files(file.path(RAW, "NAR-202606", "Locations"), "^Location_\\d+.*\\.csv$")) > 0
if (!nar_ok) {
  z <- file.path(dl_dir, "NAR_202606.zip")
  if (fetch(NAR_URL, z)) {
    n1 <- extract_to(z, file.path(RAW, "NAR-202606", "Addresses"), "^Address_\\d+.*\\.csv$")
    n2 <- extract_to(z, file.path(RAW, "NAR-202606", "Locations"), "^Location_\\d+.*\\.csv$")
    message("NAR: ", n1, " address files, ", n2, " location files")
  }
}

# ---- ODA ----------------------------------------------------------------------------
for (p in ODA_PROV) {
  csv <- file.path(RAW, "ODA", sprintf("ODA_%s_v1", p), sprintf("ODA_%s_v1.csv", p))
  if (file.exists(csv)) next
  z <- file.path(dl_dir, sprintf("ODA_%s_v1.zip", p))
  if (fetch(sprintf(ODA_URL, p), z)) extract_to(z, dirname(csv), sprintf("^ODA_%s_v1\\.csv$", p))
}

# ---- Manual downloads ----------------------------------------------------------------
manual <- list(
  list(file = file.path(RAW, "GAF-2021", "2021_92-151_X.csv"),
       what = "Geographic Attribute File, 2021 Census (catalogue 92-151-X), CSV",
       where = "https://www12.statcan.gc.ca/census-recensement/2021/geo/aip-pia/attribute-attribs/index-eng.cfm"),
  list(file = file.path(RAW, "lda_000b21a_e", "lda_000b21a_e.shp"),
       what = "Dissemination areas, cartographic boundary file, 2021, Shapefile, English (lda_000b21a_e.zip)",
       where = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?year=21"),
  list(file = file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en", "1. EquivalenceTableCanada2021_en.xlsx"),
       what = "INSPQ Material and Social Deprivation Index, Canada 2021 (.zip), English",
       where = "https://www.inspq.qc.ca/en/deprivation/material-and-social-deprivation-index"),
  list(file = file.path(RAW, "IOR-2021", "2021IR_DATABASE.csv"),
       what = "Index of Remoteness 2021 (included in the repository; restore data/raw/IOR-2021 from git)",
       where = "https://www150.statcan.gc.ca/n1/pub/17-26-0001/172600012020001-eng.htm"))
missing <- Filter(function(m) !file.exists(m$file), manual)
for (m in missing) message("\nMANUAL DOWNLOAD NEEDED: ", m$what, "\n  from: ", m$where,
                           "\n  so that this file exists: ", m$file)

# ---- Summary ---------------------------------------------------------------------------
status <- c(NAR = length(list.files(file.path(RAW, "NAR-202606"), "\\.csv$", recursive = TRUE)) > 0,
            ODA = all(file.exists(file.path(RAW, "ODA", sprintf("ODA_%s_v1", ODA_PROV), sprintf("ODA_%s_v1.csv", ODA_PROV)))),
            setNames(vapply(manual, function(m) file.exists(m$file), TRUE), c("GAF", "DA boundary", "MSDI", "IOR")),
            MODIS = length(list.files(file.path(RAW, "MODIS"), "\\.nc$")) > 0,
            `PCCF+ (licensed)` = file.exists(file.path(PCCF_DIR, "DATA_R")))
message("\nInput data status:"); for (k in names(status)) message(sprintf("  %-18s %s", k, if (status[k]) "ok" else "MISSING"))
if (!status["MODIS"]) message("  MODIS: run scripts/00_modis_appeears.R")
if (!status["PCCF+ (licensed)"]) message("  PCCF+: see data/licensed/README.md")
message("\nDownloaded zip files can be deleted: ", dl_dir)
