# =============================================================================
# P1 Data assembly: NAR residential addresses geocoded with PCCF+ 8D
# PCCF+ geocoding validation paper (preprint target 2026-10-23)
#
# Builds one record per residential NAR address with:
#   - a reference location (EPSG:3347) from a hierarchy, flagged in ref_source:
#       "building"  NAR building point (BG_X/BG_Y)
#       "oda"       Open Database of Addresses point matched on the civic address
#                   (R/oda_match.R; no ODA for NL, YT, NU)
#       "blockface" NAR blockface representative point (BF_REPPOINT_X/Y)
#     Addresses with none of the three are dropped. All three points are kept
#     (bg_*, oda_*, bf_*) so analyses can be run by source or restricted to buildings.
#   - the PCCF+ 8D residential-routine point and geography for its postal code
#   - delivery mode (DMT) from PCCF+ (the NAR has no delivery-mode field)
#   - 2021 DA, CSD, SAC type and Index of Remoteness at the reference location
#   - INSPQ MSDI 2021 scores and quintiles at the reference DA and the PCCF+ DA
#
# MODIS vegetation index extraction is done in P2. This script only checks the netCDF.
#
# Steps (each is resumable: outputs that already exist are skipped unless FORCE = TRUE)
#   1. nar   Extract residential NAR addresses, match ODA, assign reference point
#   2. pccf  Run PCCF+ 8D (residential routine) in chunks, each in a fresh R process.
#            Only addresses not already in a finished chunk are run, so adding
#            addresses later does not repeat earlier chunks.
#   3. link  Join PCCF+ output, DA point-in-polygon for reference points, GAF, IOR, MSDI
#   4. ndvi  Check the MODIS NDVI netCDF
#   5. combine  Stack provinces into one analysis file
#
# Run from a terminal:  Rscript scripts/p1_data_assembly.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(readxl)
})

# ---- Settings ---------------------------------------------------------------

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
OUT       <- DERIVED          # all outputs go here

NAR_DIR   <- file.path(RAW, "NAR-202606", "Addresses")
LOC_DIR   <- file.path(RAW, "NAR-202606", "Locations")
# PCCF_DIR (licensed PCCF+ files) is set in config.R
PATCH_FILE <- file.path(ROOT, "R", "pccf_engine_patch.R")  # fixes for the PCCF+ R engine
GAF_FILE  <- file.path(RAW, "GAF-2021", "2021_92-151_X.csv")
IOR_FILE  <- file.path(RAW, "IOR-2021", "2021IR_DATABASE.csv")
DA_SHP    <- file.path(RAW, "lda_000b21a_e", "lda_000b21a_e.shp")
MSDI_XLSX <- file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en",
                       "1. EquivalenceTableCanada2021_en.xlsx")
NDVI_DIR  <- file.path(RAW, "MODIS")             # AppEEARS MOD13Q1 netCDF

STEPS      <- c("nar", "pccf", "link", "ndvi", "combine")
PROVINCES  <- NULL      # NULL = all available; or pass codes: Rscript p1_data_assembly.R 11
CHUNK_SIZE <- 1e6       # addresses per PCCF+ run (keeps memory under ~8 GB)
FORCE      <- FALSE     # TRUE re-runs steps whose outputs already exist
BU_USE_RES <- c(1L)     # NAR building use: 1L residential, 2L partial residential

CRS_LAMBERT <- 3347     # Statistics Canada Lambert (NAR BG_X/BG_Y, DA boundaries)
CRS_NAD83   <- 4269     # PCCF+ LAT/LONG

PC_REGEX <- "^[ABCEGHJ-NPRSTVXY][0-9][A-Z][0-9][A-Z][0-9]$"

dir.create(file.path(OUT, "nar"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "pccf"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "linked"), recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
need <- function(f) FORCE || !file.exists(f)

source(file.path(ROOT, "R", "oda_match.R"))

# Provinces and territories available in the NAR Address files
nar_files <- list.files(NAR_DIR, pattern = "^Address_\\d+.*\\.csv$", full.names = TRUE)
nar_pr    <- sub("^Address_(\\d+).*$", "\\1", basename(nar_files))
args <- commandArgs(trailingOnly = TRUE)   # e.g. Rscript p1_data_assembly.R 10 12 13
if (length(args)) PROVINCES <- args
if (is.null(PROVINCES)) PROVINCES <- sort(unique(nar_pr))
msg("Provinces:", paste(PROVINCES, collapse = " "))

# ---- Step 1: NAR residential extract and reference point ---------------------

nar_cols <- c("ADDR_GUID", "LOC_GUID", "PROV_CODE", "APT_NO_LABEL", "CIVIC_NO", "CIVIC_NO_SUFFIX",
              "OFFICIAL_STREET_NAME", "OFFICIAL_STREET_TYPE", "OFFICIAL_STREET_DIR",
              "MAIL_POSTAL_CODE", "BU_N_CIVIC_ADD", "BU_USE",
              "BG_X", "BG_Y", "BF_REPPOINT_X", "BF_REPPOINT_Y")

extract_nar <- function(pr) {
  f_out <- file.path(OUT, "nar", paste0("nar_ref_", pr, ".rds"))
  if (!need(f_out)) return(invisible(NULL))
  files <- nar_files[nar_pr == pr]
  msg("NAR", pr, ": reading", length(files), "file(s)")
  x <- rbindlist(lapply(files, fread, select = nar_cols, encoding = "UTF-8",
                        colClasses = list(character = c("ADDR_GUID", "LOC_GUID", "PROV_CODE",
                                                        "MAIL_POSTAL_CODE", "CIVIC_NO",
                                                        "CIVIC_NO_SUFFIX"))))
  n_read <- nrow(x)
  x[, PCODE := gsub("\\s+", "", toupper(MAIL_POSTAL_CODE))]
  x <- unique(x[BU_USE %in% BU_USE_RES], by = "ADDR_GUID")
  n_res <- nrow(x)
  x <- x[grepl(PC_REGEX, PCODE)]
  n_res_valid_pc <- nrow(x)

  # 2026 interim CSD code from the Location file (used for the ODA match)
  lf <- list.files(LOC_DIR, paste0("^Location_", pr, "(_part_\\d+)?\\.csv$"),
                   full.names = TRUE, recursive = TRUE)
  loc <- unique(rbindlist(lapply(lf, fread, select = c("LOC_GUID", "CSD_CODE"),
                                 colClasses = "character")), by = "LOC_GUID")
  x <- merge(x, loc, by = "LOC_GUID", all.x = TRUE)

  setnames(x, c("BG_X", "BG_Y", "BF_REPPOINT_X", "BF_REPPOINT_Y"), c("bg_x", "bg_y", "bf_x", "bf_y"))
  x <- match_oda(x, pr)

  # Reference point hierarchy: building -> ODA -> blockface
  x[, ref_source := fcase(!is.na(bg_x) & !is.na(bg_y),   "building",
                          !is.na(oda_x) & !is.na(oda_y), "oda",
                          !is.na(bf_x) & !is.na(bf_y),   "blockface",
                          default = NA_character_)]
  x[, `:=`(ref_x = fcase(ref_source == "building", bg_x, ref_source == "oda", oda_x,
                         ref_source == "blockface", bf_x),
           ref_y = fcase(ref_source == "building", bg_y, ref_source == "oda", oda_y,
                         ref_source == "blockface", bf_y))]

  val <- ref_validation(x, pr)
  fwrite(val, file.path(OUT, "nar", paste0("ref_validation_", pr, ".csv")))

  counts <- data.table(pr = pr, n_read = n_read, n_res = n_res, n_res_valid_pc = n_res_valid_pc,
                       n_ref_building  = x[, sum(ref_source == "building", na.rm = TRUE)],
                       n_ref_oda       = x[, sum(ref_source == "oda", na.rm = TRUE)],
                       n_ref_blockface = x[, sum(ref_source == "blockface", na.rm = TRUE)],
                       n_no_ref        = x[, sum(is.na(ref_source))])
  x <- x[!is.na(ref_source)]
  x[, multi_unit := APT_NO_LABEL != "" & !is.na(APT_NO_LABEL)]
  x[, c("MAIL_POSTAL_CODE", "CIVIC_NO", "CIVIC_NO_SUFFIX", "OFFICIAL_STREET_NAME",
        "OFFICIAL_STREET_TYPE", "OFFICIAL_STREET_DIR") := NULL]
  setnames(x, "CSD_CODE", "nar_CSD_2026")
  counts[, n_extract := nrow(x)]
  saveRDS(x, f_out)
  fwrite(counts, file.path(OUT, "nar", paste0("nar_counts_", pr, ".csv")))
  msg("NAR", pr, ":", format(n_read, big.mark = ","), "read ->", format(nrow(x), big.mark = ","),
      "residential | building", counts$n_ref_building, "| ODA", counts$n_ref_oda,
      "| blockface", counts$n_ref_blockface, "| none (dropped)", counts$n_no_ref)
  if (nrow(val) && "median_m" %in% names(val))
    msg("NAR", pr, ": ODA vs building point median", val[point == "ODA", median_m], "m;",
        "blockface vs building median", val[point == "NAR blockface", median_m], "m")
}

# ---- Step 2: PCCF+ 8D in chunks --------------------------------------------

pccf_vars <- c(
  # IDs
  "CODER", "Version", "ID", "PCODE",
  # geography
  "PR", "DAUID", "CSDuid", "CSDtype", "CMA", "CMAtype", "SACcode", "SACtype",
  "LAT", "LONG", "CSIZEMIZ",
  # accuracy and delivery mode
  "SLI", "Rep_Pt_Type", "RPF", "PCtype", "DMT", "H_DMT", "DMTDIFF", "PO", "QI",
  "Source", "Link_Source", "LINK", "nCD", "nCSD", "PREC",
  # supplementary
  "Comm_Name", "ResFlag", "InstFlag")

check_pccf_packages <- function() {
  pk <- c("dplyr", "readr", "stringr", "data.table", "gtsummary", "gt")
  miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("PCCF+ needs these R packages: ", paste(miss, collapse = ", "),
         "\nInstall with: install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))",
         call. = FALSE)
}

# The PCCF+ engine runs in the global environment and deletes objects by prefix,
# so each chunk runs in its own Rscript process from a generated control file.
run_pccf_chunk <- function(tag, ids, seed) {
  f_out <- file.path(OUT, "pccf", paste0(tag, ".rds"))
  if (!need(f_out)) return(invisible(NULL))
  f_in  <- file.path(OUT, "pccf", paste0(tag, "_input.rds"))
  saveRDS(ids, f_in)
  ctrl <- file.path(OUT, "pccf", paste0(tag, "_control.R"))
  writeLines(c(
    "suppressPackageStartupMessages({library(dplyr); library(readr); library(stringr);",
    "  library(data.table); library(gtsummary); library(gt)})",
    sprintf('installDir <- "%s"', path.expand(PCCF_DIR)),
    sprintf('input_data <- readRDS("%s")', path.expand(f_in)),
    sprintf('outData <- "%s"', path.expand(file.path(OUT, "pccf"))),
    sprintf('outName <- "%s"', tag),
    "problem_file <- 1",
    # diagnostics as HTML: gt needs pandoc to write .docx
    sprintf('file_path  <- file.path(outData, "%s_geocoded_tables.html")', tag),
    sprintf('file_pathx <- file.path(outData, "%s_problem_tables.html")', tag),
    "codeVersion <- 0  # residential routine",
    sprintf("seedVal <- %d", seed),
    'Version <- "8D"',
    "set.seed(seedVal)",
    sprintf("selectedVars <- c(%s)", paste0('"', pccf_vars, '"', collapse = ", ")),
    "# timestamp the engine's step messages in the log",
    'cat <- function(...) base::cat(format(Sys.time(), "%H:%M:%S"), ...)',
    "# run the engine with random_select_by_pcode() replaced by the patched version",
    'eng <- readLines(file.path(installDir, "DATA_R", "PCCFplus_8D_do_not_edit.R"))',
    's <- grep("^random_select_by_pcode <- function", eng)',
    'e <- grep("^# 2\\\\. Function to assemble", eng) - 1',
    'stopifnot(length(s) == 1, length(e) == 1, e > s)',
    sprintf('eng <- c(eng[seq_len(s - 1)], readLines("%s"), eng[e:length(eng)])', PATCH_FILE),
    "source(exprs = parse(text = eng))",
    'base::cat("Patched selection: records with all-zero/NA weights (equal-probability fallback):",',
    '          .rsp_state$n_fallback, "| records with some NA weights:", .rsp_state$n_na_weight, "\\n")'
  ), ctrl)
  msg("PCCF+", tag, ":", format(nrow(ids), big.mark = ","), "addresses")
  rc <- system2(file.path(R.home("bin"), "Rscript"), shQuote(ctrl),
                stdout = file.path(OUT, "pccf", paste0(tag, ".log")),
                stderr = file.path(OUT, "pccf", paste0(tag, ".log")))
  if (rc != 0 || !file.exists(f_out))
    stop("PCCF+ failed for ", tag, "; see ", file.path(OUT, "pccf", paste0(tag, ".log")),
         call. = FALSE)
  unlink(f_in)
  unlink(file.path(OUT, "pccf", paste0(tag, ".csv")))   # rds copy is kept
}

# IDs already geocoded for a province: each chunk's input IDs are kept in <tag>_ids.rds
done_ids <- function(pr) {
  f <- list.files(file.path(OUT, "pccf"), pattern = paste0("^pccf_", pr, "_\\d+_ids\\.rds$"),
                  full.names = TRUE)
  f <- f[file.exists(sub("_ids\\.rds$", ".rds", f))]          # only finished chunks
  unique(unlist(lapply(f, readRDS), use.names = FALSE))
}

run_pccf <- function(pr) {
  x <- readRDS(file.path(OUT, "nar", paste0("nar_ref_", pr, ".rds")))
  ids <- x[!ADDR_GUID %chin% done_ids(pr), .(ID = ADDR_GUID, PCODE)]
  if (!nrow(ids)) return(invisible(NULL))
  done <- list.files(file.path(OUT, "pccf"), pattern = paste0("^pccf_", pr, "_\\d+\\.rds$"))
  k0 <- if (length(done)) max(as.integer(sub(".*_(\\d+)\\.rds$", "\\1", done))) else 0L
  chunks <- split(seq_len(nrow(ids)), ceiling(seq_len(nrow(ids)) / CHUNK_SIZE))
  for (k in seq_along(chunks)) {
    tag <- sprintf("pccf_%s_%02d", pr, k0 + k)
    saveRDS(ids[chunks[[k]], ID], file.path(OUT, "pccf", paste0(tag, "_ids.rds")))
    # seed is fixed per province and chunk number so the random WCF allocation is reproducible
    run_pccf_chunk(tag, ids[chunks[[k]]], seed = as.integer(pr) * 100L + k0 + k)
  }
}

# ---- Step 3: reference layers and linkage -----------------------------------

load_refs <- function() {
  msg("Loading reference layers")
  gaf <- fread(GAF_FILE, colClasses = "character",
               select = c("PRUID_PRIDU", "CSDUID_SDRIDU", "SACTYPE_CSSGENRE",
                          "SACCODE_CSSCODE", "CMAUID_RMRIDU", "DAUID_ADIDU"))
  setnames(gaf, c("PRUID", "CSDUID", "SACTYPE", "SACCODE", "CMAUID", "DAUID"))
  gaf_da <- unique(gaf, by = "DAUID")   # DAs nest within CSDs

  ior <- fread(IOR_FILE, colClasses = list(character = "CSDuid"))
  ior <- ior[, .(CSDUID = CSDuid, IOR = suppressWarnings(as.numeric(Index_of_remoteness)))]

  msdi <- as.data.table(read_excel(MSDI_XLSX, sheet = "Data"))
  msdi <- msdi[, .(DAUID = as.character(DA), SCOREMAT, SCORESOC, QuintMat, QuintSoc,
                   QuintMatCR, QuintSocCR, Centmat, CentSoc)]

  da <- st_read(DA_SHP, quiet = TRUE)[, c("DAUID", "PRUID")]
  if (!isTRUE(st_crs(da)$epsg == CRS_LAMBERT)) da <- st_transform(da, CRS_LAMBERT)

  list(gaf_da = gaf_da, ior = ior, msdi = msdi, da = da)
}

sac_labels <- c("1" = "CMA", "2" = "CA (tracted)", "3" = "CA (non-tracted)",
                "4" = "Strong MIZ", "5" = "Moderate MIZ", "6" = "Weak MIZ",
                "7" = "No MIZ", "8" = "Territories (outside CA)")

dmt_group <- function(dmt, dmtdiff) {
  d <- fifelse(dmt == "Z" & !is.na(dmtdiff) & dmtdiff != "", dmtdiff, dmt)
  fcase(d == "A",               "Urban letter carrier: house",
        d == "B",               "Urban letter carrier: apartment",
        d %chin% c("E", "G"),   "Business / large volume receiver",
        d %chin% c("H", "T"),   "Rural route / suburban service",
        d == "W",               "Rural postal code",
        d %chin% c("J", "K", "M", "R", "X"), "PO box / general delivery / other",
        default = "Unknown / retired")
}

link_province <- function(pr, refs) {
  f_out <- file.path(OUT, "linked", paste0("linked_", pr, ".rds"))
  if (!need(f_out)) return(invisible(NULL))
  msg("Link", pr)
  x <- readRDS(file.path(OUT, "nar", paste0("nar_ref_", pr, ".rds")))

  # PCCF+ output (one row per ID from the residential routine)
  pf <- list.files(file.path(OUT, "pccf"), pattern = paste0("^pccf_", pr, "_\\d+\\.rds$"),
                   full.names = TRUE)
  pp <- list.files(file.path(OUT, "pccf"), pattern = paste0("^pccf_", pr, "_\\d+_problem\\.rds$"),
                   full.names = TRUE)
  pc <- rbindlist(lapply(pf, function(f) as.data.table(readRDS(f))), fill = TRUE)
  prob_ids <- unique(rbindlist(lapply(pp, function(f) as.data.table(readRDS(f))[, .(ID, LINK)])),
                     by = "ID")
  pc <- unique(pc, by = "ID")
  setnames(pc, setdiff(names(pc), "ID"), paste0("pccf_", setdiff(names(pc), "ID")))
  pc[, pccf_PCODE := NULL]
  x <- merge(x, pc, by.x = "ADDR_GUID", by.y = "ID", all.x = TRUE)
  x[, pccf_problem := ADDR_GUID %chin% prob_ids$ID]

  # PCCF+ point to Lambert
  x[, `:=`(pccf_lat = suppressWarnings(as.numeric(pccf_LAT)),
           pccf_lon = suppressWarnings(as.numeric(pccf_LONG)))]
  x[, c("pccf_LAT", "pccf_LONG") := NULL]
  ok <- x[, !is.na(pccf_lat) & !is.na(pccf_lon)]
  xy <- st_coordinates(st_transform(
    st_as_sf(x[ok, .(pccf_lon, pccf_lat)], coords = c("pccf_lon", "pccf_lat"), crs = CRS_NAD83),
    CRS_LAMBERT))
  x[ok, `:=`(pccf_x = xy[, 1], pccf_y = xy[, 2])]
  x[, pccf_coded := ok & !is.na(pccf_DAUID) & pccf_DAUID != ""]

  # Reference point to 2021 DA (point in polygon; nearest DA for points just outside
  # the cartographic boundary, e.g. on shorelines)
  da_pr <- refs$da[refs$da$PRUID == pr, ]
  pts <- st_as_sf(x[, .(ref_x, ref_y)], coords = c("ref_x", "ref_y"), crs = CRS_LAMBERT)
  hit <- st_intersects(pts, da_pr)
  idx <- vapply(hit, function(h) if (length(h)) h[1] else NA_integer_, integer(1))
  x[, ref_da_method := fifelse(is.na(idx), "nearest", "within")]
  if (anyNA(idx)) idx[is.na(idx)] <- st_nearest_feature(pts[is.na(idx), ], da_pr)
  x[, ref_DAUID := da_pr$DAUID[idx]]
  rm(pts, hit, idx); gc()

  # Geography and rurality at the reference (true) location
  g <- refs$gaf_da[, .(ref_DAUID = DAUID, ref_CSDUID = CSDUID, ref_SACTYPE = SACTYPE,
                       ref_SACCODE = SACCODE, ref_CMAUID = CMAUID)]
  x <- merge(x, g, by = "ref_DAUID", all.x = TRUE)
  x <- merge(x, refs$ior[, .(ref_CSDUID = CSDUID, ref_IOR = IOR)], by = "ref_CSDUID", all.x = TRUE)
  x[, rurality_sac := factor(sac_labels[ref_SACTYPE], levels = sac_labels)]

  # MSDI at both DAs
  m <- refs$msdi
  x <- merge(x, setnames(copy(m), paste0("ref_", names(m))), by = "ref_DAUID", all.x = TRUE)
  x <- merge(x, setnames(copy(m), paste0("pccf_", names(m))), by = "pccf_DAUID", all.x = TRUE)

  # Delivery mode and agreement flags
  x[, delivery_mode := dmt_group(pccf_DMT, pccf_DMTDIFF)]
  x[, `:=`(same_da  = pccf_coded & pccf_DAUID == ref_DAUID,
           same_csd = pccf_coded & pccf_CSDuid == ref_CSDUID)]

  saveRDS(x, f_out)

  cnt <- fread(file.path(OUT, "nar", paste0("nar_counts_", pr, ".csv")),
               colClasses = list(character = "pr"))
  cnt[, `:=`(n_pccf_coded     = x[, sum(pccf_coded)],
             n_pccf_problem   = x[, sum(pccf_problem)],
             n_ref_da_nearest = x[, sum(ref_da_method == "nearest")],
             n_msdi_ref       = x[, sum(!is.na(ref_QuintMat))],
             n_msdi_pccf      = x[, sum(!is.na(pccf_QuintMat))],
             pct_same_da      = round(100 * x[pccf_coded == TRUE, mean(same_da)], 1))]
  fwrite(cnt, file.path(OUT, "linked", paste0("counts_", pr, ".csv")))
  msg("Link", pr, ": coded", format(cnt$n_pccf_coded, big.mark = ","), "of",
      format(nrow(x), big.mark = ","), "| same DA", cnt$pct_same_da, "%")
}

# ---- Step 4: NDVI check -----------------------------------------------------

check_ndvi <- function() {
  nc <- list.files(NDVI_DIR, pattern = "\\.(nc|nc4)$", full.names = TRUE, recursive = TRUE)
  if (!length(nc)) {
    msg("MODIS: no netCDF found in", NDVI_DIR)
    return(invisible(NULL))
  }
  for (f in nc) {
    msg("NDVI:", basename(f), sprintf("(%.0f MB)", file.size(f) / 1e6))
    if (requireNamespace("ncdf4", quietly = TRUE)) {
      h <- ncdf4::nc_open(f)
      msg("  variables:", paste(names(h$var), collapse = ", "))
      msg("  dimensions:", paste(sprintf("%s[%d]", names(h$dim),
                                         vapply(h$dim, `[[`, numeric(1), "len")), collapse = ", "))
      ncdf4::nc_close(h)
    }
  }
  if (!requireNamespace("terra", quietly = TRUE))
    msg("NDVI: install the 'terra' package before P2 extraction")
}

# ---- Step 5: combine --------------------------------------------------------

combine_all <- function() {
  # Always rebuilt from every finished province, so partial runs accumulate
  f_out <- file.path(OUT, "p1_linked_all.rds")
  lf <- list.files(file.path(OUT, "linked"), pattern = "^linked_\\d+\\.rds$", full.names = TRUE)
  msg("Combining", length(lf), "province file(s)")
  all <- rbindlist(lapply(lf, readRDS), fill = TRUE)
  saveRDS(all, f_out)
  cf <- list.files(file.path(OUT, "linked"), pattern = "^counts_\\d+\\.csv$", full.names = TRUE)
  cnt <- rbindlist(lapply(cf, fread, colClasses = list(character = "pr")), fill = TRUE)
  fwrite(cnt, file.path(OUT, "p1_counts.csv"))
  print(cnt)
  msg("Combined:", format(nrow(all), big.mark = ","), "addresses ->", f_out)
}

# ---- Run --------------------------------------------------------------------

t0 <- Sys.time()
if ("nar" %in% STEPS)  for (pr in PROVINCES) extract_nar(pr)
if ("pccf" %in% STEPS) { check_pccf_packages(); for (pr in PROVINCES) run_pccf(pr) }
if ("link" %in% STEPS) { refs <- load_refs(); for (pr in PROVINCES) link_province(pr, refs) }
if ("ndvi" %in% STEPS) check_ndvi()
if ("combine" %in% STEPS) combine_all()
msg("Done in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
