# =============================================================================
# P2 Error and exposure variables
# PCCF+ geocoding validation paper (preprint target 2026-10-23)
#
# Input: data/derived/linked/linked_<pr>.rds from P1 (one row per residential NAR address)
# Output (data/derived/p2/):
#   p2_analysis.rds     one row per address: positional error, DA/CSD agreement,
#                       EVI at the reference and PCCF+ points, MSDI agreement, flags
#   p2_exclusions.csv   records available / excluded for each outcome, by reason
#   evi_points.rds      EVI extraction cache (point-level, both years)
#
# Variables
#   dist_m          geodesic distance, reference point -> PCCF+ point (haversine,
#                   mean Earth radius; < 0.5% from ellipsoidal distance)
#   same_da/csd     PCCF+ DA/CSD equals the DA/CSD of the reference point (from P1)
#   evi_<pt>_<w>_<yr>  growing-season mean EVI at pt = ref | pccf, window w =
#                   pix (pixel), b1 (3x3, ~700 m), b2 (5x5, ~1.2 km), yr = 2025 | 2026 | mean
#   evi_diff_<w>    PCCF+ minus reference, 2025-26 mean
#   msdi_*          quintile difference (PCCF+ minus reference) and misclassification
#
# EVI: MODIS MOD13Q1 v061 (AppEEARS netCDF, native sinusoidal 231.66 m grid).
#   Growing season = 16-day composites starting 1 Jun - 31 Aug. Observations kept when
#   MODLAND QA is 00/01 and no mixed clouds, possible snow/ice or possible shadow.
#   Pixel composite = mean of kept observations. The grid's northern edge is 64.4 N;
#   points beyond it, or with no kept observations, get NA (flagged).
#
# Run: Rscript scripts/p2_variables.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ncdf4)
})

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets RAW, PCCF_DIR and DERIVED
P1_DIR  <- Sys.getenv("P2_IN",  file.path(DERIVED, "linked"))  # env overrides for testing
OUT     <- Sys.getenv("P2_OUT", file.path(DERIVED, "p2"))
NC_FILE <- file.path(RAW, "MODIS", "MOD13Q1.061_250m_aid0001.nc")

EVI_VAR  <- "_250m_16_days_EVI"
QA_VAR   <- "_250m_16_days_VI_Quality"
GS_START <- "06-01"           # growing season: composites starting 1 Jun ...
GS_END   <- "08-31"           # ... to 31 Aug
YEARS    <- c(2025, 2026)
BUFFERS  <- c(b1 = 1L, b2 = 2L)   # half-widths in cells: 3x3 and 5x5 windows
STRIP    <- 633L              # rows per read (= netCDF chunk height)
FORCE_EVI <- FALSE            # TRUE re-extracts EVI even if evi_points.rds exists

CRS_LAMBERT <- 3347
CRS_SINU <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
R_EARTH  <- 6371008.8

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
t0 <- Sys.time()

# ---- 1. Assemble analysis columns from P1 -----------------------------------

keep <- c("ADDR_GUID", "PROV_CODE", "ref_source", "multi_unit", "delivery_mode", "pccf_DMT",
          "pccf_Link_Source", "pccf_RPF", "pccf_SLI", "pccf_coded", "pccf_problem",
          "rurality_sac", "ref_SACTYPE", "ref_IOR", "ref_DAUID", "pccf_DAUID",
          "ref_CSDUID", "pccf_CSDuid", "same_da", "same_csd",
          "ref_x", "ref_y", "pccf_x", "pccf_y", "pccf_lat", "pccf_lon",
          paste0(rep(c("ref_", "pccf_"), each = 6),
                 c("SCOREMAT", "SCORESOC", "QuintMat", "QuintSoc", "QuintMatCR", "QuintSocCR")))

lf <- list.files(P1_DIR, "^linked_\\d+\\.rds$", full.names = TRUE)
msg("Reading", length(lf), "P1 province files")
x <- rbindlist(lapply(lf, function(f) { d <- readRDS(f); d[, intersect(keep, names(d)), with = FALSE] }),
               fill = TRUE)
msg(format(nrow(x), big.mark = ","), "addresses")

# ---- 2. Positional error -----------------------------------------------------

ll <- st_coordinates(st_transform(st_as_sf(x[, .(ref_x, ref_y)], coords = c("ref_x", "ref_y"),
                                           crs = CRS_LAMBERT), 4269))
x[, `:=`(ref_lon = ll[, 1], ref_lat = ll[, 2])]
rm(ll)
haversine <- function(lat1, lon1, lat2, lon2) {
  p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 + cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * R_EARTH * asin(pmin(1, sqrt(a)))
}
x[, dist_m := fifelse(pccf_coded, haversine(ref_lat, ref_lon, pccf_lat, pccf_lon), NA_real_)]

# ---- 3. EVI extraction -------------------------------------------------------

nc <- nc_open(NC_FILE)
gx <- nc$dim$xdim$vals; gy <- nc$dim$ydim$vals          # cell centres
res <- abs(gx[2] - gx[1]); nx <- length(gx); ny <- length(gy)
dates <- as.Date(as.POSIXct("2000-01-01", tz = "UTC") + nc$dim$time$vals)

qa_ok <- function(q) {
  q <- as.integer(q)
  !is.na(q) & bitwAnd(q, 3L) <= 1L & bitwAnd(q, 1024L) == 0L &
    bitwAnd(q, 16384L) == 0L & bitwAnd(q, 32768L) == 0L
}

# grid row/col for a set of Lambert points
to_cell <- function(px, py) {
  ok <- !is.na(px) & !is.na(py)
  xy <- matrix(NA_real_, length(px), 2)
  xy[ok, ] <- st_coordinates(st_transform(
    st_as_sf(data.frame(x = px[ok], y = py[ok]), coords = c("x", "y"), crs = CRS_LAMBERT), CRS_SINU))
  col <- as.integer(round((xy[, 1] - gx[1]) / res)) + 1L
  row <- as.integer(round((gy[1] - xy[, 2]) / res)) + 1L
  out <- !is.na(col) & (col < 1L | col > nx | row < 1L | row > ny)
  col[out] <- NA_integer_; row[out] <- NA_integer_
  list(col = col, row = row)
}

extract_evi <- function(pts) {
  # pts: data.table(pid, col, row) for unique point cells; returns pix and buffer means per year
  res_all <- pts[, .(pid)]
  for (yr in YEARS) {
    slices <- which(dates >= as.Date(paste0(yr, "-", GS_START)) &
                    dates <= as.Date(paste0(yr, "-", GS_END)))
    msg("EVI", yr, ": composites", paste(format(dates[slices], "%b %d"), collapse = ", "))
    vals <- matrix(NA_real_, nrow(pts), 1L + length(BUFFERS),
                   dimnames = list(NULL, c("pix", names(BUFFERS))))
    halo <- max(BUFFERS)
    for (r0 in seq(1L, ny, by = STRIP)) {
      r1 <- min(ny, r0 + STRIP - 1L)
      in_strip <- which(!is.na(pts$row) & pts$row >= r0 & pts$row <= r1)
      if (!length(in_strip)) next
      ra <- max(1L, r0 - halo); rb <- min(ny, r1 + halo); nr <- rb - ra + 1L
      s <- matrix(0, nx, nr); n <- matrix(0L, nx, nr)
      for (t in slices) {
        e <- ncvar_get(nc, EVI_VAR, start = c(1L, ra, t), count = c(nx, nr, 1L))
        q <- ncvar_get(nc, QA_VAR,  start = c(1L, ra, t), count = c(nx, nr, 1L))
        ok <- qa_ok(q) & !is.na(e)
        s[ok] <- s[ok] + e[ok]; n[ok] <- n[ok] + 1L
      }
      comp <- s / n                      # NaN where no kept observation
      comp[n == 0L] <- NA_real_
      rm(s, n, e, q, ok)
      cc <- pts$col[in_strip]; rr <- pts$row[in_strip] - ra + 1L
      vals[in_strip, "pix"] <- comp[cbind(cc, rr)]
      for (b in names(BUFFERS)) {
        w <- BUFFERS[[b]]; sm <- numeric(length(cc)); ct <- integer(length(cc))
        for (dx in -w:w) for (dy in -w:w) {
          c2 <- cc + dx; r2 <- rr + dy
          inb <- c2 >= 1L & c2 <= nx & r2 >= 1L & r2 <= nr
          v <- rep(NA_real_, length(cc)); v[inb] <- comp[cbind(c2[inb], r2[inb])]
          g <- !is.na(v); sm[g] <- sm[g] + v[g]; ct[g] <- ct[g] + 1L
        }
        vals[in_strip, b] <- fifelse(ct > 0L, sm / ct, NA_real_)
      }
      rm(comp); gc(verbose = FALSE)
      if ((r0 - 1L) %/% STRIP %% 5L == 0L) msg("  rows", r0, "-", r1, "of", ny)
    }
    colnames(vals) <- paste0(colnames(vals), "_", yr)
    res_all <- cbind(res_all, as.data.table(vals))
  }
  res_all
}

f_evi <- file.path(OUT, "evi_points.rds")
if (FORCE_EVI || !file.exists(f_evi)) {
  rc <- to_cell(x$ref_x, x$ref_y)
  pc <- to_cell(x$pccf_x, x$pccf_y)
  x[, `:=`(ref_cell = fifelse(is.na(rc$col), NA_real_, (rc$row - 1) * nx + rc$col),
           pccf_cell = fifelse(is.na(pc$col), NA_real_, (pc$row - 1) * nx + pc$col))]
  cells <- unique(c(x$ref_cell, x$pccf_cell)); cells <- cells[!is.na(cells)]
  pts <- data.table(pid = cells, row = as.integer((cells - 1) %/% nx) + 1L,
                    col = as.integer((cells - 1) %% nx) + 1L)
  msg("EVI: extracting at", format(nrow(pts), big.mark = ","), "unique MODIS cells")
  evi <- extract_evi(pts)
  saveRDS(evi, f_evi)
  x[, c("ref_cell", "pccf_cell") := NULL]
} else {
  msg("EVI: using cached", f_evi)
  evi <- readRDS(f_evi)
}
nc_close(nc)

# attach to both points
rc <- to_cell(x$ref_x, x$ref_y); pc <- to_cell(x$pccf_x, x$pccf_y)
x[, `:=`(ref_cell = fifelse(is.na(rc$col), NA_real_, (rc$row - 1) * nx + rc$col),
         pccf_cell = fifelse(is.na(pc$col), NA_real_, (pc$row - 1) * nx + pc$col),
         ref_in_modis = !is.na(rc$col), pccf_in_modis = !is.na(pc$col))]
rm(rc, pc)
ev <- copy(evi); setnames(ev, c("ref_cell", paste0("evi_ref_", names(evi)[-1])))
x <- merge(x, ev, by = "ref_cell", all.x = TRUE)
ev <- copy(evi); setnames(ev, c("pccf_cell", paste0("evi_pccf_", names(evi)[-1])))
x <- merge(x, ev, by = "pccf_cell", all.x = TRUE)
rm(ev)
for (pt in c("ref", "pccf")) for (w in c("pix", names(BUFFERS))) {
  cols <- paste0("evi_", pt, "_", w, "_", YEARS)
  x[, (paste0("evi_", pt, "_", w, "_mean")) := rowMeans(.SD, na.rm = TRUE), .SDcols = cols]
  x[is.nan(get(paste0("evi_", pt, "_", w, "_mean"))), (paste0("evi_", pt, "_", w, "_mean")) := NA_real_]
}
for (w in c("pix", names(BUFFERS)))
  x[, (paste0("evi_diff_", w)) := get(paste0("evi_pccf_", w, "_mean")) - get(paste0("evi_ref_", w, "_mean"))]

# ---- 4. MSDI agreement ---------------------------------------------------------

x[, `:=`(msdi_mat_qdiff = pccf_QuintMat - ref_QuintMat,
         msdi_soc_qdiff = pccf_QuintSoc - ref_QuintSoc,
         msdi_mat_scorediff = pccf_SCOREMAT - ref_SCOREMAT,
         msdi_soc_scorediff = pccf_SCORESOC - ref_SCORESOC)]
x[, `:=`(msdi_mat_miscl = msdi_mat_qdiff != 0, msdi_soc_miscl = msdi_soc_qdiff != 0)]

# ---- 5. Exclusion audit --------------------------------------------------------

aud <- rbindlist(list(
  x[, .(outcome = "positional error", reason = fcase(!pccf_coded, "PCCF+ not coded",
                                                     default = "included"))],
  x[, .(outcome = "DA agreement", reason = fcase(!pccf_coded, "PCCF+ not coded",
                                                 default = "included"))],
  x[, .(outcome = "EVI (2025-26 mean, pixel)",
        reason = fcase(!pccf_coded, "PCCF+ not coded",
                       !ref_in_modis, "reference point outside MODIS grid",
                       !pccf_in_modis, "PCCF+ point outside MODIS grid",
                       is.na(evi_ref_pix_mean), "no valid EVI at reference",
                       is.na(evi_pccf_pix_mean), "no valid EVI at PCCF+",
                       default = "included"))],
  x[, .(outcome = "MSDI material quintile",
        reason = fcase(!pccf_coded, "PCCF+ not coded",
                       is.na(ref_QuintMat), "no MSDI at reference DA",
                       is.na(pccf_QuintMat), "no MSDI at PCCF+ DA",
                       default = "included"))],
  x[, .(outcome = "MSDI social quintile",
        reason = fcase(!pccf_coded, "PCCF+ not coded",
                       is.na(ref_QuintSoc), "no MSDI at reference DA",
                       is.na(pccf_QuintSoc), "no MSDI at PCCF+ DA",
                       default = "included"))]))
aud <- aud[, .N, by = .(outcome, reason)][, pct := round(100 * N / sum(N), 2), by = outcome]
setorder(aud, outcome, -N)
fwrite(aud, file.path(OUT, "p2_exclusions.csv"))
msg("Exclusion audit:"); print(aud)

x[, c("ref_cell", "pccf_cell") := NULL]
saveRDS(x, file.path(OUT, "p2_analysis.rds"))
msg("Saved", format(nrow(x), big.mark = ","), "rows ->", file.path(OUT, "p2_analysis.rds"))
msg("Done in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
