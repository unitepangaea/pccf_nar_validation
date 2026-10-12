# =============================================================================
# Match NAR addresses to the Open Database of Addresses (ODA v1, 2021)
# Sourced by p1_data_assembly.R; builds the ODA tier of the reference-point
# hierarchy (NAR building point -> ODA point -> NAR blockface point).
#
# Match: civic number + street name + street type + direction, standardized the same
# way in both sources, searched across the province. A candidate ODA address is
# accepted when it is in the same CSD (code) OR within ODA_MAX_DIST metres of the
# NAR blockface point; among accepted candidates the one nearest the blockface point
# (or the only same-CSD candidate) is kept. Accuracy is validated against NAR
# building points in each province (data/derived/nar/ref_validation_<pr>.csv).
#   - NAR CSD codes are interim 2026 (Location file); ODA CSD codes are 2016. Where
#     municipalities were restructured (e.g. New Brunswick, 2023) codes differ, so
#     the blockface distance gate does the work there.
#   - ODA points sharing a key are averaged; keys whose points spread more than
#     ODA_MAX_SPREAD metres are not used.
#   - ODA does not state a datum; WGS84 is assumed (< 2 m from NAD83 in Canada).
# =============================================================================

ODA_DIR        <- file.path(RAW, "ODA")
ODA_MAX_SPREAD <- 50     # m: ODA points sharing an address must lie this close together
ODA_MAX_DIST   <- 1000   # m: distance gate between ODA point and NAR blockface point
CRS_ODA        <- 4326

oda_abbr <- c("10" = "NL", "11" = "PE", "12" = "NS", "13" = "NB", "24" = "QC", "35" = "ON",
              "46" = "MB", "47" = "SK", "48" = "AB", "59" = "BC", "60" = "YT", "61" = "NT",
              "62" = "NU")

# ---- Standardization (applied identically to NAR and ODA) --------------------

to_ascii_upper <- function(x) {
  x <- toupper(iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = ""))
  trimws(gsub("\\s+", " ", gsub("[^A-Z0-9 ]", " ", x)))
}

# Street name: drop English/French particles, standardize Saint/ordinals, remove spaces
key_name <- function(x) {
  x <- paste0(" ", to_ascii_upper(x), " ")
  x <- gsub(" ([0-9]+)(ST|ND|RD|TH|E|ER|IEME) ", " \\1 ", x)
  x <- gsub(" (DE|DU|DES|LA|LE|LES|D|L|THE|OF|AUX|AU)(?= )", "", x, perl = TRUE)
  x <- gsub(" SAINTE ", " STE ", x)
  x <- gsub(" SAINT ", " ST ", x)
  x <- gsub(" MOUNT ", " MT ", x)
  gsub(" ", "", x, fixed = TRUE)
}

# Street type: one code for English and French spellings and abbreviations
type_map <- c(
  STREET = "ST", ST = "ST", RUE = "RUE", R = "RUE",
  ROAD = "RD", RD = "RD", AVENUE = "AV", AVE = "AV", AV = "AV",
  DRIVE = "DR", DR = "DR", CRESCENT = "CRES", CRES = "CRES", CROISSANT = "CRES", CROIS = "CRES",
  COURT = "CRT", CRT = "CRT", CT = "CRT", BOULEVARD = "BLVD", BLVD = "BLVD", BOUL = "BLVD", BD = "BLVD",
  PLACE = "PL", PL = "PL", CIRCLE = "CIR", CIR = "CIR", TERRACE = "TERR", TERR = "TERR",
  TERRASSE = "TERR", TSSE = "TERR", SQUARE = "SQ", SQ = "SQ", CARRE = "CAR", CAR = "CAR",
  PARKWAY = "PKY", PKY = "PKY", PKWY = "PKY", SIDEROAD = "SIDERD", SIDERD = "SIDERD",
  HIGHWAY = "HWY", HWY = "HWY", GARDENS = "GDNS", GDNS = "GDNS", HEIGHTS = "HTS", HTS = "HTS",
  CONCESSION = "CONC", CONC = "CONC", CIRCUIT = "CIRCT", CIRCT = "CIRCT",
  TOWNLINE = "TLINE", TLINE = "TLINE", PRIVATE = "PVT", PVT = "PVT", PARK = "PK", PK = "PK",
  POINT = "PT", PT = "PT", TRAIL = "TRAIL", TRL = "TRAIL", LANE = "LANE", LN = "LANE",
  CHEMIN = "CH", CH = "CH", RANG = "RANG", RG = "RANG", MONTEE = "MONTEE", MTEE = "MONTEE",
  ALLEE = "ALLEE", ALL = "ALLEE", IMPASSE = "IMP", IMP = "IMP", ROUTE = "RTE", RTE = "RTE",
  COTE = "COTE", PROMENADE = "PROM", PROM = "PROM", PASSAGE = "PASS", PASS = "PASS",
  RUELLE = "RLE", RLE = "RLE", SENTIER = "SENT", SENT = "SENT", GROVE = "GROVE",
  GATE = "GATE", MEWS = "MEWS", PATH = "PATH", WALK = "WALK", WAY = "WAY", LINE = "LINE",
  CLOSE = "CLOSE", COMMON = "COMMON", ROW = "ROW", RUN = "RUN", RIDGE = "RIDGE",
  ISLAND = "ISLAND", QUAY = "QUAY", HILL = "HILL", VIEW = "VIEW", LANDING = "LANDNG",
  LANDNG = "LANDNG", ESTATES = "ESTATE", ESTATE = "ESTATE", GREEN = "GREEN", HEATH = "HEATH",
  BAY = "BAY", BEND = "BEND", COVE = "COVE", GLEN = "GLEN", MANOR = "MANOR", ACRES = "ACRES")
key_type <- function(x) {
  x <- to_ascii_upper(x)
  fifelse(x %chin% names(type_map), unname(type_map[x]), x)
}

dir_map <- c(NORTH = "N", NORD = "N", N = "N", SOUTH = "S", SUD = "S", S = "S",
             EAST = "E", EST = "E", E = "E", WEST = "W", OUEST = "W", O = "W", W = "W",
             NORTHEAST = "NE", NE = "NE", NORTHWEST = "NW", NW = "NW", NO = "NW",
             SOUTHEAST = "SE", SE = "SE", SOUTHWEST = "SW", SW = "SW", SO = "SW")
key_dir <- function(x) {
  x <- to_ascii_upper(x)
  fifelse(x %chin% names(dir_map), unname(dir_map[x]), x)
}

key_civic <- function(x) {
  x <- toupper(gsub("\\s+", "", x))
  sub("^0+(?=[0-9])", "", x, perl = TRUE)
}

ODA_KEYS  <- c("k_csd", "k_civ", "k_name", "k_type", "k_dir")
ADDR_KEYS <- c("k_civ", "k_name", "k_type", "k_dir")

# ---- ODA for one province: one point per address key ---------------------------

read_oda <- function(pr) {
  f <- file.path(ODA_DIR, sprintf("ODA_%s_v1", oda_abbr[pr]), sprintf("ODA_%s_v1.csv", oda_abbr[pr]))
  if (!file.exists(f)) return(NULL)
  o <- fread(f, colClasses = "character", encoding = "UTF-8",
             select = c("latitude", "longitude", "street_no", "str_name", "str_name_pcs",
                        "str_type", "str_type_pcs", "str_dir", "str_dir_pcs", "csduid"))
  o <- o[street_no != "" & (str_name != "" | str_name_pcs != "")]
  # Provider-supplied fields first; ODA's processed (_pcs) fields only where blank,
  # because the processing truncates some names (e.g. "Woods" -> "WOOD")
  o[, `:=`(k_csd  = csduid,
           k_civ  = key_civic(street_no),
           k_name = key_name(fifelse(str_name != "", str_name, str_name_pcs)),
           k_type = key_type(fifelse(str_type != "", str_type, str_type_pcs)),
           k_dir  = key_dir(fifelse(str_dir != "", str_dir, str_dir_pcs)))]
  xy <- st_coordinates(st_transform(
    st_as_sf(o[, .(lon = as.numeric(longitude), lat = as.numeric(latitude))],
             coords = c("lon", "lat"), crs = CRS_ODA), CRS_LAMBERT))
  o[, `:=`(x = xy[, 1], y = xy[, 2])]
  o[, .(oda_x = mean(x), oda_y = mean(y), oda_n = .N,
        spread = sqrt(max((x - mean(x))^2 + (y - mean(y))^2))), by = ODA_KEYS
  ][spread <= ODA_MAX_SPREAD][, spread := NULL][]
}

# ---- Add oda_x / oda_y / oda_n to NAR records ------------------------------------
# x needs CSD_CODE, CIVIC_NO, CIVIC_NO_SUFFIX, OFFICIAL_STREET_NAME/TYPE/DIR

match_oda <- function(x, pr) {
  oda <- read_oda(pr)
  if (is.null(oda)) {
    msg("ODA", pr, ": no ODA file for", oda_abbr[pr])
    return(x[, `:=`(oda_x = NA_real_, oda_y = NA_real_, oda_n = NA_integer_, oda_gate = NA_character_)][])
  }
  k <- x[, .(ADDR_GUID, nar_csd = CSD_CODE, bf_x, bf_y,
             k_civ  = key_civic(paste0(CIVIC_NO, CIVIC_NO_SUFFIX)),
             k_name = key_name(OFFICIAL_STREET_NAME),
             k_type = key_type(OFFICIAL_STREET_TYPE),
             k_dir  = key_dir(OFFICIAL_STREET_DIR))][k_civ != "" & k_name != ""]
  cand <- merge(k, oda, by = ADDR_KEYS, allow.cartesian = TRUE)
  cand[, `:=`(same_csd = !is.na(nar_csd) & nar_csd == k_csd,
              d_bf = sqrt((oda_x - bf_x)^2 + (oda_y - bf_y)^2))]
  cand <- cand[same_csd | (!is.na(d_bf) & d_bf <= ODA_MAX_DIST)]
  # nearest to blockface; without a blockface point, accept only a single same-CSD candidate
  cand[, n_cand := .N, by = ADDR_GUID]
  cand <- cand[!is.na(d_bf) | n_cand == 1L]
  setorder(cand, ADDR_GUID, d_bf, na.last = TRUE)
  best <- cand[!duplicated(ADDR_GUID)]
  best[, oda_gate := fifelse(same_csd, "same CSD", "blockface distance")]
  msg("ODA", pr, ":", format(nrow(best), big.mark = ","), "NAR addresses matched (",
      best[oda_gate == "same CSD", .N], "same CSD,", best[oda_gate != "same CSD", .N], "by distance gate )")
  merge(x, best[, .(ADDR_GUID, oda_x, oda_y, oda_n, oda_gate)], by = "ADDR_GUID", all.x = TRUE)
}

# ---- Accuracy of ODA and blockface points against NAR building points ------------

ref_validation <- function(x, pr) {
  v <- x[!is.na(bg_x)]
  d_oda <- v[!is.na(oda_x), sqrt((oda_x - bg_x)^2 + (oda_y - bg_y)^2)]
  d_bf  <- v[!is.na(bf_x),  sqrt((bf_x - bg_x)^2 + (bf_y - bg_y)^2)]
  acc <- function(d, n_base) data.table(
    N = length(d), pct_of_building_records = round(100 * length(d) / n_base, 1),
    median_m = round(median(d)), p90_m = round(quantile(d, 0.9)),
    within_50m = round(100 * mean(d <= 50), 1), within_250m = round(100 * mean(d <= 250), 1))
  rbind(cbind(pr = pr, point = "ODA", if (length(d_oda)) acc(d_oda, nrow(v)) else data.table(N = 0L)),
        cbind(pr = pr, point = "NAR blockface", acc(d_bf, nrow(v))), fill = TRUE)
}
