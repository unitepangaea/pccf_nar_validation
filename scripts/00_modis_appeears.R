# =============================================================================
# 00  MODIS MOD13Q1 v061 EVI for Canada through the NASA AppEEARS API
#
# Resubmits the area request used in the paper (data/modis_request.json: EVI,
# 250 m, 16-day composites, 1 January 2025 to 30 September 2026, native sinusoidal
# projection, netCDF-4; AppEEARS adds the VI Quality layer), waits for processing
# and downloads the result into data/raw/MODIS. Processing by AppEEARS can take
# several hours; the script can be stopped and rerun with the task id.
#
# Requires a free NASA Earthdata login (https://urs.earthdata.nasa.gov). Set
# EARTHDATA_USERNAME and EARTHDATA_PASSWORD in ~/.Renviron (never in this repository).
#
# Run:  Rscript scripts/00_modis_appeears.R             submit a new request and wait
#       APPEEARS_TASK=<task id> Rscript scripts/00_modis_appeears.R   resume a request
# =============================================================================

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))
suppressPackageStartupMessages({ library(httr2); library(jsonlite) })

API <- "https://appeears.earthdatacloud.nasa.gov/api"
OUT <- file.path(RAW, "MODIS"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
user <- Sys.getenv("EARTHDATA_USERNAME"); pass <- Sys.getenv("EARTHDATA_PASSWORD")
if (!nzchar(user) || !nzchar(pass)) stop("Set EARTHDATA_USERNAME and EARTHDATA_PASSWORD in ~/.Renviron")

token <- request(file.path(API, "login")) |> req_auth_basic(user, pass) |> req_method("POST") |>
  req_perform() |> resp_body_json() |> (\(x) x$token)()
auth <- function(req) req_auth_bearer_token(req, token)

task <- Sys.getenv("APPEEARS_TASK")
if (!nzchar(task)) {
  body <- fromJSON(file.path(ROOT, "data", "modis_request.json"), simplifyVector = FALSE)
  task <- request(file.path(API, "task")) |> auth() |> req_body_json(body) |> req_perform() |>
    resp_body_json() |> (\(x) x$task_id)()
  message("Submitted AppEEARS task ", task, " (rerun with APPEEARS_TASK=", task, " to resume)")
}

repeat {
  st <- request(file.path(API, "task", task)) |> auth() |> req_perform() |> resp_body_json()
  message(format(Sys.time(), "%H:%M"), " status: ", st$status)
  if (identical(st$status, "done")) break
  if (st$status %in% c("error", "expired")) stop("AppEEARS task ", task, " ended with status ", st$status)
  Sys.sleep(600)
}

files <- request(file.path(API, "bundle", task)) |> auth() |> req_perform() |> resp_body_json()
for (f in files$files) {
  dest <- file.path(OUT, basename(f$file_name))
  if (file.exists(dest) && file.size(dest) == f$file_size) next
  message("Downloading ", f$file_name, " (", round(f$file_size / 1e9, 2), " GB)")
  request(file.path(API, "bundle", task, f$file_id)) |> auth() |> req_timeout(36000) |>
    req_perform(path = dest)
}
nc <- list.files(OUT, "\\.nc$", full.names = TRUE)
message("MODIS netCDF: ", paste(basename(nc), collapse = ", "))
if (length(nc) == 1 && basename(nc) != "MOD13Q1.061_250m_aid0001.nc")
  message("Note: scripts/p2_variables.R expects MOD13Q1.061_250m_aid0001.nc; rename the file if needed.")
