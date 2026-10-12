# Installs the R packages used by the pipeline (developed with R 4.6.1).
# The correction tool alone (R/geoerror.R) needs only data.table.
pkgs <- c(
  "data.table", "readxl", "writexl", "sf", "ncdf4", "ggplot2", "patchwork", "scales",  # pipeline
  "dplyr", "readr", "stringr", "gtsummary", "gt",                                       # PCCF+ R engine
  "httr2", "jsonlite")                                                                    # MODIS download
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
message("All packages installed: ", all(pkgs %in% rownames(installed.packages())))
