# Public input data

Filled by `scripts/00_download_data.R` and `scripts/00_modis_appeears.R`. Only the Index
of Remoteness (`IOR-2021/`, Statistics Canada Open Licence) is stored in the repository.

| Folder | Source | Licence | How to obtain |
|---|---|---|---|
| `NAR-202606/` | Statistics Canada, National Address Register, June 2026 (46-26-0002) | Statistics Canada Open Licence | automatic |
| `ODA/` | Statistics Canada, Open Database of Addresses v1 (46-26-0001) | Municipal and provincial open-data licences, listed with attribution statements in each `Data_Sources.csv` | automatic |
| `IOR-2021/` | Statistics Canada, Index of Remoteness 2021 (17-26-0001) | Statistics Canada Open Licence | included |
| `GAF-2021/` | Statistics Canada, Geographic Attribute File 2021 (92-151-X) | Statistics Canada Open Licence | manual (link printed by the script) |
| `lda_000b21a_e/` | Statistics Canada, 2021 dissemination area cartographic boundary file | Statistics Canada Open Licence | manual |
| `Canada2021_en/` | INSPQ, Material and Social Deprivation Index 2021, Canada | INSPQ terms of use | manual |
| `MODIS/` | NASA LP DAAC, MOD13Q1 v061 via AppEEARS | NASA open data policy | `scripts/00_modis_appeears.R` (Earthdata login) |

Statistics Canada Open Licence: https://www.statcan.gc.ca/en/terms-conditions/open-licence
