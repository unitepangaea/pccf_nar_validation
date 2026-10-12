# Licensed data: PCCF+ Version 8D

The Postal Code Conversion File Plus (PCCF+) is distributed by Statistics Canada under
a licence agreement (for example, through the Data Liberation Initiative) and **must not
be uploaded to this or any other public repository**. Nothing in this folder except this
README is tracked by git.

To run the pipeline, place your own licensed copy of PCCF+ Version 8D (December 2025
postal codes, 2021 census geography) in `data/licensed/PCCF8D/`, so that these exist:

```
data/licensed/PCCF8D/PCCFplus_8D.R
data/licensed/PCCF8D/DATA_R/PCCFplus_8D_do_not_edit.R
data/licensed/PCCF8D/DATA_R/...        (the R data files from DATA.zip)
```

To keep the files elsewhere, set `PCCF_DIR` in `config_local.R` (see
`config_local_example.R`).

`R/pccf_engine_patch.R` contains two run-time changes to one helper function of the
PCCF+ R engine (a faster candidate lookup that gives identical output, and equal-probability
allocation where every weight is zero or missing). The distributed PCCF+ files are not
modified.

## What is published from PCCF+ processing

Only stratum-level aggregates (error measures and classification probabilities by
province, SAC type and delivery mode, validation summaries and the paper's tables) are
published in `outputs/` and `correction/tables/`. No address-, postal-code- or DA-level
file containing PCCF+ fields (PCCF+ DA, delivery mode, SAC type or coordinates) is
published. The postal-code tables in `correction/tables/` are built from the National
Address Register alone.
