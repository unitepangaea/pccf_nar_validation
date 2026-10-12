# Patched replacement for random_select_by_pcode() in PCCF+ 8D's R engine
# (DATA_R/PCCFplus_8D_do_not_edit.R). The engine file itself is not modified:
# p1_data_assembly.R swaps this definition in when it runs each chunk.
#
# Two changes from the StatCan version:
#   1. Speed. The original filters the whole lookup table once per input record,
#      so run time grows with records x table size. This version splits the table
#      once per lookup and reuses it. It makes the same sample() call on the same
#      rows in the same order, so with the same seed the draws are identical.
#   2. Zero or missing weights. When every weight in a group is 0 or NA (e.g. a DA
#      whose blocks all have zero or missing 2021 population), the original computes
#      0/0 and sample() stops with "NA in probability vector". Here those groups get
#      equal-probability selection, NA weights elsewhere are treated as 0, and each
#      such record is counted in .rsp_state$n_fallback for reporting.

.rsp_state <- new.env()
.rsp_state$n_fallback <- 0L
.rsp_state$n_na_weight <- 0L

random_select_by_pcode <- function(pcode_value, data, PCODEn, with_prob = FALSE, wgt = NULL) {
  key_col <- rlang::as_name(rlang::enquo(PCODEn))
  sig <- paste(rlang::obj_address(data), key_col, nrow(data),
               paste(names(data), collapse = ","))
  if (!identical(.rsp_state$sig, sig)) {
    .rsp_state$sig <- sig
    .rsp_state$idx <- split(seq_len(nrow(data)), data[[key_col]])
  }

  # filter(key == value) drops NA keys, so an NA value matches nothing
  rows <- if (is.na(pcode_value)) NULL else .rsp_state$idx[[as.character(pcode_value)]]
  if (length(rows) == 0) {
    return(NULL)  # No matching rows found
  }
  filtered_rows <- data[rows, , drop = FALSE]

  probs <- NULL
  if (with_prob) {
    w <- filtered_rows[[rlang::as_name(rlang::enquo(wgt))]]
    if (anyNA(w)) {
      .rsp_state$n_na_weight <- .rsp_state$n_na_weight + 1L
      w[is.na(w)] <- 0
    }
    if (sum(w) > 0) {
      probs <- w / sum(w)
    } else {
      .rsp_state$n_fallback <- .rsp_state$n_fallback + 1L  # equal probability
    }
  }

  filtered_rows[sample(seq_len(nrow(filtered_rows)), size = 1, prob = probs), , drop = FALSE]
}
