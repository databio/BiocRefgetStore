#' Reverse complement a DNA sequence
#'
#' @param seq DNA sequence string
#' @return Reverse complement string
#' @keywords internal
.reverse_complement <- function(seq) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq)))
}

#' Null coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
