#' Function to output a summary table of precision and recall
#' 
#' @param true a character vector denoting pairs of overlapping/duplicated
#' characters (e.g., "1--2")
#' @param test a list of pairs of potential duplicated characters (e.g., each
#' element of the list might represent a different method of locating duplicates
#' with the `duplicated.nex` function).
#' 
#' @return a data frame object with false negatives (FN), true positives (TP),
#' false positive (FP), along with precision (TP/(TP+FP)), recall (TP/(TP+FN)),
#' and F1 score (2 * (precision * recall) / (precision + recall))
#' 
dupsum <- function(true, test = list()) {
  fneg <- sapply(seq_along(test), function(x) { sum(!true %in% test[[x]]) })
  tpos <- sapply(seq_along(test), function(x) { sum(test[[x]] %in% true) })
  fpos <- sapply(seq_along(test), function(x) { sum(!test[[x]] %in% true) })
  precision <- tpos / (tpos + fpos)
  recall <- tpos / (tpos + fneg)
  f1 <- 2 * (precision * recall / (precision + recall))
  if (is.null(names(test))) {
    names(test) <- paste0("method", seq_along(test))
  }
  res <- data.frame(method = names(test), fneg, tpos, fpos, precision, recall, F1=f1)
  res
}
