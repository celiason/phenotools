#' function to output a summary table of precision and recall
#' 
#' @param true unique pairs of true duplicated characters
#' @param test list of pairs of potential duplicated characters (each member of list is a different method)
#' 
dupsum <- function(true, test = list()) {
  fneg <- sapply(seq_along(test), function(x) { sum(!true %in% test[[x]]) })
  tpos <- sapply(seq_along(test), function(x) { sum(test[[x]] %in% true) })
  fpos <- sapply(seq_along(test), function(x) { sum(!test[[x]] %in% true) })
  precision <- round(tpos / (tpos + fpos), 6)
  recall <- round(tpos / (tpos + fneg), 6)
  if (is.null(names(test))) {
    names(test) <- paste0("method", seq_along(test))
  }
  res <- data.frame(method = names(test), fneg, tpos, fpos, precision, recall)
  res
}

