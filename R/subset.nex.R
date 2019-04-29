#' Subset a nexus file
#' 
#' @param x nexus file
#' @param condition subset condition
#' @param ... other arguments (not used)
#' 
#' @export
#' 
subset.nex <- function(x, condition) {
  res <- x
  cond <- substitute(condition)
  ss <- eval(substitute(condition), envir=x)
  if (any(grep('charpartition|charlabels|charset|file', cond))) {
    res$charlabels <- x$charlabels[ss & !is.na(ss)]
    res$charpartition <- x$charpartition[ss & !is.na(ss)]
    res$charset <- x$charset[ss & !is.na(ss)]
    res$charnums <- x$charnums[ss & !is.na(ss)]
    res$statelabels <- x$statelabels[ss & !is.na(ss)]
    res$file <- x$file[ss & !is.na(ss)]
    res$data <- x$data[, ss & !is.na(ss)]
    return(res)
  }
  if (any(grep('taxlabels', cond))) {
    res$taxlabels <- x$taxlabels[ss & !is.na(ss)]
    res$data <- x$data[ss & !is.na(ss), ]
    return(res)
  }
}

#' @param i rows (taxa) to keep
#' @param j columns (characters) to keep
#' @param drop unused
#' 
#' @rdname subset.nex
#' @export
#' 
'[.nex' <- function(x, i, j=NULL, drop = FALSE, ...) {
  res <- x
  if (is.null(j)) {
    j <- seq_along(x$charlabels)
  }
  res$charlabels <- x$charlabels[j]
  res$charpartition <- x$charpartition[j]
  res$charset <- x$charset[j]
  res$charnums <- x$charnums[j]
  res$statelabels <- x$statelabels[j]
  res$file <- x$file[j]
  res$data <- x$data[i, j, drop = drop]
  res$taxlabels <- x$taxlabels[i]
  # print(as.matrix(x$data[i, j, drop = drop, ...]))
  return(res)
  # rownames(res) <- x$taxlabels[i]
  # colnames(res) <- sapply(x$charlabels[j], substring, first=1, last=8)
}
