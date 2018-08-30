# utility function to print metadata for nexus file

print.nex <- function(x) {
  nchar <- ncol(x$data)
  ntax <- nrow(x$data)
  cat('NEXUS data file with', nchar, 'morphological characters and', ntax, 'taxa\n')
}


# TODO printneat() function for character pairs
# x = nexus object
# y = list of characters (numeric)
printneat <- function(x, y) {
  charnames <- x$charlabels[y]
  statenames <- x$statelabels[y]
  charnames <- str_replace_all(charnames, "Note\\:.*?$", "")
  charnames <- str_replace_all(charnames, "\\[.*?\\]", "")
  statenames <- str_replace_all(statenames, "\\[.*?\\]", "")
  statenames <- str_replace_all(statenames, "Note\\:.*?$", "")
  paste0(charnames, statenames)
}


# sort nexus file by something

sort.nex <- function(x, by = c('taxlabels', 'charlabels'), ...) {
  res <- x
  by <- match.arg(by)
  if (by == 'taxlabels') {
    ord <- order(x$taxlabels)
    res$data <- x$data[ord, ]
    res$taxlabels <- sort(x$taxlabels)
    return(res)  # it is important to have return() for some reason
  }
  if (by == 'charlabels') {
    ord <- order(x$charlabels, ...)
    res$data <- x$data[, ord]
    res$charlabels <- sort(x$charlabels, ...)
    return(res)
  }
}

# subset a nexus file

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


# remove NAs

na.omit.nex <- function(x) {
  res <- x
  x$data <- as.matrix(x$data)
  id.col <- sapply(1:ncol(x$data), function(z) { all(is.na(x$data[,z])) } )
  id.row <- sapply(1:nrow(x$data), function(z) { all(is.na(x$data[z,])) } )
  res$data <- x$data[!id.row, !id.col]
  res$charlabels <- x$charlabels[!id.col]
  res$statelabels <- x$statelabels[!id.col]
  res$taxlabels <- x$taxlabels[!id.row]
  res$charpartition <- x$charpartition[!id.col]
  res$charset <- x$charset[!id.col]
  res$charnums <- x$charnums[!id.col]
  res$file <- x$file[!id.col]
  res
}

head.nex <- function(x) {
  utils::head(x$data)
}



# subset a nexus file

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

# is a nexus file
is.nex <- function(x) {
  inherits(x, "nex")
}


# function to output a summary table of precision and recall
# true = unique pairs of true duplicated characters
# test = list of pairs of potential duplicated characters (each member of list is a different method)
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

