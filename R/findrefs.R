# function to find instances of citations to other datasets/authors within character/state names
# x = `nex` object
findrefs <- function(x) {
  statelabels <- x$statelabels
  tmp <- str_match_all(statelabels, '([A-Z][a-z]+\\s(?:et\\sal\\.\\s|(?:\\&|and)\\s[A-Z][a-z]+\\s)?)[\\(]?(\\d{4})[\\)]?.{1,5}[Ch]ar(?:\\.|acter)\\s(\\d+)')
  names(tmp) <- seq_along(statelabels)
  # convert to matrix with only non-missing cases
  dd <- lapply(seq_along(tmp), function(z) {cbind(df_name = z, tmp[[z]])})
  res <- do.call('rbind', dd)
  res <- setNames(data.frame(res), c('charnum', 'text', 'author', 'year', 'charmatch'))
  res$id <- as.numeric(as.character(res$charnum))
  res$infile <- x$file[res$id]
  res$charnum <- x$charnum[res$id]
  res$matchfile <- paste0(tolower(str_extract(res$author, "[A-Z][a-z]+")), "_", res$year)
  res <- res[, c("infile", "charnum", "author", "year", "charmatch", "matchfile")]
  res
}

# new version to work with just text (character names, state names, etc.) and
# what filename/dataset the text correponds to
findrefs2 <- function(x, filenames) {
  # this is looking for the pattern: AUTHOR(s) YEAR CHARACTER CHARNUM
  tmp <- str_match_all(x, '([A-Z][a-z]+\\s(?:et\\sal\\.\\s|(?:\\&|and)\\s[A-Z][a-z]+\\s)?)[\\(]?(\\d{4})[\\)]?.{1,5}[Ch]ar(?:\\.|acter)\\s(\\d+)')
  names(tmp) <- seq_along(x)
  # convert to matrix with only non-missing cases
  dd <- lapply(seq_along(tmp), function(z) {cbind(df_name = z, tmp[[z]])})
  res <- do.call('rbind', dd)
  res <- setNames(data.frame(res), c('charnum', 'text', 'author', 'year', 'charmatch'))
  res$id <- as.numeric(as.character(res$charnum))
  res$infile <- filenames[res$id]
  res$charnum <- res$id
  res$matchfile <- paste0(tolower(str_extract(res$author, "[A-Z][a-z]+")), "_", res$year)
  res <- res[, c("infile", "charnum", "author", "year", "charmatch", "matchfile")]
  res
}