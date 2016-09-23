cleantext <- function(x, comma=TRUE) {
  # remove comments in square brackets
  x <- str_replace_all(x, "[\\s]*[\\[\\(].*?[\\]\\)][\\s]*", "")
  # make lowercase
  x <- sapply(x, tolower)
  # remove extra spaces
  x <- gsub("\\s{2,}", " ", x)
  # take out beginning spaces
  x <- gsub("^\\s", "", x)
  x <- gsub("\\:", " ", x)
  x <- gsub("\\s$", "", x)
  # x <- removePunctuation(x)  # might work better not doing this??
  x <- gsub("\\/", " ", x)
  if (comma) {
    x <- gsub("\\,", "", x)
  }
  x <- gsub("\\bthe\\b", "", x)
  x <- gsub("\\bin\\b", "", x)
  x <- gsub("\\bat\\b", "", x)
  x <- gsub("\\ba\\b", "", x)
  # take out common words (if, to, the...)
  # x <- removeWords(x, stopwords("english"))
  names(x) <- NULL
  x
}
