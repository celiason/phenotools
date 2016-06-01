cleantext <- function(x, comma=TRUE) {
  x <- str_replace_all(x, "[\\s]*[\\[\\(].+[\\]\\)][\\s]*", "")  # remove comments in square brackets
  x <- sapply(x, tolower)  # make lowercase
  x <- gsub("\\s{2,}", " ", x)  # remove extra spaces
  x <- gsub("^\\s", "", x)  # take out beginning spaces
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
  # x <- removeWords(x, stopwords("english"))  # take out common words (if, to, the...)
  names(x) <- NULL
  x
  }