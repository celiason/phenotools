#' Function to remove spaces, punctuation from phenome terms
#' comma = whether 
#' fast = whether to use fast version (for duplicated analysis) or slow version (for printing characters)
#' schinke = whether to use Latin token algorithm from Schinke 1996 (?)
#' example:
#' cleantext("dorsalmost part of the processus ligamentus cranii")
cleantext <- function(x, comma=TRUE, fast=TRUE, latin=TRUE, cuts=TRUE, comments=TRUE) {
  removeSpaces <- function(x) {
    x <- gsub("\\s{2,}", " ", x)  # remove extra spaces
    x <- gsub("^\\s", "", x)  # take out beginning spaces
    x <- gsub("\\s$", "", x)  # remove end whitespace
    x
  }
  # remove comments in square brackets
  if (comments) {
    x <- str_replace_all(x, "[\\s]*[\\[\\(].*?[\\]\\)][\\s]*", " ")
    x <- str_replace_all(x, "Note\\:.*?$", " ")
    x <- str_replace_all(x, "\\([Ff]ig(\\.|ure).*\\)", " ")
  }
  # if (comma) {
  #   x <- gsub("\\,", "", x)
  # }
  tocut1 <- c("with", "than", "then", "those", "with", "to", "the", "and", "an", "a", "or", "of", "for", "not", "along", "length", "less", "below", "above", "around", "longer", "shorter", "sits", "absent", "present", "form", "process", "state", "view", "margin", "shape", "placed", "recess", "(un)?ordered", "external")
  tocut2 <- c("along", "less", "below", "above", "around", "longer", "shorter", "sits", "absent", "present")
  tocut3 <- c("lateral", "distal", "ventral", "posterior", "anterior", "medial", "dors.*?", "external")
  tocut <- c(tocut1, tocut2, tocut3)
  # tocut <- paste0("\\b", tocut, "(\\w*)?", collapse="|")
  # tocut <- paste0("\\b", tocut, "\\b", collapse="|")
  if (fast) {
    x <- Corpus(VectorSource(x))
    x <- tm_map(x, content_transformer(tolower))  # Convert the text to lower case
    x <- tm_map(x, removeNumbers)  # Remove numbers
    x <- tm_map(x, removeWords, stopwords("english"))  # Remove english common stopwords
    if (cuts) {
      x <- tm_map(x, removeWords, tocut)  # specify stopwords
    }
    x <- tm_map(x, gsub, pattern="[[:punct:]]", replacement=" ")
    if (latin) {
      x <- tm_map(x, schinke)  # Latinized tokens
    }
    x <- tm_map(x, removeSpaces)  # Eliminate extra white spaces
    # return(lapply(x, as.character))
    return(x)
  }
  x <- sapply(x, strsplit, " ")
  x <- lapply(x, tolower) 
  x <- lapply(x, removeNumbers)
  x <- lapply(x, removeWords, stopwords("english"))
  if (cuts) {
    x <- lapply(x, removeWords, tocut)
  }
  x <- lapply(x, gsub, pattern="[[:punct:]]", replacement=" ")
  if (latin) {
    x <- lapply(x, schinke)
  }
  x <- lapply(x, removeSpaces)
  return(lapply(x, as.character))
}
