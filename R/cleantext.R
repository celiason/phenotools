#' Function to remove spaces, punctuation from phenome terms
#' 
#' Given a character vector, remove comments, superfluous characters, extra spaces
#' and non-anatomical terms
#' 
#' @param x a character vector
#' @param fast whether to use fast version (for duplicated analysis) or slow version (for printing characters)
#' @param latin whether to use Latin token algorithm (see Schinke 1996)
#' @param cuts whether to cut terms from a list
#' @param comments logical whether to remove test in square brackets
#' 
#' @examples \dontrun{
#' cleantext("dorsalmost part of the processus ligamentus cranii")
#' }
#' 
#' @importFrom stringr str_replace_all
#' @importFrom tm Corpus
#' @importFrom tm VectorSource
#' @importFrom tm tm_map
#' @importFrom tm content_transformer
#' @importFrom tm removeNumbers
#' @importFrom tm removeWords
#' @importFrom tm stopwords
#' 
#' @export
#' 
#' @author Chad M. Eliason
#' 
cleantext <- function(x, fast=TRUE, latin=TRUE, cuts=TRUE, comments=TRUE) {
  require(tm)
  if (comments) {
    x <- str_replace_all(x, "[\\s]*[\\[\\(].*?[\\]\\)][\\s]*", " ")  # remove comments in square brackets
    x <- str_replace_all(x, "Note\\:.*?$", " ")  # remove note comments
    x <- str_replace_all(x, "\\([Ff]ig(\\.|ure).*\\)", " ")  # remove figure references
  }
  tocut1 <- c("with", "than", "then", "those", "with", "to", "the", "and", "an", "a", "or", "of", "for", "not", "along", "length", "less", "below", "above", "around", "longer", "shorter", "sits", "absent", "present", "form", "process", "state", "view", "margin", "shape", "placed", "recess", "(un)?ordered", "external")
  tocut2 <- c("along", "less", "below", "above", "around", "longer", "shorter", "sits", "absent", "present")
  tocut3 <- c("lateral", "distal", "ventral", "posterior", "anterior", "medial", "dors.*?", "external")
  tocut <- c(tocut1, tocut2, tocut3)
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


removeSpaces <- function(x) {
    x <- gsub("\\s{2,}", " ", x)  # remove extra spaces
    x <- gsub("^\\s", "", x)  # take out beginning spaces
    x <- gsub("\\s$", "", x)  # remove end whitespace
    x
}

# stemmer function based on Schinke et al. 1996
# TODO make it so things like -ity, -ed will be removed from end of word
schinke <- function(x) {
  x <- tolower(x)
  x <- gsub("m\\.", "muscul", x)
  x <- gsub("proc\\.", "processus", x)
  x <- gsub("(lig|ligg)\\.", "ligamentos", x)
  x <- gsub("n\\.", "nervos", x)
  # x <- gsub("j", "i", x)
  # x <- gsub("v", "u", x)
  # x <- str_replace_all(x, "(ibus|ius|ae|am|as|em|es|ia|is|nt|os|ud|um|us|a|e|i|o|u)\\b", "")
    # NEW VERSIon  (not exactly Schinke)
  x <- str_replace_all(x, "(ity|ed|al|ibus|ius|ae|am|as|em|es|ia|is|nt|os|ud|um|us|a|e|i|o|u)\\b", "")
  x
}

