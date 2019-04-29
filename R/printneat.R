#' Print character names and states
#' 
#' @param x nexus object
#' @param y numeric vectors of character numbers to print
#' 
#' @export
#' 
printneat <- function(x, y) {
  charnames <- x$charlabels[y]
  statenames <- x$statelabels[y]
  charnames <- str_replace_all(charnames, "Note\\:.*?$", "")
  charnames <- str_replace_all(charnames, "\\[.*?\\]", "")
  statenames <- str_replace_all(statenames, "\\[.*?\\]", "")
  statenames <- str_replace_all(statenames, "Note\\:.*?$", "")
  paste0(charnames, statenames)
}
