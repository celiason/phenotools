# x <- part1
# fast = TRUE
cleanchars <- function(x, fast=TRUE) {
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
		x <- tm_map(x, removeWords, tocut)  # specify stopwords
		# x <- tm_map(x, removePunctuation, preserve_intra_word_contractions = TRUE, preserve_intra_word_dashes = TRUE, ucp=TRUE)  # Remove punctuations
		x <- tm_map(x, gsub, pattern="[[:punct:]]", replacement=" ")
		x <- tm_map(x, stripWhitespace)  # Eliminate extra white spaces
		x <- tm_map(x, schinke)  # Latinized tokens
		return(lapply(x, as.character))
	}
	x <- sapply(x, strsplit, " ")
	x <- lapply(x, tolower)	
	x <- lapply(x, removeNumbers)
	x <- lapply(x, removeWords, stopwords("english"))
	x <- lapply(x, removeWords, tocut)
	x <- lapply(x, gsub, pattern="[[:punct:]]", replacement=" ")
	x <- lapply(x, gsub, pattern="\\s+", replacement="")
	x <- lapply(x, schinke)
	return(lapply(x, as.character))
}
