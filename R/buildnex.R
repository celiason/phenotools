#' Read a nexus data file
#' 
#' A function that reads data stored in nexus files
#' NOTE: it is necessary to first install poppler with, e.g., `brew install poppler`
#' 
#' @param file (required) path to either a PDF file or a TXT file
#' @param charlabels optional path to text file with character statements
#' @param ntax number of taxa (needed for PDF file reading)
#' @param nchar number of characters (needed for PDF file reading)
#' @param first first page for data matrix in a PDF file
#' @param last last page for data matrix in a PDF file
#' @param missing character representing missing data
#' @param gap character representing inapplicable/incomporable data
#' @param charnums optional vector with character numbers
#' @param statelabels optional character vector with state labels for characters
#' @param taxlabels optional vector with names of taxa
#' @param filename vector with names of files associated with characters
#' 
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' 
#' @examples \dontrun{
#' # Build a nex object from raw text files:
#' charpath <- system.file("extdata", "brusatte2014_charlist.txt", package = "phenotools")
#' matpath <- system.file("extdata", "brusatte2014_matrix.txt", package = "phenotools")
#' x <- buildnex(matpath, charlabels=charpath)
#' plot(x, legend.pos = "top")
#' }
#' 
#' @export
#' 
#' @import stringr
#' @importFrom stats sd
#' @importFrom stats na.omit
#' 
#' @author Chad Eliason \email{chad_eliason@@utexas.edu}
#'
buildnex <- function(file, charlabels=NULL, charnums=NULL, statelabels=NULL, taxlabels=NULL, filename=NULL, ntax=NULL, nchar=NULL, first=NULL, last=NULL, missing="?", gap="-") {

	# identify input file type
	if (grepl("pdf$", file)) {
		filetype <- "pdf"
	} else if (grepl("txt$", file)) {
		filetype <- "txt"
	} else {
		stop("Unrecognized file type")
	}

	# text file input
	if (filetype=="txt") {
		data <- text2mat(file)
		if (is.null(taxlabels)) {
			if (is.null(rownames(data))) {
				taxlabels <- rep('', nrow(data))
			} else {
			taxlabels <- rownames(data)
			}
		}
		if (is.null(statelabels)) {
			statelabels <- rep('', ncol(data))
		}
		if (!is.null(charlabels)) {
			charlabels <- text2charlabels(charlabels)
		} else {
			charlabels <- rep('', ncol(data))
		}
		if (is.null(charnums)) {
			charnums <- 1:ncol(data)
		}
		if (is.null(filename)) {
			filename <- "file"
		}
		file <- rep(filename, ncol(data))
		dimnames(data) <- NULL
		res <- list(data = data, file = file, taxlabels = taxlabels, charlabels = charlabels, statelabels = statelabels, charnums = charnums, missing = "?", gap="-")
		class(res) <- c("nex", "list")
		res
	}

	# pdf input

	if (filetype=="pdf") {

		# convert pdf to text file
		syscall <- paste("pdftotext -layout -f ", first, " -l ", last, " '", file, "'", sep="")

		# might need to give option for this in case it doesn't read well
		# turning off `-layout` can help
		# syscall <- paste("pdftotext -f ", first, " -l ", last, " '", file, "'", sep="")

		# removing layout was a better option for Livezey and Zusi (2006)

		system(syscall)

		# load and scan newly created text file

		txtfile <- gsub('pdf', 'txt', file)

		raw <- scan(txtfile, what="", sep="\n")

		# remove whitespace at beginning of lines
		raw <- gsub("^\\s+", "", raw)

		# remove/clean up exported text files
		system(paste("rm '", txtfile, "'", sep=""))

		# find first character label
		# charlabelstart <- grep('^[\\s]*\\s1\\s', x) 
		# charlabelend <- grep(paste0('^[\\s]*', nchar, '\\s'), x[charlabelstart:length(x)]) + charlabelstart
		# charlabels <- x[charlabelstart:charlabelend]
		# charlabels <- paste0(charlabels, collapse="\n")
		# charlabels <- gsub('[^\\.]\\n\\s', '', charlabels)
		# charlabels <- gsub('^\\s', '\n', charlabels, perl=TRUE)
		# charmatches <- str_match_all(charlabels, regex('\\n[\\s]*(\\d{1,3})(.+)'))
		# charnums <- charmatches[[1]][,2]
		# charlabels <- charmatches[[1]][,3]

		####### START OF THIS NEW STUFF ############

		raw2 <- do.call(paste0, list(raw, collapse="\n"))

		# how to find the data matrix?

		# patterns to search for:

		# Genus species 00110101011001000010010002003020 (continuous string of data)
		tmp <- str_match_all(raw, "([A-Z][a-z]+\\s[a-z]+)[\\s\\t]*([0-9\\-\\?]+.+)$")

		tmp <- str_match_all(raw2, "([A-Za-z]+[\\._\\s]*[a-z]*[\\.]*)\\s+([0-9\\?\\-]{1}(?![a-z]).+)")

		# SEARCH PATTERN:
		# letters/numbers/periods/single spaces separated by tab/multiple spaces than
		# contiguous blocks of numbers/brackets/dashes/?
		tmp <- str_match_all(raw, "(^[A-Z].*?)\\s{2,}([0-9\\?\\-]{1}(?![a-z]).+)")

		# Bertelli et al. (2014) issue
		# problem is that there is a new line after the taxon name (Crypturellus undulatus)
		# want to be able to "eat" through text in front of this until reaching a pattern, like 00011010 100100220 2003000

		# Fix??
		# maybe - if subsequent lines have all text (no numbers) remove one line??

		# str_match_all(raw2, regex("\\n([A-Z][a-z]+\\s[a-z]+).*?([0-9\\-\\?]+.+)", multiline=F, dotall=F))



		# Genus a b a a b (spaces between alphanumeric data)
		# tmp <- str_match_all(raw, "([A-Z][a-z]+)[\\s\\t]+([a-z0-9\\-\\?^,]+.+)$")

		tmp <- do.call(rbind, tmp)

		# Genus 0 1 1 1 0 1 0 0 1 (spaces between numeric data)

		# Genus 011101001 01111210 (chunks of data greater than 3 long, separated by spaces)

		# really want to find all consecutive numbers/sep. by up to one space, and with/without
		# period after words, at start of line, proceeded by numbers:

		# in progress..
		# tmp <- str_match_all(raw2, regex("([A-Za-z\\.]*).*?(([0-9\\?\\[\\]\\-]{3,}\\s*)+)", dotall=TRUE, multiline=FALSE))
		# tmp <- str_match_all(raw2, regex("([A-Z].*?(?=\\s[A-Z])).*?(([0-9\\?\\[\\]\\-]{3,}\\s*)+)", dotall=TRUE))

		# TODO need to be able to look in a few following lines to see if data are there

		# or maybe look for data first, and then find label in preceding lines?

		# 011010101Taxon

		# Genus species
		# 1101012201010101 (numeric data on new line compared to taxon label)

		# this works! don't change!
		# tmp <- str_match_all(raw2, regex("([A-Z][a-z]+[\\.]?\\s[a-z]+[\\.]?).*?(([0-9\\?\\[\\]\\-]{3,}\\s*)+)", dotall=TRUE, multiline=FALSE))

		# not sure what this one does
		# tmp <- str_match_all(raw, "([A-Z][a-z]+\\s[a-z]+[\\.]?)[\\s\\t]+(([0-9\\-\\?\\[\\]]{3,}[\\s\\t\\z]*)+)")

		# extract taxon labels
		taxlabels <- tmp[, 2]

		# get data matrix
		datamatrix <- tmp[, 3]

		names(datamatrix) <- taxlabels

		# taxlabels <- lapply(tmp, "[", 4)

		# datamatrix <- lapply(tmp, "[", 3)

		# unique taxon labels
		untaxlabels <- unique(na.omit(unlist(taxlabels)))

		# nums <- table(unlist(taxlabels))  # these should all be equal
		# find species with less data
		# if (dim(table(nums)) > 1) {
		# 	names(nums)[[which(nums %in% which.min(table(nums)))]]
		# }

		# remove spaces in data
		datamatrix <- gsub("\\s", "", datamatrix)

		# merge data for same species in multiple rows/lines of the text
		datamatrix <- sapply(untaxlabels, function(x) {paste0(datamatrix[which(names(datamatrix) %in% x)], collapse="")})

		# extract character states
		datamatrix <- str_extract_all(datamatrix, '\\d{1}|[\\(\\[\\{]\\d{1,4}[\\)\\]\\}]|\\-|\\–|\\?|\\w')  # extract scorings
		nchars <- sapply(datamatrix, length)

		# check - all same number of characters?
		all(diff(nchars)==0)

		# same number of characters as specified in input?
		# nchars==nchar

		names(datamatrix) <- taxlabels

		# only keep correct number of characters
		datamatrix <- datamatrix[nchars==nchar]


		# get first 1:nchar characters
		datamatrix <- sapply(seq_along(datamatrix), function(x) {datamatrix[[x]][1:nchar]})

		colnames(datamatrix) <- untaxlabels[nchars==nchar]

		# remove cases of NAs (probably not real characters)
		# datamatrix <- t(datamatrix[, !apply(datamatrix, 2, anyNA)])

		# checks
		if (! ntax == nrow(datamatrix) & nchar == ncol(datamatrix) ) {
			warning("Number of rows in data matrix does not equal number of input taxa")
		}

		# rownames(datamatrix) <- gsub(" ", "_", rownames(datamatrix))

		# dim(datamatrix)

		datamatrix <- t(datamatrix)

		# res <- list(data = datamatrix, taxlabels = taxlabels, missing = missing, gap = gap, symbols = symbols)
		res <- list(data = datamatrix, taxlabels = taxlabels, missing = missing, gap = gap)

		class(res) <- c('nex', 'list')

		res

	}

	res

}



# convert string of scorings to a matrix
text2mat <- function(x) {
	x <- readLines(x)
	require(stringr)
	if (any(str_detect(x, "[Mm]atrix"))) {
		x <- x[!grepl("[Mm]atrix", x)]
	}
	# scorings have to be separated from taxon names by either a tab or 2+ spaces (as in Xu et al. 2011)
	x <- str_match(x, "(.*?)(\t|\\s)(.*)")
	# scorings have to be separated from taxon names by either a tab or 2+ spaces (as in Brusatte et al. 2014)
	# x <- str_match(x, "(.*?)(\t|\\s{2,})(.*)")
	nms <- x[, 2]
	scores <- x[, 4]
	res <- str_extract_all(scores, '\\d{1}|[\\(\\[\\{]\\d{1,4}[\\)\\]\\}]|\\-|\\?')
	names(res) <- nms
	if (!sd(sapply(res, length))==0) {
		stop("There is a problem with reading the text. Not an equal number of characters for each taxon.")
	}
	res <- do.call("rbind", res)
	res
}

# get character labels from text
text2charlabels <- function(x) {
	x <- readLines(x)
	# matches <- str_match(x, "^(Character|[\\s\\t]*)(\\d+)[\\.]?(.*)(\\[\\[\\]\\])(.*)")
	# DONE: make it so this selects all text on multiple lines
	# TODO: make it so we can extract characters and character states separately (e.g., separated by a ":")
	# for Brusatte 2014
	if (any(str_detect(x, "^[Cc]haracter"))) {
		x <- paste0(x, collapse="\n")
		# matches <- str_match_all(x, regex("Character\\s*(\\d*)\\:\\s(.*?)", dotall=TRUE))
		matches <- str_match_all(x, regex("Character\\s*(\\d*)\\:\\s(.+?)(?=Character\\s*\\d*\\:\\s|$)", dotall=TRUE))
		matches <- matches[[1]]
		charnums <- as.numeric(matches[, 2])
		# id <- which(!is.na(as.numeric(matches[, 2])))
		# charnums <- na.omit(as.numeric(matches[, 2]))
		# seqs <- sapply(seq_along(id), function(i) {
		# 	if (id[i] == max(id) ) {
		# 		id[i]:length(x)
		# 	} else {
		# 		id[i]:(id[i+1]-1)
		# 	}
		# })
		# charlabs <- sapply(seq_along(seqs), function(i) {
		# 	paste0(x[seqs[[i]]], collapse=" ")
		# })
		charlabs <- matches[, 3]
	} else {
		# for Xu 2011
		matches <- str_match(x, regex("(Character\\s*)?(\\d*)[\\.]?[\\:]?[\\s]?(.*?)$", dotall=TRUE))
		charnums <- as.numeric(matches[, 3])
		charlabs <- matches[, 4]
	}
	res <- stats::setNames(charlabs, charnums)
	res <- res[!is.na(res)]
	res <- gsub("\n", " ", res)
	res
}

# setwd("/Users/chadeliason/github/nexustools/")
# buildnex(file = 'example/Bertelli_Chiappe_2005.pdf', ntax=34, nchar=63, first=22, last=22)
# buildnex(file = "example/brusatte2014_matrix.txt", charlabels="example/brusatte2014_charlist.txt")
