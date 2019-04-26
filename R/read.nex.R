#' Read in nexus data
#' 
#' This function reads NEXUS data stored in nexus, text, or pdf files
#' NOTE: for pdf files, it is necessary to first install poppler with, e.g.,
#' `brew install poppler`
#' 
#' @param file (required) path to either a .nex, .txt, or .pdf file containin data
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
#' @return an object of class \code{nex} for use in further \code{phenotools} functions
#' 
#' @examples \dontrun{
#' # Build a `nex` object from text files:
#' charpath <- system.file("extdata", "brusatte2014_charlist.txt",
#' package = "phenotools")
#' matpath <- system.file("extdata", "brusatte2014_matrix.txt",
#' package = "phenotools")
#' x <- read.nex(matpath, charlabels=charpath)
#' plot(x, legend.pos = "top")
#' 
#' # Build a `nex` object from a PDF file (Bertelli & Chiappe 2005):
#' x <- read.nex(file = system.file("extdata", "Bertelli_2005.pdf",
#' package = "phenotools"), ntax = 34, nchar = 63, first = 20, last = 22)
#' x
#' plot(x)
#' 
#' # Read in a nexus file:
#' x <- read.nex(file = system.file("extdata", "clarke_2006.nex",
#' package = "phenotools"))
#' }
#' 
#' @export
#' 
#' @import stringr
#' @importFrom stats sd
#' @importFrom stats na.omit
#' 
#' @references Bertelli, S., & Chiappe, L. M. (2005). Earliest Tinamous (Aves: Palaeognathae)
#' from the Miocene of Argentina and Their Phylogenetic Position. Contributions
#' in Science, 502, 1–20.
#' 
#' @author Chad Eliason \email{celiason@@fieldmuseum.org}
#'
# old:
# read.nex <- function(file, filetype=c('nexus', 'txt', 'pdf'), charlabels=NULL,
read.nex <- function(file, charlabels=NULL, charnums=NULL, statelabels=NULL,
	taxlabels=NULL, filename=NULL, ntax=NULL, nchar=NULL, first=NULL, last=NULL,
	missing="?", gap="-") {

	# filetype <- match.arg(filetype, choices=c('nexus', 'txt', 'pdf'))
	
	# auto-identify input file type:
	# if (length(filetype)==3) {
		if (grepl("pdf$", file)) {
			filetype <- "pdf"
		} else if (grepl("txt$", file)) {
			filetype <- "txt"
		} else if (grepl("nex$", file)) {
			filetype <- "nexus"
		} else {
			stop("Unrecognized file type")
		}
	# }

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
		# all(diff(nchars)==0)
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

	if (filetype == "nexus") {
		x <- scan(file = file, what = "", sep = "\n")
		
		# find number of characters
		nchar <- as.numeric(na.omit(stringr::str_extract(x, regex('(?<=NCHAR=)\\d+', ignore_case=TRUE))))
		# throws an error if ntax specified multiple times
		ntax <- as.numeric(na.omit(stringr::str_extract(x, regex('(?<=NTAX=)\\d+', ignore_case=TRUE)))[1])
		
		# get taxon labels
		taxlabelsstart <- grep('TAXLABELS', x, ignore.case=TRUE) + 1
		# lines that had END + semicolon
		taxlabelsend <- taxlabelsstart + ntax - 1
		taxlabels <- stringr::str_match(x[taxlabelsstart:taxlabelsend], '[\\t]*(.*)')[,2]
		# remove space at start of taxon label
		taxlabels <- gsub('^\\s*', '', taxlabels)
		# remove puncutation at end of taxon label
		taxlabels <- gsub('(\\s|\\;)$', '', taxlabels)
		# this makes it possible to read mesquite saved nexus files with taxa in a single line separated by spaces
		if (length(taxlabels)!=ntax){
			taxlabels <- strsplit(taxlabels, split="\\s")[[1]]
		}
		
		# extract data matrix
		matstart <- grep('MATRIX$', x, ignore.case=TRUE) + 1
		ends <- grep('\\;', x)
		matend <- ends[which(ends > matstart)[1]] - 1
		# mesquite saves spaces between polymorphic characters (annoying)
		# convert to single string of text (causing problems downstream?)
		mat <- x[matstart:matend]
		mat <- gsub("^(\t|\\s)+", "", mat)
		mat <- paste0(mat, collapse="\n")
		mat <- paste0("\n", mat)
		taxlabels0 <- taxlabels
		taxlabels <- gsub("[^_'A-Za-z0-9]", " ", taxlabels)
		taxlabels <- gsub("\\s{2,}", " ", taxlabels)
		
		# replace "bad" characters in taxon names
		for (i in seq_along(taxlabels)) {
			mat <- stringr::str_replace_all(string=mat, pattern = stringr::fixed(taxlabels0[i]), replacement=taxlabels[i])
		}
		locs <- stringr::str_locate(mat, paste0("\n", taxlabels, "(\\b|\\t)"))
		
		# Two formats:
		# Genus_species
		# 'Genus species'
		# there's a problem if some OTUs are just genus and others genus + species (with spaces between)

		# break up matrix by start/end locations of taxon labels
		mat <- stringr::str_sub(mat, start = locs[,1], end = c(locs[2:nrow(locs), 1] - 1, stringr::str_length(mat)))

		# remove taxon labels from matrix
		mat <- stringr::str_replace_all(mat, stringr::fixed(taxlabels), "")

		# remove white space
		mat <- gsub("\\s", "", mat)

		# remove commas in data matrix
		mat <- gsub(',', '', mat)

		# this didn't work with theropod working matrix (Turner 2012 AMNH version)
		# so fixed with:
		mat <- stringr::str_extract_all(mat, '\\d{1}|[\\(\\[\\{]\\d{1,4}[\\)\\]\\}]|\\-|\\?')  # extract scorings

		# convert '{ }', '[ ]' --> '( )'
		mat <- lapply(mat, gsub, pattern='\\{|\\[', replacement='\\(')
		mat <- lapply(mat, gsub, pattern='\\}|\\]', replacement='\\)')

		taxlabels <- gsub(' ', '_', taxlabels)
		taxlabels <- gsub("'", "", taxlabels)
		taxlabels <- gsub('"', '', taxlabels)

		# checks
		if (length(mat)!=ntax) {
			warning('Number of rows in data matrix not equal to number of taxa.')
		}
		if (any(sapply(mat, length) != nchar)) {
			warning('Number of characters for some taxa does not equal to that in defined by `nchar`')
		}

		mat <- do.call(rbind, mat)
		symbols <- paste(unique(unlist(strsplit(gsub('[\\(\\)\\??\\-]', '', sort(unique(as.vector(mat)))), ""))),collapse="")
		symbols <- gsub('[\\?\\-]','',symbols)

		mat <- ifelse(mat==missing, NA, mat)
		
		res <- list(taxlabels = taxlabels, data = mat, symbols = symbols, gap = gap, missing = missing)

		# extract charpartitions
		# use this format for character partition by body region: CHARPARTITION bodyparts=head: 1-4 7, body:5 6, legs:8-10;
		if (length(grep('\\bCHARPARTITION\\b', x, ignore.case=TRUE)) > 0) {
			charstart <- grep('CHARPARTITION', x, ignore.case=TRUE)
			charparts <- lapply(charstart, function(i) {
			  charpart <- x[i]
				# charpartname <- str_match(charpart, '^\\w+')
				charmatch <- stringr::str_match_all(charpart, "(\\w+):([0-9\\s\\-]*)")[[1]]
				# charmatch <- str_match_all(charpart, '(\\w+):[\\s]*(\\d+[\\-\\s]*[\\d+]*)')[[1]]
				charpartsets <- charmatch[,2]
				charpartranges <- charmatch[,3]
				ids <- sapply(charpartranges, text2numeric)
				if (is.matrix(ids)) {
					ids <- ids[,1]
					names(ids) <- rep(charpartsets, length(ids))
				} else {
					names(ids) <- charpartsets
				}
				# now create a vector and set names at ids according to charpartsets labels
				charparts <- rep(NA, nchar)
				for (i in 1:length(ids)) {
					charparts[ids[[i]]] <- names(ids)[i]
				}
				# res$charpartition <- charparts
				charparts
			})
			# look for file partition
			if (any(grep('CHARPARTITION file', x, ignore.case=TRUE))) {
				id <- grep('CHARPARTITION file', x, ignore.case=TRUE)
				id1 <- match(id, charstart)
				id2 <- match(setdiff(charstart, id), charstart)
				res$file <- charparts[[id1]]
				res$charpartition <- charparts[[id2]]
			} else {
				res$file <- rep(stringr::str_extract(file, '\\w+[\\s\\w]*(?=\\.nex)'), ncol(mat))
				res$charpartition <- charparts[[1]]
			}
		} else {
			res$file <- rep(stringr::str_extract(file, '\\w+[\\s\\w]*(?=\\.nex)'), ncol(mat))
			res$charpartition <- rep("''", ncol(mat))
		}

		# get character names and numbers
		if (length(grep('\\bCHARLABELS', x, ignore.case=TRUE)) > 0) {
			charlabelsstart <- grep('\\bCHARLABELS', x, ignore.case=TRUE) + 1
			charlabelsend <- grep('\\;$', x[charlabelsstart:length(x)])[1] + charlabelsstart - 2	
			charnames <- na.omit(stringr::str_match(x[charlabelsstart:charlabelsend], '\\[(\\d{1,4}|[A-Z]+)[.]*\\]\\s(.+)'))
			charnames <- gsub("^'|'$", "", charnames)
			charnums <- as.numeric(charnames[,2])
			charlabels <- charnames[,3]
			res$charlabels <- charlabels
			res$charnums <- charnums
		} else {
			res$charlabels <- rep("''", ncol(mat))
			res$charnums <- 1:ncol(mat)
		}

		# get state labels
		if (length(grep('\\bSTATELABELS', x, ignore.case=TRUE)) > 0) {
			statelabelsstart <- grep('\\bSTATELABELS', x, ignore.case=TRUE) + 1
			statelabelsend <- grep('\\;$', x[statelabelsstart:length(x)])[1] + statelabelsstart - 2
			statelabels <- stringr::str_match(x[statelabelsstart:statelabelsend], '[\\t]*\\d+\\s*(.+)$')[,2]
			# sometimes states are on multiple lines:
			if (length(statelabels) != nchar) {
				statelabelsend <- grep('\\;$', x[statelabelsstart:length(x)])[1] + statelabelsstart - 1
				statelabels <- x[statelabelsstart:statelabelsend]
				statelabels <- paste(statelabels, collapse="")
				states_start <- stringr::str_locate_all(statelabels, "\\d+[\\t\\s]+[\\']")  # start
				states_end <- stringr::str_locate_all(statelabels, "\\t+\\,|\\t+\\;")  # end
				statelabels <- sapply(1:nchar, function(x) {substr(statelabels, start=states_start[[1]][x,1], stop=states_end[[1]][x,1])})
				statelabels <- gsub('\\t|^\\d+', '', statelabels)
			}
			statelabels <- gsub('\\,$|\\,(?=\\s\\[)', '', statelabels, perl=TRUE)  # cut out the commas
			if (length(statelabels) != nchar) {
				warning('state labels found not equal to number of characters')
			}
			res$statelabels <- statelabels
		} else {
			res$statelabels <- rep("''", ncol(mat))
		}

		if (length(grep('\\bCHARSTATELABELS', x, ignore.case=TRUE)) > 0) {
			charlabelsstart <- grep('\\bCHARSTATELABELS', x, ignore.case=TRUE) + 1
			charlabelsend <- grep('\\;', x[charlabelsstart:length(x)])[1] + charlabelsstart - 1
			# charstatelabels <- str_match_all(x[charlabelsstart:charlabelsend], "(\\t|\\,)\\s(\\d{1,})\\s'(.*?)'(\\s\\/\\s)?(.*?)(?=\\;|(\\,\\s\\d))")
			
			charstatelabels <- stringr::str_match_all(x[charlabelsstart:charlabelsend],
				"(\\t|\\,)\\s(\\d{1,})\\s'(.*?)('\\s\\/(.*?))?(?=\\;(\\s)?$|(\\,\\s\\d{1,}\\s'))")

			charnums <- as.numeric(charstatelabels[[1]][,3])
			charlabels <- charstatelabels[[1]][,4]
			statelabels <- charstatelabels[[1]][,6]
			res$charlabels <- charlabels
			res$charnums <- charnums
			res$statelabels <- gsub("^\\s|\\s$", "", statelabels)
		}
		
		res$statelabels <- gsub("''", "' '", res$statelabels)

		class(res) <- c('nex', 'list')

	}

	res

}
