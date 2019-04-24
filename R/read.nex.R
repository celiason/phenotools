#' Read a nexus data file
#'
#' A function that reads data stored in nexus files
#'
#' @param file (required) path to nexus file
#' @param missing character representing missing data
#' @param gap character representing inapplicable/incomporable data
#' 
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' 
#' @examples \dontrun{
#' x <- read.nex(file = system.file("extdata", "clarke_2006.nex", package = "phenotools"))
#' x
#' plot(x)
#' }
#' 
#' @author Chad Eliason \email{celiason@@fieldmuseum.org}
#' 
#' @export
#'
read.nex <- function(file, missing = '?', gap = '-') {

  	x <- scan(file = file, what = "", sep = "\n")
	
	# find number of characters
	nchar <- as.numeric(na.omit(stringr::str_extract(x, regex('(?<=NCHAR=)\\d+', ignore_case=TRUE))))
	
	# throws an error if ntax specified multiple times
	ntax <- as.numeric(na.omit(stringr::str_extract(x, regex('(?<=NTAX=)\\d+', ignore_case=TRUE)))[1])
	
	# find taxon labels
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
	
	res
}


# Function to go from this "2,5-7,10,12-15" to this "c(2,5,6,7,10,12,13,14,15)"
# see http://r.789695.n4.nabble.com/convert-delimited-strings-with-ranges-to-numeric-td4673763.html
text2numeric <- function(xx) {
	xx <- gsub("^\\s*|\\s*$", "", xx)
  xx <- gsub('\\s|,\\s', ',', xx)
	xx <- gsub('\\-', ':', xx)
	eval(parse(text = paste("c(", xx, ")")))
}
