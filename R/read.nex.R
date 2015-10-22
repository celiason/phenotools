#' Read a nexus data file
#'
#' a function that reads data stored in nexus files
#'
#' @param file (required) path to nexus file
#' @param missing character representing missing data
#' @param gap character representing inapplicable/incomporable data
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' @examples \dontrun{
#'
#' x <- read.nex(file='example/toy1.nex')
#'
#' @author Chad Eliason \email{chad_eliason@@utexas.edu}
#'
#'
read.nex <- function(file, missing = '?', gap = '-') {

  	x <- scan(file = file, what = "", sep = "\n")
	
	# find number of characters
	nchar <- as.numeric(na.omit(str_extract(x, regex('(?<=NCHAR=)\\d+', ignore_case=TRUE))))
	
	# problem if ntax specified multiple times
	ntax <- as.numeric(na.omit(str_extract(x, regex('(?<=NTAX=)\\d+', ignore_case=TRUE)))[1])
	
	# where are taxon labels
	taxlabelsstart <- grep('TAXLABELS', x, ignore.case=TRUE) + 1
	taxlabelsend <- grep('\\;', x[taxlabelsstart:length(x)])[1] + taxlabelsstart - 2
	# taxlabels <- str_match(x[taxlabelsstart:taxlabelsend], '[\\t]*(\\w+)')[,2]
	taxlabels <- str_match(x[taxlabelsstart:taxlabelsend], '[\\t]*(.+)')[,2]
	taxlabels <- gsub('^ ', '', taxlabels)
	# this makes it possible to read mesquite saved nexus files with taxa in a single line separated by spaces
	if (length(taxlabels)!=ntax){
		taxlabels <- strsplit(taxlabels, split="\\s")[[1]]
	}
	# extract data matrix
	matstart <- grep('MATRIX$', x, ignore.case=TRUE) + 1
	matend <- grep('\\;', x[matstart:length(x)])[1] + matstart - 2
	# start edit
	# mesquite saves spaces between polymorphic characters (annoying)
	mat <- paste0(x[matstart:matend], collapse="")

	# add word boundaries for taxa without 
	if (is.na(any(str_locate(pattern="\\s", taxlabels)))) {
		locs <- str_locate(mat, paste0("\\b", taxlabels, "\\b", sep=""))
	} else {
		locs <- str_locate(mat, taxlabels)
	}

	mat <- str_sub(mat, start=locs[,1], end=c(locs[2:nrow(locs),1]-1, str_length(mat)))

	mat <- str_replace_all(mat, taxlabels, "")

	mat <- gsub("\\s", "", mat)

	# stop edit

	# mat <- str_replace_all(x[matstart:matend], taxlabels, "")  # remove taxon labels
	mat <- gsub(',', '', mat)  # remove commas in data matrix
	mat <- str_extract_all(mat, '\\d{1}|[\\(\\[\\{]\\d{1,4}[\\)\\]\\}]|\\-|\\?|\\w')  # extract scorings
	mat <- lapply(mat, gsub, pattern='\\{|\\[', replacement='\\(')
	mat <- lapply(mat, gsub, pattern='\\}|\\]', replacement='\\)')

	taxlabels <- gsub(' ', '_', taxlabels)
	taxlabels <- gsub("'", "", taxlabels)
	# taxlabels <- gsub('"', '', taxlabels)


	# checks
	if (length(mat)!=ntax)
		warning('Number of rows in data matrix not equal to number of taxa.')
	if (all(sapply(seq_along(mat), function(x) {length(mat[[x]]==nchar)})!=nchar))
		warning('Number of characters for some taxa does not equal to that in defined by `nchar`')
	# list(taxlabels=taxlabels, nchar=nchar, ntax=ntax, data=setNames(mat, taxlabels))
	mat <- do.call(rbind, mat)
	symbols <- paste(unique(unlist(strsplit(gsub('[\\(\\)\\??\\-]', '', sort(unique(as.vector(mat)))), ""))),collapse="")
	
	# symbols <- paste(sort(unique(as.vector(mat))), collapse="")
	symbols <- gsub('[\\?\\-]','',symbols)

	mat <- ifelse(mat==missing, NA, mat)

	# mat <- as.data.frame(mat)  # Maybe remove this??? slows things down a bit
	
	# charset <- rep(file, ncol(mat))
	# charset <- rep(str_extract(file, '\\w+[\\s\\w]*\\.nex'), ncol(mat))

	file <- rep(str_extract(file, '\\w+[\\s\\w]*\\.nex'), ncol(mat))

	res <- list(taxlabels = taxlabels, data = mat, symbols = symbols, gap = gap, missing = missing, file = file)

	# extract charpartitions
	# use this format for character partition by body region:
	# CHARPARTITION bodyparts=head: 1-4 7, body:5 6, legs:8-10;
	if (length(grep('CHARPARTITION', x, ignore.case=TRUE)) > 0) {
		charstart <- grep('CHARPARTITION', x, ignore.case=TRUE) + 1
		charpart <- na.omit(str_match(x, regex('CHARPARTITION\\s(.+)\\;', ignore_case=TRUE)))[,2]
		charpartname <- str_match(charpart, '^\\w+')
		charmatch <- str_match_all(charpart, '(\\w+):[\\s]*(\\d+[\\-\\s]*[\\d+]*)')[[1]]
		charpartsets <- charmatch[,2]
		charpartranges <- charmatch[,3]
		# I got this here - http://r.789695.n4.nabble.com/convert-delimited-strings-with-ranges-to-numeric-td4673763.html
		text2numeric <- function(xx) {
			xx <- gsub('\\s|,\\s', ',', xx)
			xx <- gsub('\\-', ':', xx)
			eval(parse(text = paste("c(", xx, ")")))
		}
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
		res$charpartition <- charparts
	} else {
		res$charpartition <- rep("''", ncol(mat))
	}

		# get character names and numbers
	if (length(grep('\\bCHARLABELS', x, ignore.case=TRUE)) > 0) {
		charlabelsstart <- grep('\\bCHARLABELS', x, ignore.case=TRUE) + 1
		charlabelsend <- grep('\\;$', x[charlabelsstart:length(x)])[1] + charlabelsstart - 2	
		charnames <- na.omit(str_match(x[charlabelsstart:charlabelsend], '\\[(\\d{1,4}|[A-Z]+)[.]*\\]\\s(.+)'))
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
		statelabels <- str_match(x[statelabelsstart:statelabelsend], '[\\t]*\\d+\\s*(.+)$')[,2]
		# sometimes states are on multiple lines:
		if (length(statelabels) != nchar) {
			statelabelsend <- grep('\\;$', x[statelabelsstart:length(x)])[1] + statelabelsstart - 1
			statelabels <- x[statelabelsstart:statelabelsend]
			statelabels <- paste(statelabels, collapse="")
			states_start <- str_locate_all(statelabels, "\\d+[\\t\\s]+[\\']")  # start
			states_end <- str_locate_all(statelabels, "\\t+\\,|\\t+\\;")  # end
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
		
		charstatelabels <- str_match_all(x[charlabelsstart:charlabelsend],
			"(\\t|\\,)\\s(\\d{1,})\\s'(.*?)('\\s\\/(.*?))?(?=\\;(\\s)?$|(\\,\\s\\d{1,}\\s'))")

# charstatelabels[[1]][3, ]
# charstatelabels[[1]][19, ]
# charstatelabels[[1]][1119, ]
# charstatelabels[[1]][2504, ]

		charnums <- as.numeric(charstatelabels[[1]][,3])
		charlabels <- charstatelabels[[1]][,4]
		statelabels <- charstatelabels[[1]][,6]
		res$charlabels <- charlabels
		res$charnums <- charnums
		res$statelabels <- gsub("^\\s|\\s$", "", statelabels)
	}

	class(res) <- c('nex', 'list')
	
	res
}

# system.time(tmp <- read.nex(file = '~/Documents/UT/projects/phenome/data/concatenated.nex'))  # 0.4 seconds

# problem with bertelli_2002.nex

# library(phylobase)

# system.time(tmp2 <- readNexus(file='~/Documents/UT/projects/phenome/data/concatenated.nex', type='data', check.names=FALSE, return.labels=FALSE))  # 2 seconds

# ~4Xs faster than readNexus, and read.nexus.data doesn't work with morphological characters

# should I work on making this integrate with existing programs, like MrBayes, PAUP, etc. 
