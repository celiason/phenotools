#' Read a nexus data file
#'
#' a function that reads data stored in nexus files
#' NOTE: it is necessary to first install poppler with, e.g., `brew install poppler`
#'
#' @param file (required) path to nexus file
#' @param missing character representing missing data
#' @param gap character representing inapplicable/incomporable data
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' @examples \dontrun{
#'
#' tmp1 <- buildnex(pdffile = 'example/Johnston_2011.pdf', ntax = 13, nchar = 57, first = 21, last = 24)
#' tmp2 <- buildnex(pdffile = 'example/Bertelli_Chiappe_2005.pdf', ntax = 13, nchar = 63, first = 20, last = 22)
#' plot(tmp1, legend.pos = 'top')
#' plot(tmp2, legend.pos = 'top')
#'
#' @author Chad Eliason \email{chad_eliason@@utexas.edu}
#'
buildnex <- function(pdffile, ntax, nchar, first, last, missing = '?', gap = '-') {

# pdffile <- 'example/Bertelli_Chiappe_2005.pdf'
# pdffile <- 'example/Johnston_2011.pdf'
# pdffile <- 'example/McKitrick 1991.pdf'
# pdffile <- '/User/chadeliason/Desktop/Bledsoe.1988.pdf' #, nchar=83, ntax=9, first=17, last=17) 
# pdffile <- '/Users/chadeliason/Desktop/Bledsoe.1988.pdf'
# nchar <- 83
# ntax <- 9
# first <- 17
# last <- 17

# first = 71
# last = 76
# pdffile <- 'example/Bertelli_2014.pdf' #, ntax = 56, nchar = 157, first = 24, last = 25)
# first <- 24
# last <- 25

# pdffile = 'example/Mayr_2003.pdf' #, ntax = 46, nchar = 148, first = 20, last = 24)
# ntax = 46
# nchar = 148
# first = 20
# last = 24

# pdffile <- '/Users/chadeliason/Downloads/Mayr_2011.pdf'
# first = 11
# last = 13
# nchar = 153
# ntax = 47


	syscall <- paste("pdftotext -layout -f ", first, " -l ", last, " '", pdffile, "'", sep="")

	# convert pdf to text file
	system(syscall)

	# load newly created text file

	txtfile <- gsub('pdf', 'txt', pdffile)

	x <- scan(txtfile, what="", sep="\n")

	# clean up text file
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

	matstart <- grep('table|data [set]|character|matrix', x, ignore.case=TRUE)[1]

	# matend <- 

	mat <- x[(matstart+1):length(x)]

	mat <- gsub('^\\s+', '', mat)

	mat <- paste0(mat, collapse="\n")

	# str_match_all(mat, regex('\\n([A-Za-z][A-Za-z\\.\\-0-9]+[\\s]*[A-Za-z\\.\\-0-9]+)(.+)'))


	# data might be in continuous rows (tax1 tax2 tax3)
	# or discontinuous (e.g., tax1 tax2 tax3; tax1 tax2 tax3; tax1 tax2 tax3; ...)

	# this has problems if taxa labels aren't in the form 'H. sapiens'

	# should read contiguous lines

	# and how do we make sure these are valid taxon names?

	# some datasets have taxa at end of matrix row

	datamatches <- str_match_all(mat, regex('\\n([A-Za-z]+[\\._\\s]*[a-z]*[\\.]*)\\s+([A-Z0-9\\?\\-]{1}(?![a-z]).+)'))


	# Option 2: inputting taxon names maybe as the function argument (e.g., Bledsoe 1988 PDF document)
	# datamatches <- str_match_all(mat, regex('\\n([A-Za-z]+[\\._\\s]*[a-z]*[\\.]*)\\s+([A-Z0-9\\?\\-]{1}(?![a-z]).+)'))
	# taxlabels <- c('Ancestor', 'Dinomithidae', 'Apteryx', 'Casuarius', 'Dromaius', 'Rhea', 'Struthio', 'Aepyornis', 'Dromo(m|rn)ithidae')
	# sapply(taxlabels, grep, x)
	# mat[c(8, 96, 242)]
	# some of the data are on separate rows, but column position should be sequential...
	# str_locate_all(mat, '([a-z])')


	# Option 3: getting data from Mayr (2011): data on multiple lines after taxon listed
	# this specifies the taxon labels should be like 'Passeriformes' without spaces, capitalized
	# and that the data matrix cannot have capital letters and must be at least nchar long

	# tomatch <- paste0('([A-Z][a-z\\_]*)([0-9\\?a-z\\-\\n]{', nchar, ',})', collapse='')
	# datamatches <- str_match_all(mat, regex(tomatch))
	# dat <- gsub('\n', '', datamatches[[1]][,3])

	# get 1-nchar columns

	# tmp <- str_match_all(dat, '(\\d{1}|[\\[\\(]\\d{1,3}[\\]\\)]|\\-|\\?|\\w|\\–\\s)')
	# prob <- sapply(tmp, nrow) != nchar
	# if (any(prob)) {
	# 	warning('Found extra characters for taxon: ', taxlabels[which(prob)], ', trimming extras--check PDF')
	# }
	# dat <- t(sapply(tmp, '[', i=1:nchar, j=2))



	taxlabels <- datamatches[[1]][,2]
	taxlabels <- gsub('\\s+$', '', taxlabels, perl=TRUE)
	# taxlabels <- gsub(' ', '_', taxlabels)
	taxlabels <- unique(taxlabels)

	# matches <- paste0('(', paste0(taxlabels, collapse='|'), ')')

	x2 <- x[matstart:length(x)]

	id <- lapply(taxlabels, grep, x2)

	names(id) <- taxlabels

	# this replaces taxon names with same number of spaces to keep column locations the same (in the case of missing data - white space)
	
	strlens <- str_length(names(id))

	mat <- lapply(seq_along(id), function(i) { gsub(names(id)[i], paste0(rep(' ', strlens[i]), collapse=''), x2[id[[i]]]) } )

	names(mat) <- taxlabels

	# TODO: WORK ON

	mat <- lapply(mat, paste0, collapse='')

	# if (is.matrix(mat)) {
		# 
	# } else {
	mat <- do.call(rbind, mat)
	# }


# sapply(mat, paste0, collapse='')

	# mat <- do.call(rbind, mat)

	# names(mat) <- taxlabels

	# extract data/character scores

	# datamatches <- str_match_all(mat, '(\\d{1}|[\\s\\[\\(]\\d{1,3}[\\]\\)\\s]|\\-|\\?|\\w|\\–\\s)')

	# datamatches <- str_match_all(mat, '(\\d{1}|\\[\\d{1,3}\\]|\\-|\\?|\\w|\\–\\s)')

	# check if any columns lengths differ (i.e. if different number of scored characters for a taxon)

	# find columns positions


# str_length(mat)

# tmp <- str_match_all(mat, '(\\d{1,3}|[\\[\\(]\\d{1,3}[\\]\\)]|\\-|\\?|\\w|\\–\\s)')

# sapply(tmp, nrow)

	colpos <- str_locate_all(mat, '(\\d{1,3}|[\\[\\(]\\d{1,3}[\\]\\)]|\\-|\\?|\\w|\\–\\s)')

	# colpos <- str_locate_all(mat, '(\\d{1}|[\\s\\[\\(]\\d{1,3}[\\]\\)\\s]|\\-|\\?|\\w|\\–\\s)')

	# colpos <- str_locate_all(mat, '(\\d{1}|[\\[\\(]\\d{1,3}[\\]\\)]|\\-|\\?|[A-Z]|\\–)')

	# nchar <- 157

	# only 3 possible missing characters (e.g., spaces)
	
	charlens <- sapply(colpos, nrow)

	keep <- abs(charlens - nchar) < 3

	colpos <- colpos[keep]
	mat <- mat[keep]
	taxlabels <- taxlabels[keep]

	charlens <- sapply(colpos, nrow)

	id <- which(abs(charlens - nchar) > 0)

	# problem: if page title gets between data matrix on separate pages
	# maybe look for a contiguous block of character scorings (e.g., 00??---)?

	# problem: if some characters have space as a scoring
	# does this work correctly?
	if (any(id)) {
		# switchpoint <- which(abs(apply(cbind(colpos[[id]][, 2], colpos[[id+1]][, 2]), 1, diff)) > 1)[1]
		# df <- expand.grid(colpos[[8]][,2], colpos[[9]][,2])
		# df[,3] <- abs(df[,2]-df[,1])
		# head(df)
		# tmp <- str_locate_all(mat2, '\\|(.*?)(?=\\||$)')
		# idmaxcols <- which.max(sapply(colpos, nrow))
		# compare to row above, set equal to?
		# start <- setdiff(colpos[[id+1]][,'start'], colpos[[id]][,'start'])
		# id
		# dim(colpos[[id]])
		# dim(colpos[[id+1]])
		# mat[[id]][switchpoint, 2]
		# end <- setdiff(colpos[[id+1]][,'end'], colpos[[id]][,'end'])
		# start
		# end
		# colpos[[id+1]]
		colpos[[id]] <- colpos[[id+1]]
	}

	# extract scores

	# df <- expand.grid(x=1:19, y=seq_along(taxlabels))


	# dim(df)

	# head(df)

	# rep(taxlabels, ncol(mat))

	# , tax = rep(taxlabels, ncol(mat))

	# head(df)

	data <- matrix(NA, nrow=length(taxlabels), ncol=nrow(colpos[[1]]))

	for (i in seq_along(taxlabels)) {
		for (j in 1:nrow(colpos[[1]])) {
			data[i, j] <- substr(mat[i], colpos[[i]][j,1], colpos[[i]][j,2])
		}
	}

	# sapply(1:ncol(mat), function(x) { grep('(\\d{1}|\\[\\d{1,3}\\]|\\-|\\?|\\w|\\–)', mat[, x])  } )


	# problem character is coded as a space (probably rare?)


	# mat

	# 	data <- t(do.call(cbind, data))

	if (nrow(data) != ntax) {
		warning('Number of rows in dataset not equal to specified number of taxa')
	}

	if (ncol(data) != nchar) {
		warning('Number of columns in dataset not equal to specified number of characters')
	}

	symbols <- paste(unique(unlist(strsplit(gsub('[\\[\\]\\(\\)\\?\\-\\–]', '', sort(unique(as.vector(data))), perl=TRUE), ""))),collapse="")

	# change blanks to gaps
	data <- gsub('^ $', '-', data)

	# strip whitespace
	data <- gsub(' ', '', data)

	data <- ifelse(data == '–', '-', data)
	
	data <- ifelse(data == missing, NA, data)

	taxlabels <- gsub(' ', '_', taxlabels)

	res <- list(data = data, taxlabels = taxlabels, missing = missing, gap = gap, symbols = symbols)

	if (exists('charnums')) {
		res$charnums <- charnums
	}
	if (exists('charlabels')) {
		res$charlabels <- charlabels
	}

	class(res) <- c('nex', 'list')

	res

}
