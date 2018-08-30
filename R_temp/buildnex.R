#' Read a nexus data file
#' 
#' A function that reads data stored in nexus files
#' NOTE: it is necessary to first install poppler with, e.g., `brew install poppler`
#' 
#' @param file (required) path to nexus file
#' @param missing character representing missing data
#' @param gap character representing inapplicable/incomporable data
#' 
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' 
#' @examples \dontrun{
#' tmp <- buildnex('/Users/chadeliason/github/nexustools/example/Bertelli_Chiappe_2005.pdf', ntax = 34,
#'	nchar = 63, first = 22, last = 22)
#' plot(tmp, legend.pos = 'top')
#' }
#' 
#' @author Chad Eliason \email{chad_eliason@@utexas.edu}
#'
buildnex <- function(pdffile, ntax, nchar, first, last, missing = '?', gap = '-') {

# convert pdf to text file
syscall <- paste("pdftotext -layout -f ", first, " -l ", last, " '", pdffile, "'", sep="")

# might need to give option for this in case it doesn't read well
# turning off `-layout` can help
# syscall <- paste("pdftotext -f ", first, " -l ", last, " '", pdffile, "'", sep="")

# removing layout was a better option for Livezey and Zusi (2006)

system(syscall)

# load and scan newly created text file

txtfile <- gsub('pdf', 'txt', pdffile)

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



# [ ] need to be able to look in a few following lines to see if data are there

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
nchars==nchar

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

res <- list(data = datamatrix, taxlabels = taxlabels, missing = missing, gap = gap, symbols = symbols)

class(res) <- c('nex', 'list')

res

}


# DOESN'T WORK

# buildnex(pdffile = '/Users/chadeliason/github/nexustools/example/Bledsoe_1988.pdf', ntax=9, nchar=83, first = 17, last = 17)
# buildnex(pdffile = '/Users/chadeliason/github/nexustools/example/Bertelli_2014.pdf', ntax = 56, nchar = 157, first = 24, last = 25)
# buildnex(pdffile <- "/Users/chadeliason/Dropbox/Articles/LIVEZEY 2006_PHYLOGENY OF NEORNITHES_matrix.pdf", first=28, last=441, ntax=188, nchar=2954)

# WORKS

# buildnex(pdffile = '/Users/chadeliason/github/nexustools/example/Bertelli_Chiappe_2005.pdf', ntax = 34, nchar = 63, first = 22, last = 22)


