# read characters from a PDF or TXT file

pdffile <- "/Users/chadeliason/Dropbox/PDFs/LIVEZEY 2006_PHYLOGENY OF NEORNITHES_matrix.pdf"
first <- 28
last <- 441
ntax <- 188
nchar <- 2954

syscall <- paste("pdftotext -f ", first, " -l ", last, " '", pdffile, "'", sep="")

# NB removing layout was a better option for Livezey and Zusi (2006)

system(syscall)

# load and scan newly created text file

txtfile <- gsub('pdf', 'txt', pdffile)

raw <- scan(txtfile, what="", sep="\n")

# remove/clean up exported text files

system(paste("rm '", txtfile, "'", sep=""))

raw2 <- do.call(paste0, list(raw, collapse="\n"))

raw2

library(stringi)
library(stringr)

?str_match_all

tmp <- str_match_all(raw2, regex("\n(\\d{4,4})\\.(.*?)(?=\n\\d{4,4}\\.)", dotall=TRUE, multiline=TRUE))

setdiff(as.numeric(tmp[[1]][, 2]), 1:2954)

charnums <- as.numeric(tmp[[1]][, 2])

charlabs <- tmp[[1]][, 3]

charlabs <- charlabs[match(1:2954, charnums)]

charnums <- charnums[match(1:2954, charnums)]

diff(charnums)
