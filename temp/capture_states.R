# should be able to extract state labels from a chunk of text

x <- allnex.genome$statelabels

# strip out comments
xnew <- gsub('(\\s\\[.*?\\])', '', x)

# capture text between '' or separated by spaces
matches <- str_match_all(xnew, "([^\\s'']+)|'([^']*)'")  #http://stackoverflow.com/questions/366202/regex-for-splitting-a-string-using-space-when-not-surrounded-by-single-or-double

# extract matches
x <- lapply(matches, '[', , 2:3)

# collapse to get list of character states
states <- sapply(1:length(x), function(z) {tryCatch(apply(x[[z]], 1, paste, collapse=''), error=function(e) NULL)})

# check lengths
sapply(states, length)

# find characters with an 'absent' token
states[grep('absent', states)]

# check that number of states found isn't greater than number of states in the data matrix?

# convert absent-multistate to multiple characters

sapply(states, length)[grep('absent', states)]

