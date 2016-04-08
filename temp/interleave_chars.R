# TODO: interleave data based on bone/body region

# i need to write a read.chars function for reading from a PDF

# pattern:

# Character 1: CHARACTER NAME
# 	STATE 1: STATE 1 TEXT
# 	STATE 2: STATE 2 TEXT
#	COMMENT/NOTE TEXT
# Character 2: ...

pdffile <- "/Users/chadeliason/Downloads/monograph\ charlist\ public.pdf"

first = 1
last = 74
syscall <- paste("pdftotext -layout -f ", first, " -l ", last, " '", pdffile, "'", sep="")
system(syscall)
# load and scan newly created text file
txtfile <- gsub('pdf', 'txt', pdffile)
raw <- scan(txtfile, what="", sep="\n")
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

raw2 <- do.call(paste0, list(raw, collapse="\n"))


tmp <- str_match_all(raw2, regex("Character\\s(\\d+)(.*)", multiline=TRUE, dotall=FALSE))

charnums <- as.numeric(tmp[[1]][, 2])

charnames <- tmp[[1]][, 3]

# charnames <- final$charlab


# remove comments
charnames_clean <- str_replace(charnames, "(\\[|\\().*?(\\]|\\))", "")

# remove starting line punctuation, whitespace
charnames_clean <- str_replace(charnames_clean, "^(\\:)?[\\s\\t]*", "")


words <- str_extract_all(charnames_clean, "\\w+")

sort(table(unlist(words)))


# find connections among characters/with similar words in them

pairids <- combn(seq_along(words), m=2)

overlaps <- list()

words_clean <- lapply(seq_along(words), function(x) {
	words[[x]][!words[[x]] %in% c("of", "and", "the", "to", "on", "at", "as")]
})

words_clean <- lapply(words_clean, tolower)

# find and remove one letter words
words_clean <- lapply(seq_along(words_clean), function(x) {
	ss <- !str_detect(words_clean[[x]], "\\b[A-Za-z]{1,1}\\b|\\b\\d{1,4}\\b")
	words_clean[[x]][ss]
})


uniwords <- sort(unique(unlist(words_clean)))

uniwords

# match words using regular expressions:

setwd("/Users/chadeliason/Documents/UT/projects/phenome")

termlist <- read.csv("data/phenome_terms.csv")

head(termlist)

# grep(as.character(tomatch$search.term)[15], uniwords, value=TRUE)

tomatch <- lapply(as.character(termlist$search.term), grep, uniwords, perl=TRUE, value=TRUE)

# tomatch <- list('pelvis' = grep("ili|isch", uniwords, value=TRUE),
# 	'axial' = grep("verte|fem", uniwords, value=TRUE),
# 	'forelimb' = grep("hum|rad", uniwords, value=TRUE))

res <- matrix(nrow=length(tomatch), ncol=length(words_clean))

for (i in seq_along(tomatch)) {
	for (j in seq_along(words_clean)) {
		res[i, j] <- as.numeric(any(tomatch[[i]] %in% words_clean[[j]]))
	}
}

rownames(res) <- paste0("term", seq_len(nrow(termlist)))

colnames(res) <- paste0("char", seq_along(words_clean))

dim(res)

image(res)






# option B: figure out all unique words, and which characters have the words in the name

res <- matrix(NA, nrow=length(uniwords), ncol=length(words_clean))

for (i in seq_along(uniwords)) {
	for (j in seq_along(words_clean)) {
		res[i, j] <- uniwords[i] %in% words_clean[[j]]
	}
}



res <- matrix(as.numeric(res), nrow=nrow(res))

dim(res)

rownames(res) <- paste0("term",seq_len(nrow(termlist)))

paste0("char",seq_along(words_clean)))

rownames(res)

colnames(res)



# i <- sample(1:nrow(res), 150)
# j <- sample(1:ncol(res), 150)
# g <- graph_from_incidence_matrix(res[i, j], directed=FALSE)

library(igraph)

# create graph
g <- graph_from_incidence_matrix(res, directed=FALSE)

# # decompose graph... (??)
# graphs <- decompose.graph(g)

# # figure out how decomposed network graphs are
# graph.sizes <- sapply(seq_along(graphs), function(x) { length(E(graphs[[x]])) })

# plot(graphs[[which.max(graph.sizes)]])

# graphs[[which.max(graph.sizes)]]

# wc <- cluster_fast_greedy(graphs[[18]])

wc <- cluster_fast_greedy(g)

# wc <- cluster_walktrap(g, steps=8)
# wc <- cluster_spinglass(g) # vertex=1
modularity(wc)

barplot(sort(table(membership(wc)),decr=T), xlab="clique name", ylab="number of connections", horiz=F, las=1)


pdf(file="~/Desktop/phenome_communities.pdf")
plot(wc, g, vertex.size=1, vertex.label.cex=.25, edge.width=.25)
dev.off()


mem <- membership(wc)

table(mem)

g

which(mem==25)

plot(induced_subgraph(g, names(which(mem==3))), vertex.size=0)

# max connections
plot(induced_subgraph(g, names(which(mem==which.max(table(mem))))), vertex.size=0, layout=layout_in_circle)
# min connections
plot(induced_subgraph(g, names(which(mem==which.min(table(mem))))), vertex.size=0)


# subset characters within a given cluster

ss <- names(which(mem==which.max(table(mem))))

plot(induced_subgraph(g, ss))

# show list of characters in a give group/subgraph:
charnames_clean[na.omit(as.numeric(str_extract(ss, "(?<=char)\\d+")))]





# option A:

for (i in 1:ncol(pairids)) {
	x1 <- pairids[1, i]
	x2 <- pairids[2, i]
	numboth <- length(unique(c(words_clean[[x1]], words_clean[[x2]])))
	numsame <- length(intersect(words_clean[[x1]], words_clean[[x2]]))
	d <- 1 - (numsame/numboth)
	overlaps[[i]] <- intersect(words_clean[[x1]], words_clean[[x2]])
}


# plot(sort(sapply(overlaps, length), decreasing=T), type='l')

# names(overlaps) <- 1:ncol(pairids)

# numoverlap <- sapply(overlaps, length)

# M <- matrix(NA, ncol=length(words), nrow=length(words))
# M[lower.tri(M, diag=F)] <- numoverlap
# M[upper.tri(M, diag=F)] <- rev(numoverlap)
# diag(M) <- sapply(words, length)
# M[1:5, 1:5]


# I want distance to be number of unique words out of total number of unique words

length(numoverlap)

dim(pairids)

sort(unique(unlist(words_clean)))

# I want similarity to be the number of shared words divided by total number of words

sort(unique(c(words_clean[[15]], words_clean[[20]])))


plot(hclust(as.dist(M[1:10, 1:10])))



ov <- overlaps[sapply(overlaps, length) > 0]


set1 <- pairids[1, as.numeric(names(ov))]
set2 <- pairids[2, as.numeric(names(ov))]

library(igraph)

verts <- sort(unique(c(set1, set2)))

edgs <- data.frame(from=set1, to=set2)

edgs


words_clean[[472]]
words_clean[[475]]

sort(unique(unlist(words_clean)))

g <- graph_from_data_frame(edgs, vertices=verts, directed=FALSE)

grps <- max_cliques(g)

com <- components(g)

# plot(induced_subgraph(g, v=com$membership==which.max(com$csize)))
plot(induced_subgraph(g, v=c(1, 3, 5, 8, 10, 15)))

words_clean[1]
words_clean[15]




# visualize connections with a network graph...



# group by connectivity, also by a set of predefine morphometric categories/terms (ontology)



# another way to visualize things

library(wordcloud)

wordcloud(charnames_clean)
