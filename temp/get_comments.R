# work with markup comments on chars

devtools::load_all('~/github/nexustools')

# how characters are characterized
# DD = duplicate
# KK = keep
# KKM = keep and modify wording
# XX = cut

setwd("/Users/chadeliason/Documents/UT/projects/phenome")

oldfile <- read.nex("~/Dropbox/phenome dataset/data/2015-09-02/original/final_reordered.nex")
newfile <- read.nex("/Users/chadeliason/Dropbox/phenome dataset/matrix updated 6 8 16.nex")

# read.nex("/Users/chadeliason/Desktop/matrix updated 6 8 16 TEST.nex")

commentfile <- "/Users/chadeliason/Dropbox/phenome dataset/data/2015-09-02/modified/final_reorderedJAC+CME.txt"

# read text
txt <- readLines(commentfile)
txt <- paste0(txt, collapse="\n")
txt <- str_trim(txt)

# find actions, comments associated with characters
tmp <- str_match_all(txt, regex("\n([A-Z\\?]+)(\\d+)", multiline=TRUE, dotall=TRUE))
tmp <- tmp[[1]]
todo <- tmp[,2]
charnum <- as.numeric(tmp[,3])
matches <- tmp[, 1]

# find comments
locs <- str_locate(txt, fixed(matches))
comments <- list()
for (i in seq_along(matches)) {
	start <- locs[i, 1]
	if (i == length(matches)){
		end <- str_locate(txt, "$")[2]
	} else {
		end <- locs[(i+1), 1]	
	}
	newtext <- substr(txt, start, end)	
	comments[[i]] <- str_extract_all(newtext, "(?<=\\{)(.*?)(?=\\})")[[1]]
}
names(comments) <- charnum
comments <- setNames(unlist(comments, use.names=F), rep(charnum, times = sapply(comments, length)))

# merge
res1 <- data.frame(charnum = as.numeric(names(comments)), comment = comments)
res2 <- data.frame(charnum = charnum, todo = todo)
res <- dplyr::left_join(res2, res1, by = "charnum")

# output
write.csv(res, file = "output/regex_extracted.csv")


# extract duplicate characters from comments

# str_match_all(as.character(res1$comment), "(\\bduplic|overlap).*?(\\d+[;,\\s]?.*?(?=[a-z]+))")[27]

# first match duplicate, overlaps
dups <- str_match_all(as.character(res1$comment), "(\\bdupl|overl).*?(\\d+(?:[;,]\\s\\d+)*)")

# duplicate or overlap, etc.

# types <- sapply(dups, "[", i=2)

dups <- lapply(dups, "[", , 3)

dups <- sapply(dups, strsplit, split = "[;,]")

dups <- sapply(dups, unlist)

dups <- sapply(dups, gsub, pattern = "^ ", replacement = "")

names(dups) <- res1$charnum

dups <- na.omit(setNames(unlist(dups, use.names=F), rep(names(dups), times = sapply(dups, length))))

dups <- data.frame(target = as.numeric(names(dups)), duplicate = as.numeric(dups))


write.csv(dups, file = "output/regex_duplicates.csv")



# plotting

# library(igraph)

# g <- graph_from_edgelist(as.matrix(dups), directed=FALSE)

# g <- simplify(g)

# par(mar=c(0,0,0,0))
# plot(g, vertex.color = "lightblue", vertex.label.cex = .8, vertex.size = 0, edge.width = .5, edge.arrow.size=.5)


# do stuff to matrix

drops <- res$charnum[grep("^XX$", res$todo)]

dups



fixed 
table(is.na(oldfile$data))
table(is.na(oldfile[, -drops]$data))
# 268662/(40838+268662)

table(is.na(duplicated(oldfile, map = dups)$data))
# 270173/(41077+270173)

