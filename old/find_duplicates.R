# check read files

library(phylobase)
library(stringdist)


setwd("/Users/chadeliason/Documents/UT/projects/phenome")

# files <- list.files("data", pattern="\\.nex", full.names=T)

# test <- lapply(files, readNexus, type='data', simplify=T, return.labels=F)

dat <- readNexus('data/concatenated.nex', type='data', return.labels=FALSE)

dat <- data.matrix(dat)



# look for pairwise correlations

# only want comparisons BETWEEN datasets, not within

dataset <- str_extract(colnames(dat), perl("\\w+\\_\\d{4}"))


# strip comment from character names
oldnames <- colnames(dat)
newnames <- str_extract(colnames(dat), perl(".+(?=\\.\\.\\w+\\_\\d{4}\\.nex\\.trait\\.)"))

pairids <- combn(1:ncol(dat), m=2)

# different datasets
id <- which(dataset[pairids[1,]] != dataset[pairids[2,]])

newpairids <- pairids[, id]


# text distances

stringdists <- numeric(length = ncol(newpairids))

# takes 200 seconds for 2.2M comparisons
# system.time(
# 	for (i in 1:ncol(newpairids)) {
# 	str1 <- newnames[newpairids[1,i]]
# 	str2 <- newnames[newpairids[2,i]]
# 	stringdists[i] <- stringdist(str1, str2, method="jw")
# }
# )

# saveRDS(stringdists, file="output/stringdists.RDS")

stringdists <- readRDS(file='output/stringdists.RDS')

hist(stringdists, breaks=50)

names(stringdists) <- 1:length(stringdists)

stringdists.sorted <- sort(stringdists)

sortednames <- as.numeric(names(stringdists.sorted))

## not sure what cutoff to use for text distances..

sset <- sortednames[stringdists.sorted < 0.25]  # try 0.3

sset.dist <- round(stringdists.sorted[stringdists.sorted < 0.25], 4)

# TODO: make it so the set difference between species names is also returned

out <- sapply(seq_along(sset), function(x) {
x <- 1
	id1 <- newpairids[1, sset[x]]
	id2 <- newpairids[2, sset[x]]
	spp1 <- names(na.omit(dat[,id1]))
	spp2 <- names(na.omit(dat[,id2]))
	not1 <- setdiff(spp1, spp2)
	not2 <- setdiff(spp2, spp1)
	both <- intersect(spp1, spp2)
	paste('\n---------\ntrait pair: ', oldnames[id1], oldnames[id2],
	'\nstring distance = ', sset.dist[x],
	'\nmissing from species 1:', not1,
	'\nmissing from species 2:', not2,
	'scored in both species:', both, sep="\n")
})

# cat(x)

x <- sapply(seq_along(sset), function(x) {
	id1 <- newpairids[1, sset[x]]
	id2 <- newpairids[2, sset[x]]
	spp1 <- names(na.omit(dat[,id1]))
	spp2 <- names(na.omit(dat[,id2]))
	not1 <- setdiff(spp2, spp1)
	not2 <- setdiff(spp1, spp2)
	both <- intersect(spp1, spp2)
	paste('--------\n pair', x,
		' (string distance = ', sset.dist[x], ') \n',
		oldnames[newpairids[1, sset[x]]], '\n',
		oldnames[newpairids[2, sset[x]]], '\n',
		'\nmissing in trait1:\n ', paste(not1, collapse='\n'),
		'\n\nmissing in trait2:\n ', paste(not2, collapse='\n'),
		'\n\nspecies overlapping:\n ', paste(both, collapse='\n'), '\n',
		sep="")
})

cat(x, fill=TRUE)

cat(x, file="~/Desktop/allmatches.txt")



## check if trait scorings are co-distributed (i.e. taxa scored same for each trait)

overlap <- sapply(seq_along(sset), function(x) {
	# compute number of taxa with scorings for the same trait
	nrow(na.omit(dat[, newpairids[, sset[x]]]))
	# id1 <- newpairids[1, sset[x]]
	# id2 <- newpairids[2, sset[2]]
	# length(which(complete.cases(cbind(dat[,id1], dat[,id1]))==TRUE))
	})

overlap

table(overlap)

id <- sset[which(overlap > 2)]

# calculate pairwise correlations
res <- sapply(seq_along(id), function(x) {
	# id1 <- newpairids[1, id[x]]
	# id2 <- newpairids[2, id[x]]
	# try(cor(dat[,id1], dat[,id2], use='complete.obs'))
	tmp <- dat[, newpairids[, id[x]]]
	vars <- diag(var(na.omit(tmp)))
	if (all(vars==0)) {
		'1'
	}
	if (any(vars==0)) {
		'cannot calculate'
	} else {
	cor(na.omit(tmp))[2,1]
}
})

id <- id[res=="1" | res=="-1"]

res <- res[res=="1" | res=="-1"]

res

id

match(id, names(stringdists.sorted))


x <- sapply(seq_along(id), function(x) {
	paste('--------\n pair', x, ' \n', oldnames[newpairids[1, id[x]]], '\n',
	oldnames[newpairids[2, id[x]]], '\n', sep="")
})

cat(x, file="~/Desktop/matches+codistr.txt")



