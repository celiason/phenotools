# merge taxa

dat2 <- matrix(dat, nrow=nrow(dat), dimnames=list(rownames(dat), NULL))

dat2 <- ifelse(!is.na(dat2), 1, NA)

head(dat2)

# find overlapping scorings
any(dat2[5,]==dat2[6,], na.rm=TRUE)


# find taxa with genus only

id.onlygenus <- grep('spp', rownames(dat2))

genus <- str_extract(rownames(dat2)[id.onlygenus], '[A-Z][a-z]+')

genus_spp <- str_extract(rownames(dat2)[id.onlygenus], '[A-Z][a-z]+\\ssp[p]*')

# see if scorings overlap


id.onlygenus

table(dataset[which(dat2[id.onlygenus[7], ]==1)])

dataset[dat2[id.onlygenus[1], ]]

id.onlygenus

# input = taxa without species info
# output = overlap (same characters scored) with congenerics in data matrix

tmp <- list()

for (i in seq_along(genus_spp)) {
	matches <- grep(genus[i], rownames(dat2))
	matches <- setdiff(matches, id.onlygenus)
	id <- id.onlygenus[i]
	if (length(matches) > 1) {
		res <- sapply(1:length(matches), function(x) {any(dat2[id, ] == dat2[matches[x], ], na.rm=TRUE)})
		names(res) <- rownames(dat2)[matches]
		tmp[[i]] <- res
	} else {tmp[[i]] <- 'only one case'}
}

names(tmp) <- genus_spp

tmp

# only Eudromia_spp and Nothura_spp overlap with other datasets


grep(genus[10], rownames(dat2), value=FALSE)

any(dat2[5,]==dat2[6,], na.rm=TRUE)

