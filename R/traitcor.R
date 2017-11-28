#' Get trait correlations for all pairwise combinations of discrete traits in nexus file
#' x = 'nexus' object 
#'
traitcor <- function(x, parallel = FALSE, cores = 1, method = c('hamming', 'polycor')) {
	library(parallel)
	library(devtools)
	library(polycor)
	library(psych)
	library(pbapply)
	library(stringdist)
	load_all('~/github/nexustools')	
	method <- match.arg(method)
	mat <- x$data
	# convert '?' or '-' scorings to NA:
	mat <- gsub("\\?|\\-", NA, mat)
	# get all pairs of characters
	pairids <- combn(1:ncol(mat), m=2)
	# calculate string distances by one of 2 methods
	if (method=="polycor") {
		# library(multicore)
		# will take ~2.6 hours (or 30 minutes using all cores)
		# (5/8000) * ncol(pairids)
		if (parallel) {
			dmat <- mclapply(1:ncol(pairids), function(i) {
				id1 <- pairids[1, i]
				id2 <- pairids[2, i]
				char1 <- mat[, id1]
				char2 <- mat[, id2]
				M <- cbind(char1, char2)
				# remove NAs
				M <- na.omit(M)
				M <- apply(M, 2, as.numeric)
				# get correlations
				try(polychoric(M)$rho[1,2])
			}, mc.cores = cores)
		} else {
			dmat <- mclapply(1:ncol(pairids), function(i) {
				id1 <- pairids[1, i]
				id2 <- pairids[2, i]
				char1 <- mat[, id1]
				char2 <- mat[, id2]
				# remove NAs
				M <- cbind(char1, char2)
				M <- na.omit(M)
				M <- apply(M, 2, as.numeric)
				# get correlation
				try(polychoric(M)$rho[1,2])
			}, mc.cores = cores)
		}
	}
	if (method == "hamming") {
		# dmat <- mclapply(1:ncol(pairids), function(i) {
		dmat <- pblapply(1:ncol(pairids), function(i) {
			id1 <- pairids[1, i]
			id2 <- pairids[2, i]
			char1 <- mat[, id1]
			char2 <- mat[, id2]
			M <- na.omit(cbind(char1,char2))
			wt <- nrow(M)
			M <- apply(M, 2, paste0, collapse="")
			# M <- c('0101', '0111')
			sdist <- as.numeric(stringdistmatrix(M, method="hamming"))
			# standardize by number of characters
			# 0011 vs 1100 will give std. distance of 1
			# 0011 vs 0011 will give std. distance of 0
			# 0101 vs 0111 will give std. distance of 0.25
			sdist <- sdist/str_length(M[1])
			sdist <- abs(sdist-0.5)
			sdist <- sdist/0.5
			# NB: might want to weight by number of co-scored characters @done
			# i.e. if only one species scored overlap, very weak evidence for similar traits
			c(sdist=sdist, weight=wt)
			}) #, mc.cores = cores)
	}
	dmat[grepl("Error", dmat)] <- NA
	cors <- sapply(dmat, "[[", "sdist")
	weights <- sapply(dmat, "[[", "weight")
	res <- data.frame(cor = cors, weight = weights, t1 = pairids[1,], t2 = pairids[2, ])
	res[res$cor!=Inf,]
	# list(cor = cors, weight = weights, pair = pairids[1:2, ], charnames = x$charlabels[pairids[1:2, ]])
}
