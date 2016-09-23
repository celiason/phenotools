#' Get trait correlations for all pairwise combinations of discrete traits in nexus
#' object 'x'
#'
traitcor <- function(x, parallel = FALSE, cores = 2) {

	library(polycor)
	library(psych)
	# library(multicore)
	library(parallel)
	devtools::load_all('~/github/nexustools')

	mat <- x$data

	pairids <- combn(1:ncol(mat), m=2)

	# will take ~2.6 hours (or 30 minutes using all cores)
	# (5/8000) * ncol(pairids)

	if (parallel) {

		res <- mclapply(1:ncol(pairids), function(x) {
			id1 <- pairids[1, x]
			id2 <- pairids[2, x]
			char1 <- mat[, id1]
			char2 <- mat[, id2]
			# remove NAs
			M <- cbind(char1, char2)
			M <- na.omit(M)
			M <- apply(M, 2, as.numeric)
			# get correlation
			try(polychoric(M)$rho[1,2])
		}, mc.cores = cores)
		} else {

		res <- mclapply(1:ncol(pairids), function(x) {
			id1 <- pairids[1, x]
			id2 <- pairids[2, x]
			char1 <- mat[, id1]
			char2 <- mat[, id2]
			# remove NAs
			M <- cbind(char1, char2)
			M <- na.omit(M)
			M <- apply(M, 2, as.numeric)
			# get correlation
			try(polychoric(M)$rho[1,2])
		}, mc.cores = 8)

	}

	res[grepl("Error", res)] <- NA

	res <- as.numeric(res)

	# plot(sort(res), type='s')
	# abline(h=1, lty=2)
	# abline(h=-1, lty=2)

	list(cor = res, trait1 = pairids[1, ], trait2 = pairids[2, ])

}
