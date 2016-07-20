library(polycor)
library(psych)
# library(multicore)
library(parallel)

devtools::load_all('~/github/nexustools')

setwd("/Users/chadeliason/Documents/UT/projects/phenome")

test <- read.nex("output/final_reordered.nex")

test$charpart

test <- subset(test, charpartition=="hindlimb")
mat <- test$data
pairids <- combn(1:ncol(mat), m=2)
dim(pairids)

# will take ~2.6 hours (or 30 minutes using all cores)
(5/8000) * ncol(pairids)

# speed up using C?
# took 106 seconds
system.time(res <- mclapply(1:ncol(pairids), function(x) {
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
}, mc.cores = 8))

res[grepl("Error", res)] <- NA

res <- as.numeric(res)

plot(sort(res), type='s')
abline(h=1, lty=2)
abline(h=-1, lty=2)

saveRDS(res, file="output/polychor_hindlimb_results.RDS")

# Test to see which characters these are
which(res < -0.9)
pairids[, 104619]

image(apply(cbind(mat[, 451], mat[, 453]), 2, as.numeric))
plot(test[, c(450, 452)], legend.pos="right")

# now what to do with polymorphisms?

# i need to see how many of these are real duplicates identified from our work



