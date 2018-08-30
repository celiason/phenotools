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

truedup <- read.csv("output/regex_duplicates.csv")
id1 <- paste0(pmin(truedup[,"target"], truedup[,"duplicate"]), "--", pmax(truedup[,"target"], truedup[,"duplicate"]))
id1 <- unique(id1)

# subset chars we have assessed thus far: 1856 - 2027
test2 <- test[, test$charnum %in% 1856:2027]

# keep very highly correlated characters
picks <- which(res > 0.99 | res < -.99)

# create unique sets of potential duplicates
sets <- pairids[, picks]
dups <- cbind(test$charnum[sets[1, ]], test$charnum[sets[2, ]])
dups <- dups[dups[,1] %in% 1856:2027 & dups[,2] %in% 1856:2027, ]
id2 <- paste0(pmin(dups[,1], dups[,2]), "--", pmax(dups[,1], dups[,2]))
id2 <- unique(id2)

# compute precision and recall
fneg <- sum(!id1 %in% id2)  # false negs
tpos <- sum(id2 %in% id1)  # true pos
fpos <- sum(!id2 %in% id1) # false pos
tpos / (tpos + fpos)  # precision
tpos / (tpos + fneg)  # recall

