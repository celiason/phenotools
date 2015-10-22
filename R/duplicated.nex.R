#' Find duplicate traits in nexus files
#'
#' a function that uses fuzzy text matching to output a list of potentially redundant characters
#' in a `nex` object, waits for user input to confirm, and then handles merging of characters and
#' outputs a new `nex` object
#'
#' @param x (required) a nexus data object
#' @param cutoff output potential matches under this value for the string distance between characters
#' @param map list of numeric vectors with original character ids mapped to their new ones
#' @param force whether to force character dropping and avoid data checks (not advised)
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' @examples \dontrun{
#'
#' x <- read.nex(file='example/toy1.nex')
#' map <- list(c(2, 4), c(156, 157))
#' duplicated.nex(x, method = 'user', map = map)
#' duplicated.nex(x, method = 'automated', map = map)
#' duplicated(x, cutoff = 0.01, method = 'auto')
#' duplicated(x, method = 'user', map=list(c(2,3)))
#'
#' @author Chad Eliason \email{chad_eliason@@utexas.edu}
#'
duplicated.nex <- function(x, cutoff = 0.25, map = NULL, force = FALSE, n = 50, train = TRUE) {

cleantext <- function(x) {
    x <- str_replace_all(x, "[\\s]*[\\[\\(].+[\\]\\)][\\s]*", "")  # remove comments in square brackets
    x <- sapply(x, tolower)  # make lowercase
    x <- gsub("\\s{2,}", " ", x)  # remove extra spaces
    x <- gsub("^\\s", "", x)  # take out beginning spaces
    x <- gsub("\\s$", "", x)
    x <- removePunctuation(x)
    x <- gsub("\\bthe\\b", "", x)
    x <- gsub("\\bin\\b", "", x)
    x <- gsub("\\bat\\b", "", x)
    x <- gsub("\\ba\\b", "", x)
    # x <- removeWords(x, stopwords("english"))  # take out common words (if, to, the...)
    names(x) <- NULL
    x
}

# workflow:
# string distance calculations 
# training/prediction
# merging data

library(stringdist)
library(tm)
library(MASS)

# provide map of duplicated characters
if (!is.null(map)) {
  pairids <- matrix(unlist(map), ncol=length(map))
}

# automated discovery of duplicate characters
if (is.null(map)) {
  oldcharnames <- x$charlabels
  oldstatenames <- x$statelabels
  newstatenames <- cleantext(oldstatenames)
  newcharnames <- cleantext(oldcharnames)

# [ ] look for duplicates within character names (e.g., 'From Bourdon et al. (2009a: char. 66)') 
# [x] paste together character labels and state labels, look for distances between them

  # generate all possible pairs of character combinations
  pairids <- combn(seq_along(newnames), m=2)

  # only want comparisons BETWEEN datasets/character types, not within:
  file <- x$file
  id <- file[pairids[1, ]] != file[pairids[2, ]]
  pairids <- pairids[, id]

  # only look within character partitions:
  if (!is.null(x$charpartition) & length(unique(x$charpartition)) > 1) {
    charpart <- x$charpartition
    id <- charpart[pairids[1, ]] != charpart[pairids[2, ]]
    pairids <- pairids[, id]
  }

  # calculate text distances (takes ~200 seconds for 2.2M comparisons):
  stringdists <- numeric(length = ncol(pairids))
  for (i in 1:ncol(pairids)) {
    str1a <- newcharnames[pairids[1,i]]
    str1b <- newstatenames[pairids[1,i]]
    str2a <- newcharnames[pairids[2,i]]
    str2b <- newstatenames[pairids[2,i]]
    stringdist1 <- stringdist(str1a, str2a, method="jw")
    stringdist2 <- stringdist(str1b, str2b, method="jw")
    stringdists[i] <- stringdist1 + stringdist2
  }
  names(stringdists) <- 1:length(stringdists)
  sset.dist <- sort(stringdists)
  sset <- as.numeric(names(stringdists.sorted))

  if (train) {
    # sample evenly over range of string distances
    sscut <- stringdists[stringdists < cutoff]
    # sscut <- stringdists[stringdists < 30]
    sscut <- sort(sscut)
    ss <- round(seq(1, length(sscut), length = n))
    sstrain <- sscut[ss]
    # sstrain <- sample(sscut, size = n, replace=F)
    sset.dist <- sstrain
    sset <- as.numeric(names(sset.dist))
  }

  # not sure what cutoff to use for text distances..
  # sset <- sortednames[stringdists.sorted < cutoff]  # try 0.3
  # sset.dist <- round(stringdists.sorted[stringdists.sorted < cutoff], 4)  # is this needed?
  if (length(sset)==0) {
    stop('No putative matches found.')
  }
  # return number of matches and wait for user input
  cat('Found ', length(sset), ' putative matches. Continue with user matching (y/n)?\n')
  continue <- scan(n=1, what='character')
  if (continue=='n') {
    stop('Function terminated by user')
  }
  
  answer <- character(length = length(sset))

  # TODO: [x] maybe use a different string distance metric (e.g., 'vomer' AND 'vomer extending laterally' give dist = 0)
  # loop over pairs, wait for user input confirming redundant traits
  # sset = integers corresponding to row of pairids as possible matches
  # sset.dist = string distances
  for (i in seq_along(sset.dist)) {
    id1 <- pairids[1, sset[i]]
    id2 <- pairids[2, sset[i]]
    # print pairs of characters:
    cat('\n-------\nTrait pair', i, ' (string distance = ', sset.dist[i], ') \n\nCHARLABELS:\n\n',
      oldcharnames[id1], ' (', file[id1], ', character ', x$charnums[id1],')\n\n',
      oldcharnames[id2], ' (', file[id2], ', character ', x$charnums[id2],')\n\nSTATELABELS:\n\n',
      x$statelabels[id1], '\n\n',
      x$statelabels[id2], '\n----------\n\n', sep="")
    cat('Are these the same traits (y/n/N/c)\n')
    # wait for user input
    answer[i] <- scan(n = 1, what = 'character')
    if (answer[i] == 'N') {
      message('User specified that remaining traits are not duplicates')
      answer[i:length(answer)] <- 'n'
      break
    }
    if (answer[i] == 'c') {
      stop('Function terminated by user')
    }
  }
  if (length(answer) == 0) {
    stop('Nothing selected to merge')
  }
  # get duplicated pairs
  ssetnew <- sset[answer=='y']
  pairids <- as.matrix(pairids[, ssetnew])
}


# penalized logistic regression (can estimate models with perfect separation)
# library(arm)
# df <- data.frame(match=ifelse(answer=="y", 1, 0)[1:n], dist=sset.dist[1:n])
# fit <- bayesglm(match ~ dist, data = df, family="binomial")
# visreg(fit)

# linear discrimant analysis
lda1 <- lda(match~dist, data=df)
pred <- predict(lda1, data.frame(dist=stringdists))


# set cutoff for probability of match

# greater than 90% probability of being a match
matchid <- as.numeric(names(pred$posterior[,2] [ pred$posterior[, 2] > 0.5 ] ))

dups <- pairids[, matchid]

# do the actual dropping and merging of characters
res <- x
drops <- numeric(length = ncol(dups))

for (i in 1:ncol(dups)) {

  id <- dups[, i]

  # get character scorings
  scores <- x$data[, id]
  scores1 <- scores[,1]
  scores2 <- scores[,2]
  
  # create reverse scores (in the case where someone scores 0 & 1 and someone else 1 & 0 for the same structure)
  scores2rev <- factor(scores2)
  levels(scores2rev) <- rev(levels(scores2rev))
  scores2rev <- as.character(scores2rev)
  scorings <- colSums(!is.na(scores))
  cc <- complete.cases(scores)
  n.overlap <- nrow(na.omit(scores))  # number of taxa scored for both traits
  
  # Case 1 - Same number of taxa scored
  if (diff(scorings) == 0 & n.overlap > 0) {
  
    # if both traits have scored the same taxa the same keep older/newer?
    if (all(scores1[cc] == scores2[cc], na.rm=TRUE)) {   
      newscores <- scores[,1]
      drops[i] <- id[1]
      message('Scores identical; dropping character ', id[1])
    }
  
    # if both traits have scored the same taxa, but scores differ either keep (if just reversed scores) or leave alone
    if (all(scores1[cc] == scores2rev[cc], na.rm=TRUE)) {
      newscores <- scores1
      drops[i] <- id[1]
      message('Scores differ systematically (reversed values); dropping trait ', id[1], '; check state labels to confirm')
    }
    if (!all(scores1[cc] == scores2[cc]) & !all(scores1[cc] == scores2rev[cc])) {
      warning('Traits differ in their scorings. Keeping both.')
    }
  }

  # ind <- which(is.na(scores), arr.ind=TRUE)  # maybe useful at some point??
  
  # Case 2 - Different number of taxa scores, some overlapping scores
  if (diff(scorings) != 0 & n.overlap > 0) {
    # if trait 1 scores more taxa, overlapping scores same - keep character with more scorings
    if (all(scores1[cc] == scores2[cc], na.rm=TRUE)) {
      drops[i] <- id[which.min(scorings)]  # drop trait with fewer scorings
      message('Score overlap identical; keeping character with more scorings and dropping character ', id[which.min(scorings)])
    }

    # if trait 1 scores more taxa, overlapping scores different - if reversed keep:
    if (all(scores1[cc] == scores2rev[cc], na.rm=TRUE)) {
      drops[i] <- id[which.min(scorings)]
      warning('Assuming characters are equivalent (reversed scorings) and dropping character with fewer scorings (', id[2], '); check state labels to confirm')
    }

    # otherwise do nothing:
    if (!all(scores[cc,1] == scores[cc,2]) & !all(scores[cc,1] == rev(scores[cc,2]))) {
      warning('Traits differ in their scorings; keeping both')
    }
  }

  # Case 3 - trait 1 scores non-overlapping with trait 2 scores
  
  # merge characters:
  if (n.overlap == 0) {
    newscores <- pmin(scores[,1], scores[,2], na.rm=TRUE)
    res$data[, id[2]] <- newscores  # replace values
    drops[i] <- id[1]  # drop trait with fewer scorings
    warning('Merging non-overlapping character scorings; dropping character ', id[1], '; check state labels to confirm')
  }  
}

# final drop of characters:
res$data <- res$data[, -na.omit(drops)]
# drop character labels:
res$charlabels <- res$charlabels[-na.omit(drops)]
# drop state labels:
res$statelabels <- res$statelabels[-na.omit(drops)]

# TODO drop character sets, charpartitions,etc...
# ....

res$charset <- res$charset[-na.omit(drops)]
res$charnums <- res$charnums[-na.omit(drops)]
res$charpartition <- res$charpartition[-na.omit(drops)]

res$dups <- t(dups)

return(res)

}
