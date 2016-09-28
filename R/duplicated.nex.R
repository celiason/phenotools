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
#' @param n number of pairs to use for training discriminant function analysis to predict matches
#' @param train whether to use training algorithm to predict duplicate matches
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

# @done should have this output the string distances as well in the dup data frame

# @done exclude positional terms during matching

# Testing zone:
# x <- read.nex("/Users/chadeliason/Documents/UT/projects/phenome/data/clarke_2002.nex")
# map <- NULL
# train = F
# cutoff=.25
# n = 25
# method = "jw"
# within_dataset=FALSE

duplicated.nex <- function(x, map = NULL, force = FALSE, n = 25, train = TRUE,
  cutoff = 0.35, method = c("jw", "cosine"), within_dataset = FALSE, plot = TRUE,
  drop = FALSE) {

library(stringdist)
library(tm)
library(MASS)
method <- match.arg(method)

matchfun <- function(sset) {
  id1 <- pairids[1, sset]
  id2 <- pairids[2, sset]
  # print pairs of characters:
  cat('\n-------\nTrait pair', i, ' (string distance = ', sset.dist[as.character(sset)], ') \n\nCHARLABELS:\n\n',
    oldcharnames[id1], ' (', file[id1], ', character ', x$charnums[id1],')\n\n',
    oldcharnames[id2], ' (', file[id2], ', character ', x$charnums[id2],')\n\nSTATELABELS:\n\n',
    x$statelabels[id1], '\n\n',
    x$statelabels[id2], '\n----------\n\n', sep="")
  cat('Are these the same traits (y/n/q)\n')
  # wait for user input
  answer <- scan(n = 1, what = 'character')
  if (answer=='q') {
    stop('Function terminated by user')
  }
  answer
}

# If provide map of duplicated characters
if (!is.null(map)) {
  if (class(map)=="list") {
    dups <- matrix(unlist(map), ncol = length(map))  
  } else {
  # pairids <- matrix(unlist(map), ncol=length(map))
  dups <- t(map)
  }
}

# IF no map, then do automated discovery of duplicate characters:

# Calculate text distances
if (is.null(map)) {
  oldcharnames <- x$charlabels
  oldstatenames <- x$statelabels

  # @done break apart character names into chunks (before comma, after comma)
  matches <- str_match(oldcharnames, "^((.*),)?(.+)")
  part1 <- matches[, 3]
  part2 <- matches[, 4]
  id <- is.na(part1)
  part1[id] <- part2[id]
  part2[id] <- NA
  part3 <- oldstatenames

  # newstatenames <- oldstatenames
  # newcharnames <- oldcharnames
  # remove short words
  tocut <- c("to", "the", "and", "an", "a", "or", "of")
  tocut <- paste0("\\b", tocut, "\\b", collapse="|")
  # newcharnames <- gsub(tocut, "", newcharnames)
  part1 <- gsub(tocut, "", part1)
  part2 <- gsub(tocut, "", part2)
  part3 <- gsub(tocut, "", part3)

  # remove vague terms
  tocut <- c("form", "process", "state", "view", "margin", "shape", "placed")
  tocut <- paste0("\\b", tocut, "\\b", collapse="|")
  # newcharnames <- gsub(tocut, "", newcharnames)
  part1 <- gsub(tocut, "", part1)
  part2 <- gsub(tocut, "", part2)
  part3 <- gsub(tocut, "", part3)

  # remove positional terms
  tocut <- c("lateral", "distal", "ventral", "posterior", "anterior", "medial", "dors", "external")
  tocut <- paste0("\\b", tocut, "(\\w*)?", collapse="|")
  # newcharnames <- gsub(tocut, "", newcharnames)
  part1 <- gsub(tocut, "", part1)
  part2 <- gsub(tocut, "", part2)
  part3 <- gsub(tocut, "", part3)
  
  # clean up character names:
  # newstatenames <- cleantext(newstatenames)
  # newcharnames <- cleantext(newcharnames)
  part1 <- cleantext(part1)
  part2 <- cleantext(part2)
  part3 <- cleantext(part3)

  # Setup weightings

  weighting <- cbind(rep(1, length(part1)), rep(1, length(part2)), rep(1, length(part3)))

  # @done Weight first few words (if comma) more heavily
  # e.g., Dorsal vertebrae, shape -- vs Caudal vertebrae, shape (shouldn't be identified as duplicates)
  id <- !is.na(part2)
  weighting[id, 2] <- 0.5
  weighting[is.na(part2), 2] <- 0

  # @todo downweight based on commonness of term in dataset


  # @done if state labels are absent/present (i.e. non-informative), then remove from duplicate matching
  id <- grep("('present'\\s*'absent')|(present\\s*absent)|('absent'\\s*'present')|(absent\\s*present)", newstatenames)
  weighting[id, 3] <- 0

  # generate all possible pairs of character combinations
  pairids <- combn(seq_along(part1), m=2)
  # only want comparisons BETWEEN datasets/character types, not within:
  file <- x$file
  if (within_dataset) {  
    id <- file[pairids[1, ]] != file[pairids[2, ]]
    pairids <- pairids[, id]
  }
  # only look within character partitions:
  if (!is.null(x$charpartition) & length(unique(x$charpartition)) > 1) {
    charpart <- x$charpartition
    id <- charpart[pairids[1, ]] == charpart[pairids[2, ]]
    pairids <- pairids[, id]
  }
  stringdists <- numeric(length = ncol(pairids))
  # create progress bar
  pb <- txtProgressBar(min = 0, max = ncol(pairids), style = 3)

  # calculate text distances (takes ~200 seconds for 2.2M comparisons):
  # @TODO: WORK ON OPTIMIZING THIS PART
  for (i in 1:ncol(pairids)) {
    id1 <- pairids[1, i]
    id2 <- pairids[2, i]
    str11 <- part1[id1]
    str12 <- part2[id1]
    str13 <- part3[id1]
    str21 <- part1[id2]
    str22 <- part2[id2]
    str23 <- part3[id2]
    wt1 <- sum(weighting[c(id1, id2), 1])
    wt2 <- sum(weighting[c(id1, id2), 2])
    wt3 <- sum(weighting[c(id1, id2), 3])
    sd1 <- (1/wt1) * stringdist(str11, str21, method=method)
    sd2 <- (1/wt2) * stringdist(str12, str22, method=method)
    sd3 <- (1/wt3) * ifelse(any(weighting[c(id1, id2), 3]==0), NA, stringdist(str13, str23, method=method))
    stringdists[i] <- sum(sd1, sd2, sd3, na.rm=TRUE)
    # str1a <- newcharnames[pairids[1,i]]
    # str1b <- newstatenames[pairids[1,i]]
    # str2a <- newcharnames[pairids[2,i]]
    # str2b <- newstatenames[pairids[2,i]]
    # stringdist1 <- stringdist(str1a, str2a, method = method)
    # stringdist2 <- stringdist(str1b, str2b, method = method)
    # stringdists[i] <- stringdist1 + stringdist2
    setTxtProgressBar(pb, i)  # update progress bar
  }
  close(pb)
  names(stringdists) <- 1:length(stringdists)
  sset.dist <- sort(stringdists)
  sset <- as.numeric(names(stringdists))
}



# Training to identify duplicates
if (train) {
  sscut <- stringdists[stringdists < cutoff]
  sscut <- sort(sscut)
  # sample evenly over range of string distances
  # @done employ a random sampling approach?
  ss <- sapply(seq(0, cutoff, length=n), function(z) {
    which.min(abs(sscut - z))
  })
  if (length(ss) < n) {
    warning("Number of pairs to assess is less than specified. Try increasing cutoff.")
  }
  ss <- unique(ss)
  sstrain <- sscut[ss]
  sset.dist <- sstrain
  sset <- as.numeric(names(sset.dist))
  # run loop to determine matches
  answer <- sapply(sset, matchfun)
  # linear discrimant analysis
  df <- data.frame(match=ifelse(answer=="y", 1, 0), dist=sset.dist)
  lda1 <- lda(match~dist, data=df)
  pred <- predict(lda1, data.frame(dist=stringdists))
  # matches with > 50% probability of being a match
  matchid <- which(pred$posterior[, 2] > 0.5)
  # matchid <- as.numeric(names(pred$posterior[, 2][pred$posterior[, 2] > 0.5]))
  dups <- pairids[, matchid]

  # check this?
  stringdists.output <- stringdists[, matchid]

  if (plot) {
    plot(predict(lda1, df)$posterior[,2]~df[,2], type='b', col=df[,1]+1, pch=16, xlab="String distance", ylab="P(duplicate)")
    legend("topright", pch=16, col=c("black", "red"), legend=c("not dup", "dup"), bty="n")
    # predict(lda1, df)$posterior[,2]>0.5
    abline(h=0.5, lty=2)
    title("LDA training results")  
  }
}

if (!train) {
  sset <- as.numeric(names(which(sset.dist < cutoff)))
  if (length(sset) == 0) {
    stop('No putative matches found.')
  }
  dups <- pairids[, sset]
  stringdists.output <- stringdists[sset]
}

  # answer <- character(length = length(sset))

  # @TODO: maybe use a different string distance metric (e.g., 'vomer' AND 'vomer extending laterally' give dist = 0)

  # loop over pairs, wait for user input confirming redundant traits
  # sset = integers corresponding to row of pairids as possible matches
  # sset.dist = string distances

# Dropping and merging characters:
# [x] need to have this output what characters are duplicated, retained, etc.
res <- x
drops2 <- rep(TRUE, ncol(res$data))
if (drop) {
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
    # Case 1: Same number of taxa scored
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
    # Case 2: Different number of taxa scores, some overlapping scores
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
    # Case 3: trait 1 scores non-overlapping with trait 2 scores; merge characters
    if (n.overlap == 0) {
      newscores <- pmin(scores[,1], scores[,2], na.rm=TRUE)
      res$data[, id[2]] <- newscores  # replace values
      drops[i] <- id[1]  # drop trait with fewer scorings
      warning('Merging non-overlapping character scorings; dropping character ', id[1], '; check state labels to confirm')
    }  
  }
drops2[drops] <- FALSE
}

# Output, drop duplicated characters, labels, etc.:
res$data <- res$data[, drops2]
res$charlabels <- res$charlabels[drops2]
res$statelabels <- res$statelabels[drops2]
res$charset <- res$charset[drops2]
res$charnums <- res$charnums[drops2]
res$charpartition <- res$charpartition[drops2]
res$dups <- data.frame("char1" = dups[1, ], "charnum1" = x$charnum[dups[1, ]],
  "char2" = dups[2, ], "charnum2" = x$charnum[dups[2, ]], stringdist = stringdists.output)
return(res)
}




# Test dataset
# devtools::load_all("~/github/nexustools")
# x <- read.nex("/Users/chadeliason/Documents/UT/projects/phenome/data/clarke_2002.nex")
# x.dup <- duplicated(x, cutoff=0.5, train=TRUE)
# head(x.dup$dups)
# tail(x.dup$dups)
# x$charlab[c(26, 27)]
# x$charlab[c(57, 59)]


