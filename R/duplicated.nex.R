#' Find duplicate traits in nexus files
#' Function uses fuzzy text matching to output a list of potentially redundant
#' characters in a `nex` object, waits for user input to confirm, and then
#' handles merging of characters and outputs a new `nex` object
#' @param x (required) a nexus data object
#' @param cutoff output potential matches under this value for the string distance between characters
#' @param map list of numeric vectors with original character ids mapped to their new ones
#' @param force whether to force character dropping and avoid data checks (not advised)
#' @param n number of pairs to use for training discriminant function analysis to predict matches
#' @param within_dataset whether to limit search to only among-dataset characters (e.g., useful if you are certain the individual matrices do not contain duplicates)
#' @param train whether to use training algorithm to predict duplicate matches
#' @param plot whether to plot the XX vs XX curve
#' @param weighted whether to use term weighting when searching for similar characters (more common downweighted)
#' @param parts whether to use weighting of locators/states or not (i.e. whole statement used in fuzzy matching)
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' @examples \dontrun{
#' x <- read.nex(file='example/toy1.nex')
#' map <- list(c(2, 4), c(156, 157))
#' duplicated.nex(x, method = 'user', map = map)
#' duplicated.nex(x, method = 'automated', map = map)
#' duplicated(x, cutoff = 0.01, method = 'auto')
#' duplicated(x, method = 'user', map=list(c(2,3)))
#' @author Chad Eliason \email{chad_eliason@@utexas.edu}
#'
#' TODO list
#' [x] should have this output the string distances as well in the dup data frame
#' [x] exclude positional terms during matching
#' [x] maybe use a different string distance metric (e.g., 'vomer' AND 'vomer extending laterally' give dist = 0)
#'
duplicated.nex <- function(x, map = NULL, force = FALSE, n = 25, train = FALSE,
  cutoff = 0.35, method = c("jw", "cosine", "terms"), within_dataset = FALSE, plot = FALSE,
  drop = FALSE, commasep = FALSE, weighted = FALSE, parts = FALSE, ...) {

  require(stringdist)
  require(tm)
  require(MASS)
  require(igraph)

  method <- match.arg(method)

  ################################################################################
  # Check if user provides "map" of duplicated characters (e.g., list(1 = c(2, 3, 5), 5 = c(2, 4, 8)))
  # char1 duplicates 2,3,5; char 5 duplicates 2,4,8
  ################################################################################

  if (!is.null(map)) {
    if (class(map)=="list") {
      dups <- matrix(unlist(map), ncol = length(map))  
    } else {
    # pairids <- matrix(unlist(map), ncol=length(map))
    dups <- t(map)
    }
  }

  # IF no map, then do automated discovery of duplicate characters:

  ################################################################################
  # Calculate text distances (method 1)
  ################################################################################

  if (is.null(map)) {
    oldcharnames <- x$charlabels
    oldstatenames <- x$statelabels
    # remove comments in square brackets:
    oldcharnames <- str_replace_all(oldcharnames, "[\\s]*[\\[\\(].+[\\]\\)][\\s]*", "")
    oldstatenames <- str_replace_all(oldstatenames, "[\\s]*[\\[\\(].+[\\]\\)][\\s]*", "")
    # remove comments in after 'Note:'
    oldcharnames <- str_replace_all(oldcharnames, regex("Note\\:.*?$", ignore_case=TRUE), "")
    oldstatenames <- str_replace_all(oldstatenames, regex("Note\\:.*?$", ignore_case=TRUE), "")
    # part3 <- str_replace_all(twig1$charlab[91], "Note\\:.*?$", "")

    if (method=="terms") {
      tocut <- c("with", "than", "then", "those", "with", "to", "the", "and", "an", "a", "or", "of", "for", "not", "along", "length", "less", "below", "above", "around", "longer", "shorter", "sits", "absent", "present", "form", "process", "state", "view", "margin", "shape", "placed", "recess", "(UN)?ORDERED", "(un)?ordered", "lateral", "distal", "ventral", "posterior", "anterior", "medial", "dors", "external")
      termlist <- paste(oldcharnames, oldstatenames)
      termlist <- Corpus(VectorSource(termlist))
      # Convert the text to lower case
      termlist <- tm_map(termlist, content_transformer(tolower))
      # Remove numbers
      termlist <- tm_map(termlist, removeNumbers)
      # Remove english common stopwords
      termlist <- tm_map(termlist, removeWords, stopwords("english"))
      # Remove your own stop word
      # specify your stopwords as a character vector
      termlist <- tm_map(termlist, removeWords, tocut)
      # Remove punctuations
      termlist <- tm_map(termlist, removePunctuation)
      # Eliminate extra white spaces
      termlist <- tm_map(termlist, stripWhitespace)
      # Latinized tokens
      termlist <- tm_map(termlist, schinke)
      # weighting or not (weighting options= 'weightTf', 'weightTfIdf', 'weightBin', 'weightSMART'
      dtm <- TermDocumentMatrix(termlist, control = list(wordLengths = c(2, Inf), weighting = function(x) weightTfIdf(x, normalize = TRUE), stemming=TRUE))
      # dtm <- TermDocumentMatrix(termlist, control = list(wordLengths = c(2, Inf)))
      # dtm <- TermDocumentMatrix(termlist, control = list(wordLengths = c(3, Inf), weighting = 'weightTfIdf'))
      # make document matrix, sort by term frequency
      m <- as.matrix(dtm)
      v <- sort(rowSums(m),decreasing=TRUE)
      d <- data.frame(word = names(v), freq = v)
      g <- graph_from_incidence_matrix(m, weighted = TRUE)
      # TODO add option to specify different clustering algorithms
      # cl <- cluster_walktrap(g, weights = E(g)$weight)
      cl <- cluster_fast_greedy(g, weights = E(g)$weight)
      # cl <- cluster_edge_betweenness(g, weights = E(g)$weight)
      # get clusters
      grps <- communities(cl)
      # make all possible combinations within groups/clusters
      dups <- sapply(sapply(grps, str_extract, "\\d+"), na.omit)
      keep <- sapply(dups, length) > 1
      dups <- dups[keep]
      grps <- grps[keep]
      dups <- do.call(cbind, sapply(dups, combn, m=2))
      # plot (for fun)
      # par(mar=c(0,0,0,0))
      # plot(g, vertex.label=NA, vertex.size=1, edge.width=2*E(g)$weight, edge.color=rgb(1,0,0,0.1), mark.col=rgb(0,0,1,.3), mark.border="black")
      # ggplot(d[1:100,], aes(x=freq, y=reorder(word, freq))) + geom_point()
      # TODO need to output dups list based on some cutoff??
      stringdists.output <- NA
    }

    # text similarity (a few methods)
    # whole character statement
    # broken up character statement (locator, variable, states)
    if (method %in% c("jw", "cosine")) {
      # remove stop words
      oldcharnames <- removeWords(oldcharnames, stopwords("en"))
      oldstatenames <- removeWords(oldstatenames, stopwords("en"))
      # comma separated?
      # commasep <- FALSE
      if (commasep) {
        matches <- str_match(oldcharnames, "^((.*?),)?(.*)")
        part1 <- matches[, 3]
        part2 <- matches[, 4]
        id <- is.na(part1)
        part1[id] <- part2[id]
        part2[id] <- NA
        part3 <- oldstatenames
      } else {
        part1 <- oldcharnames
        part2 <- NA
        part3 <- oldstatenames
      }
      # remove comments in square brackets:
      part3 <- str_replace_all(part3, "[\\s]*[\\[\\(].+[\\]\\)][\\s]*", "")
      # remove comments in after 'Note:'
      part3 <- str_replace_all(part3, "Note\\:.*?$", "")
      # remove reference to figures
      part1 <- str_replace(part1, "\\([Ff]ig(\\.|ure).*\\)", "")
      part2 <- str_replace(part2, "\\([Ff]ig(\\.|ure).*\\)", "")
      part3 <- str_replace(part3, "\\([Ff]ig(\\.|ure).*\\)", "")
      # remove short words
      tocut <- c("with", "than", "then", "those", "with", "to", "the", "and", "an", "a", "or", "of", "for", "not")
      tocut <- paste0("\\b", tocut, "\\b", collapse="|")
      part1 <- gsub("[[:punct:]]", " ", part1)
      part1 <- gsub("[[:digit:]]", " ", part1)
      part1 <- gsub(tocut, "", part1, ignore.case=TRUE)
      part2 <- gsub("[[:punct:]]", " ", part2)
      part2 <- gsub("[[:digit:]]", " ", part2)
      part2 <- gsub(tocut, "", part2, ignore.case=TRUE)
      part3 <- gsub("[[:punct:]]", " ", part3)
      part3 <- gsub("[[:digit:]]", " ", part3)
      part3 <- gsub(tocut, "", part3, ignore.case=TRUE)
      # Removing non-informative terms
      # [ ] fix it so this will remove things like "posteroventral"
      tocut <- c("along", "less", "below", "above", "around", "longer", "shorter", "sits", "absent", "present", "form", "process", "state", "view", "margin", "shape", "placed", "recess", "(UN)?ORDERED", "(un)?ordered")
      tocut <- paste0("\\b", tocut, "\\b", collapse="|")
      part1 <- gsub(tocut, "", part1, ignore.case=TRUE)
      part2 <- gsub(tocut, "", part2, ignore.case=TRUE)
      part3 <- gsub(tocut, "", part3, ignore.case=TRUE)
      # Removing positional terms
      tocut <- c("lateral", "distal", "ventral", "posterior", "anterior", "medial", "dors", "external")
      tocut <- paste0("\\b", tocut, "(\\w*)?", collapse="|")
      part1 <- gsub(tocut, "", part1, ignore.case=TRUE)
      part2 <- gsub(tocut, "", part2, ignore.case=TRUE)
      part3 <- gsub(tocut, "", part3, ignore.case=TRUE)
      # Clean up character names
      part1 <- cleantext(part1)
      part2 <- cleantext(part2)
      part3 <- cleantext(part3)
      part1[part1 == ""] <- NA
      part2[part2 == ""] <- NA
      part3[part3 == ""] <- NA
      part1 <- sapply(part1, unique)
      part2 <- sapply(part2, unique)
      part3 <- sapply(part3, unique)

    ################################################################################
    # Setup term weightings
    ################################################################################

      # [x] weight words by usage in dataset (common words, not important; rare words, more important)
      # [x] try weighting by commenness of word in the matrix

      weighting <- cbind(rep(1, length(part1)), rep(1, length(part2)), rep(1, length(part3)))

      # DONE: weight first few words (if comma) more heavily
      # e.g., Dorsal vertebrae, shape -- vs Caudal vertebrae, shape (shouldn't be identified as duplicates)
      id <- !is.na(part2)
      weighting[id, 2] <- 0.5
      weighting[is.na(part2), 2] <- 0

      # DONE: if state labels are absent/present (i.e. non-informative), then remove from duplicate matching
      # this is already being done above
      # id <- grep("('present'\\s*'absent')|(present\\s*absent)|('absent'\\s*'present')|(absent\\s*present)", part3)  
      # weighting[id, 3] <- 0

      # generate all possible pairs of character combinations
      pairids <- combn(seq_along(part1), m=2)

      # only comparisons BETWEEN datasets/character types, not within:
    
      file <- x$file

      if (within_dataset) {  
        id <- file[pairids[1, ]] != file[pairids[2, ]]
        pairids <- pairids[, id]
      }

      # only look within character partitions
      if (!is.null(x$charpartition) & length(unique(x$charpartition)) > 1) {
        charpart <- x$charpartition
        id <- charpart[pairids[1, ]] == charpart[pairids[2, ]]
        pairids <- pairids[, id]
      }

      stringdists <- numeric(length = ncol(pairids))
      
      # create progress bar
      pb <- txtProgressBar(min = 0, max = ncol(pairids), style = 3)

    ################################################################################
    # Calculate text distances (takes ~200 seconds for 2.2M comparisons):
    ################################################################################
    
      if (parts) {
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
          sd1 <- (1/wt1) * stringdist(str11, str21, method=method, ...)
          sd2 <- (1/wt2) * stringdist(str12, str22, method=method, ...)
          sd3 <- (1/wt3) * ifelse(any(weighting[c(id1, id2), 3]==0), NA, stringdist(str13, str23, method=method, ...))
          stringdists[i] <- sum(sd1, sd2, sd3, na.rm=TRUE)
          setTxtProgressBar(pb, i)  # update progress bar
        }
      } else {
        for (i in 1:ncol(pairids)) {  
          id1 <- pairids[1, i]
          id2 <- pairids[2, i]
          newcharnames <- paste(na.omit(part1, part2))
          newstatenames <- paste(na.omit(part3))
          str1a <- newcharnames[pairids[1,i]]
          str1b <- newstatenames[pairids[1,i]]
          str2a <- newcharnames[pairids[2,i]]
          str2b <- newstatenames[pairids[2,i]]
          stringdist1 <- stringdist(str1a, str2a, method = method)
          stringdist2 <- stringdist(str1b, str2b, method = method)
          # stringdists[i] <- stringdist1 + stringdist2
          stringdists[i] <- sum(stringdist1, stringdist2, na.rm=TRUE)
          setTxtProgressBar(pb, i)  # update progress bar
        }
      }

      close(pb)

      names(stringdists) <- 1:length(stringdists)

      sset.dist <- sort(stringdists)

      sset <- as.numeric(names(stringdists))

    }

  ################################################################################
  # Training to identify duplicates
  ################################################################################

  # dups doesn't exist yet
  # if (dim(dups)[2]==0) {
  #   stop("No matches found")
  # }

  if (train) {
    sscut <- stringdists[stringdists < cutoff]
    if (length(sscut)==0) {
      stop("No matches found, try increasing cutoff")
    }
    sscut <- sort(sscut)
    # sample evenly over range of string distances
    # [x] employ a random sampling approach
    ss <- sapply(seq(0, cutoff, length=n), function(z) {
      which.min(abs(sscut - z))
    })
    if (length(ss) < n) {
      stop("Number of pairs to assess is less than specified, try increasing cutoff")
    }
    ss <- unique(ss)
    sstrain <- sscut[ss]
    sset.dist <- sstrain
    sset <- as.numeric(names(sset.dist))
    
    # run loop to determine matches
    
    printpair <- function(sset) {
      id1 <- pairids[1, sset]
      id2 <- pairids[2, sset]
      # print pairs of characters:
      cat('\n-------\nTrait pair', i, ' (string distance = ', sset.dist[as.character(sset)], ') \n\nCHARLABELS:\n\n',
        oldcharnames[id1], ' (', x$file[id1], ', character ', x$charnums[id1],')\n\n',
        oldcharnames[id2], ' (', x$file[id2], ', character ', x$charnums[id2],')\n\nSTATELABELS:\n\n',
        x$statelabels[id1], '\n\n',
        x$statelabels[id2], '\n----------\n\n', sep="")
      cat('Are these the same traits (y/n/N/q)\n')
    }

    answer <- rep("", length(sset))

    for (i in 1:length(sset)) {  
      printpair(sset[i])  # print pair of characters
      ans <- scan(n = 1, what = 'character')  # wait for user input
      if (ans=='q') {
        stop('Function terminated by user')
      }
      if (ans=='N') {
        answer[i:length(answer)] <- answer[i-1]
        break
      }
      answer[i] <- ans
    }
    
    # Linear discrimant analysis
    df <- data.frame(match = ifelse(answer=="y", 1, 0), dist = sset.dist)
    if (all(diff(df$match)==0)) {
      stop("All pairs identified as matches/nonmatches, need variability for LDA training analysis")
    } else {
      lda1 <- lda(match~dist, data=df)
      pred <- predict(lda1, data.frame(dist=stringdists))
      sens <- 0.5  
      matchid <- which(pred$posterior[, 2] > sens)  # matches with 'sens' (i.e. 50%) probability of being a match
      dups <- pairids[, matchid]
    }
    if (plot) {
      plot(predict(lda1, df)$posterior[,2]~df[,2], type='b', col=df[,1]+1, pch=16, xlab="String distance", ylab="P(duplicate)")
      legend("topright", pch=16, col=c("black", "red"), legend=c("not dup", "dup"), bty="n")
      abline(h=0.5, lty=2)
      title("LDA training results")  
      }
    if (train & length(matchid)!=0)
      stringdists.output <- stringdists[matchid]
    }
    }

  if (!train & method!="terms") {
      sset <- as.numeric(names(which(sset.dist < cutoff)))
      if (length(sset) == 0) {
        warning('No matches found')
      }
      dups <- as.matrix(pairids[, sset])
      stringdists.output <- stringdists[sset]
  }

  # if (method=="terms" | length(matchid)==0) {
  #   stringdists.output <- NA
  # }


  # resultant NEXUS file for outputting later
  res <- x

  ################################################################################
  # Do the dropping and merging characters:
  ################################################################################

  # [x] need to have this output what characters are duplicated, retained, etc.

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

  ################################################################################
  # Output, drop duplicated characters, labels, etc.:
  ################################################################################

  res$data <- res$data[, drops2]
  res$charlabels <- res$charlabels[drops2]
  res$statelabels <- res$statelabels[drops2]
  res$charset <- res$charset[drops2]
  res$charnums <- res$charnums[drops2]
  res$charpartition <- res$charpartition[drops2]

  if (is.null(dups)) {
    dups <- matrix(NA, nrow=2, ncol=1)
  } else {
    dups <- apply(dups, 2, as.numeric)  
  }

  if (sum(dups)==0) {
    res$dups <- NULL
  } else {
      res$dups <- data.frame("char1" = dups[1, ], "charnum1" = x$charnum[dups[1, ]],
                             "char2" = dups[2, ], "charnum2" = x$charnum[dups[2, ]],
                             stringdist = stringdists.output)
    }

  if (method=="terms") {
    res$clusters <- grps
  }

  res

}



# testing
# 
# xx <- duplicated.nex(twig1[,1:100], train=TRUE, cutoff=0.55, method='jw', n = 10)
# duptree(xx)

