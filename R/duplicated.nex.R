#' Find duplicate or overlapping characters in nexus files
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
#' @param weighting vector for parts of character (before comma, after comma, states)
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
#' TODO write a lda.nex() function? separate the dup finding and dropping? maybe filter.nex()?
#' TODO add option to "nest" state labels in the network graph (so things like "size of distal end" wouldn't be matched across all characters, only those with state label as, say, "Humerus...")
#' TODO be able to plot "subclusters" (looking at all connections among characters, keeping only terms in common between the two)
#'

# x=twig1
# train=FALSE
# cutoff=0.15
# method="jw"
# drop=FALSE
# commasep=TRUE

duplicated.nex <- function(x, map = NULL, method = c("terms", "jw", "cosine"),
  within_dataset = FALSE, commasep = FALSE, parts = FALSE, force = FALSE, n = 25,
  train = FALSE, plot = FALSE,  cutoff = 0.35, drop = FALSE, weighting = c(1, 1, 1), K = 1, ...) {

  require(stringdist)
  require(tm)
  require(MASS)
  require(igraph)
  require(textmineR)

  method <- match.arg(method)

  ################################################################################
  # Check if user provides "map" of duplicated characters (e.g., list(1 = c(2, 3, 5), 5 = c(2, 4, 8)))
  # char1 duplicates 2,3,5; char 5 duplicates 2,4,8
  ################################################################################

  if (!is.null(map)) {
    if (class(map)=="list") {
      dups <- matrix(unlist(map), ncol = length(map))
    } else {
    dups <- t(map)
    }
  }

  # IF no map, then do automated discovery of duplicate characters:

  if (is.null(map)) {
    
    oldcharnames <- x$charlabels

    oldstatenames <- x$statelabels

  ################################################################################
  # Calculate text distances (method 1 - terms)
  ################################################################################

    if (method=="terms") {
      termlist <- cleantext(paste(oldcharnames, oldstatenames))
      # 2 options here- igraph/network modularity, or Ward clustering
      # weighting or not (weighting options= 'weightTf', 'weightTfIdf', 'weightBin', 'weightSMART'
      # make document matrix, sort by term frequency
      if (K == 1) {
        dtm <- TermDocumentMatrix(termlist, control = list(wordLengths = c(2, Inf), weighting = function(x) weightTfIdf(x, normalize = TRUE), stemming=FALSE))
        m <- as.matrix(dtm)
        g <- graph_from_incidence_matrix(m, weighted=TRUE)
        # TODO add option to specify different clustering algorithms (cluster_fast_greedy, cluster_label_prop, etc.)
        cl <- cluster_walktrap(g, weights=E(g)$weight)
        grps <- communities(cl)  # get clusters
      }
      if (K > 1) {
        dtm <- TermDocumentMatrix(termlist, control = list(wordLengths = c(2, Inf), stemming=FALSE))
        m <- as.matrix(dtm)
        tf_mat <- TermDocFreq(t(m))
        # TF-IDF and cosine similarity
        tfidf <- t(t(m)[ , tf_mat$term ]) * tf_mat$idf
        tfidf <- t(tfidf)
        csim <- tfidf / sqrt(rowSums(tfidf * tfidf))
        csim <- csim %*% t(csim)
        cdist <- as.dist(1 - csim)
        hc <- hclust(cdist, "ward.D")
        clust <- cutree(hc, K)
        grps <- split(names(clust), clust)
      }
      # make all possible combinations within groups/clusters
      dups <- lapply(lapply(grps, str_extract, "\\d+"), na.omit)
      dups <- lapply(dups, as.numeric)
      keep <- sapply(dups, length) > 1
      dups <- dups[keep]
      grps <- grps[keep]
      dups <- t(do.call(cbind, lapply(dups, combn, m=2)))
      # TODO need to output dups list based on some cutoff??
      stringdists.output <- NA
    }
    # text similarity using a few different methods (whole character statement, broken up character statement -locator, variable, states)

    # TODO lapply this stuff, e.g.-  list(part1, part2, part3)

  ################################################################################
  # Calculate text distances (method 2 - fuzzy string distance)
  ################################################################################

    if (method %in% c("jw", "cosine")) {
      
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
      
      # Clean up character names
      part1 <- sapply(cleantext(part1), as.character)
      part2 <- sapply(cleantext(part2), as.character)
      part3 <- sapply(cleantext(part3), as.character)


################################################################################
# Setup term weightings
################################################################################

      file <- x$file

      # only comparisons BETWEEN datasets/character types, not within:
      if (within_dataset) {  
        ids <- 1:ncol(x$data)
        splits <- split(ids, file)
      }

      # only look within character partitions
      if (!is.null(x$charpartition) & length(unique(x$charpartition)) > 1) {
        charpart <- x$charpartition
        # id <- charpart[pairids[1, ]] == charpart[pairids[2, ]]
        # pairids <- pairids[, id]
        split(ids, charpart)
      }

      sdist <- lapply(seq_along(splits), function(z) {
        nms <- splits[[z]]
        sd1 <- as.matrix(stringdistmatrix(part1[splits[[z]]], method=method))
        sd2 <- as.matrix(stringdistmatrix(part2[splits[[z]]], method=method))
        sd3 <- as.matrix(stringdistmatrix(part3[splits[[z]]], method=method))
        dimnames(sd1) <- dimnames(sd2) <- dimnames(sd3) <- list(nms, nms)
        list(sd1, sd2, sd3)
      })

# weighting <- c(1,1,1)

      # calculate final distances
      sdist_final <- lapply(seq_along(sdist), function(z) {
        # z=1
        d1 <- weighting[1] * sdist[[z]][[1]]
        d2 <- weighting[2] * sdist[[z]][[2]]
        d3 <- weighting[3] * sdist[[z]][[3]]
        res <- (d1+d2+d3)/3
        diag(res) <- NA
        res[upper.tri(res)] <- NA
        res
      })

      dups <- lapply(seq_along(sdist_final), function(z) {
        sset <- which(sdist_final[[z]] < cutoff, arr.ind=TRUE)
        data.frame(char1 = rownames(sdist_final[[z]])[sset[, 1]],
                   char2 = colnames(sdist_final[[z]])[sset[, 2]],
                   stringdist = sdist_final[[z]][which(sdist_final[[z]] < cutoff)])
      })

      dups <- do.call(rbind, dups)
      dups <- dups[order(dups$stringdist), ]
      dups$char1 <- as.numeric(as.character(dups$char1))
      dups$char2 <- as.numeric(as.character(dups$char2))

      stringdists <- unlist(sapply(seq_along(sdist_final), function(z) as.numeric(as.dist(sdist_final[[z]]))))

      names(stringdists) <- 1:length(stringdists)
      
      pairids <- t(combn(1:ncol(x$data), m=2))
      pairids <- lapply(splits, combn, m=2)
      pairids <- t(do.call(cbind, pairids))
      pairids <- cbind(pairids, stringdist=stringdists)

      sset.dist <- stringdists[stringdists < cutoff]

    }

  ################################################################################
  # Training to identify duplicates
  ################################################################################
  if (train) {
    
    # sscut <- setNames(dups$stringdist, 1:nrow(dups))
    names(stringdists) <- seq_along(stringdists)

    sscut <- stringdists[stringdists < cutoff]

    if (length(sscut)==0) {
      stop("No matches found, try increasing cutoff")
    }

    # sample evenly over range of string distances
    # n <- 25
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
    printpair <- function(i) {
      id1 <- pairids[i, 1]
      id2 <- pairids[i, 2]
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
      sens <- 0.25  
      matchid <- which(pred$posterior[, 2] > sens)  # matches with 'sens' (i.e. 50%) probability of being a match
      dups <- pairids[matchid, ]
    }
    if (plot) {
      plot(predict(lda1, df)$posterior[,2] ~ df[,2], type='b', col=df[,1]+1,
        pch=16, xlab="String distance", ylab="P(duplicate)")
      legend("topright", pch=16, col=c("black", "red"),
        legend=c("not dup", "dup"), bty="n")
      abline(h=sens, lty=2)
      title("LDA training results")  
      }
    
    if (train & length(matchid)!=0) {
      stringdists.output <- stringdists[matchid]
    }
    
  }

  if (!train & method!="terms") {
      sset <- as.numeric(names(which(sset.dist < cutoff)))
      if (length(sset) == 0) {
        warning('No matches found')
      }
      dups <- as.matrix(pairids[sset, ])
      stringdists.output <- stringdists[sset]
  }
}

  # resultant NEXUS file for outputting later

  res <- x

  ################################################################################
  # Do the dropping and merging characters:
  ################################################################################

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
      res$dups <- data.frame("char1" = dups[, 1], "charnum1" = x$charnum[dups[, 1]],
                             "char2" = dups[, 2], "charnum2" = x$charnum[dups[, 2]],
                             stringdist = stringdists.output)
      res$dups <- res$dups[order(res$dups$stringdist), ]
    }

  if (method=="terms") {
    res$clusters <- grps
  }

  if (method %in% c('jw', 'cosine')) {
    g <- graph_from_edgelist(t(apply(dups[,1:2], 1, as.character)), directed=F)
    grps <- communities(cluster_walktrap(g))
    res$clusters <- grps
  }
  
  res

}
