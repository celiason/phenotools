
################################################################################
# Training to identify duplicates
################################################################################
if (train) {
  names(stringdists) <- seq_along(stringdists)
  sscut <- stringdists[stringdists < cutoff]
  if (length(sscut)==0) {
    stop("No matches found, try increasing cutoff")
  }
  # sample evenly over range of string distances
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
  # function to print pairs of characters
  printpair <- function(i) {
    id1 <- pairids[i, 1]
    id2 <- pairids[i, 2]
    cat('\n-------\nTrait pair', i, ' (string distance = ', sset.dist[as.character(sset)], ') \n\nCHARLABELS:\n\n',
      oldcharnames[id1], ' (', x$file[id1], ', character ', x$charnums[id1],')\n\n',
      oldcharnames[id2], ' (', x$file[id2], ', character ', x$charnums[id2],')\n\nSTATELABELS:\n\n',
      x$statelabels[id1], '\n\n',
      x$statelabels[id2], '\n----------\n\n', sep="")
    cat('Are these the same traits (y/n/N/q)\n')
  }
  answer <- rep("", length(sset))
  # run loop to manually determine matches
  for (i in 1:length(sset)) {  
    printpair(sset[i])
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
