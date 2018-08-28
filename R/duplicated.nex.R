#' Find duplicate or overlapping characters in nexus files
#' 
#' Function uses fuzzy text matching to output a list of potentially redundant
#' characters in a `nex` object, waits for user input to confirm, and then
#' handles merging of characters and outputs a new `nex` object
#' 
#' @param x (required) a nexus data object
#' @param opt method to use for finding duplicates
#' @param method specific method arguments to pass to functions
#' @param n number of pairs to use for training discriminant function analysis to predict matches
#' @param within_dataset whether to limit search to only among-dataset characters (e.g., useful if you are certain the individual matrices do not contain duplicates)
#' @param weighting vector for parts of character (before comma, after comma, states)
#' @param cores how many cores to use (for parallel processing in traitcor option)
#' 
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' 
#' @examples \dontrun{
#' x <- read.nex(file='example/toy1.nex')
#' map <- list(c(2, 4), c(156, 157))
#' duplicated.nex(x, method = 'user', map = map)
#' duplicated.nex(x, method = 'automated', map = map)
#' duplicated(x, cutoff = 0.01, method = 'auto')
#' duplicated(x, method = 'user', map=list(c(2,3)))
#' 
#' @author Chad Eliason \email{celiason@@fieldmuseum.org}
#'
#' TODO write a lda.nex() function? separate the dup finding and dropping? maybe filter.nex()?
#' TODO add option to "nest" state labels in the network graph (so things like "size of distal end" wouldn't be matched across all characters, only those with state label as, say, "Humerus...")
#' TODO be able to plot "subclusters" (looking at all connections among characters, keeping only terms in common between the two)

# x <- dat2
# opt <- "comments"
# Error in `[.data.frame`(dups, , c("infile", "charnum", "author", "year",  : 
#   undefined columns selected

duplicated.nex <- function(x, opt=c("fuzzy", "terms", "comments", "traitcor"),
  method=NULL, within_dataset=FALSE, commasep=FALSE, weighting=c(1,1,1),
  K=1, cluster = c("infomap", "fast_greedy", "walktrap", "label_prop", "leading_eigen",
    "louvain", "optimal", "spinglass"), cores=1, ...) {

  require(stringdist)
  require(tm)
  require(MASS)
  require(igraph)
  require(textmineR)
  require(parallel)
  require(polycor)
  require(psych)
  require(pbapply)
  require(pbmcapply)

  cutoff <- Inf # I previously had this as an argument, but i think it's better in the printout function (that way all possible dups will be output and their string distances for later subsetting)

  opt <- match.arg(opt)
  
  cluster <- match.arg(cluster)  

  clustfun <- get(paste0("cluster_", cluster))
  
  # extract information from nexus file
  oldcharnames <- x$charlabels
  oldstatenames <- x$statelabels

  # Setup term weightings
  file <- x$file
  
  # split up character ids to use later
  if (within_dataset) {  
    ids <- 1:ncol(x$data)
    splits <- split(ids, file)
  } else {
    splits <- list(1:ncol(x$data))
  }

  # TODO add option to look within character partitions
  # if (!is.null(x$charpartition) & length(unique(x$charpartition)) > 1) {
  #   charpart <- x$charpartition
  #   splits <- split(ids, charpart)
  # }

  # Calculate text distances based on overlapping terms
  # 2 options here- igraph/network modularity, or Ward clustering
  # weighting or not (weighting options= 'weightTf', 'weightTfIdf', 'weightBin', 'weightSMART'
  if (opt == "terms") {
    chars <- paste(oldcharnames, oldstatenames)
    chars <- gsub("\\-\\s", "", chars)    
    termlist <- cleantext(chars)
    
    # make document matrix, sort by term frequency
    if (K == 1) {
      dtm <- TermDocumentMatrix(termlist, control = list(wordLengths = c(2, Inf), weighting = function(x) weightTfIdf(x, normalize = TRUE), stemming=FALSE))
      m <- as.matrix(dtm)
      g <- graph_from_incidence_matrix(m, weighted=TRUE)
      cl <- clustfun(g, E(g)$weight)
      grps <- communities(cl)  # get clusters
    }
    if (K > 1) {
      dtm <- TermDocumentMatrix(termlist, control = list(wordLengths = c(2, Inf), stemming=FALSE))
      m <- as.matrix(dtm)
      tf_mat <- TermDocFreq(t(m))
      # TF-IDF and cosine similarity:
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
    verts <- grepl("\\d+", V(g)$name)
    sdist <- as.dist(distances(g, v=verts, to=verts, weights=E(g)$weight))
    dups <- t(combn(1:ncol(x$data), m=2))
    dups <- cbind(dups, sdist)    
    dups <- dups[dups[,3]!=Inf, ]
    # old way:
    # dups <- lapply(lapply(grps, str_extract, "\\d+"), na.omit)
    # dups <- lapply(dups, as.numeric)
    # keep <- sapply(dups, length) > 1  # only keep clusters with >1 character
    # dups <- dups[keep]
    # grps <- grps[keep]
    # dups <- t(do.call(cbind, lapply(dups, combn, m=2)))
    # dups <- cbind(dups, NA)
    # colnames(dups) <- c('char1','char2','stringdist')
    # dmat <- distances(g, v=verts, to=verts)
    # dups[,3] <- sapply(1:nrow(dups), function(i) {
    #   dmat[as.character(dups[i,1]), as.character(dups[i,2])]
    # })
  }

  # Calculate text distances using fuzzy string distance
  if (opt == "fuzzy") {
    if (is.null(method)) {
      method <- "jw"
    }
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
    
    # individual text distance matrix
    sdist <- lapply(seq_along(splits), function(z) {
      nms <- splits[[z]]
      sd1 <- as.matrix(stringdistmatrix(part1[splits[[z]]], part1[splits[[z]]], method=method))
      sd2 <- as.matrix(stringdistmatrix(part2[splits[[z]]], part1[splits[[z]]], method=method))
      sd3 <- as.matrix(stringdistmatrix(part3[splits[[z]]], part1[splits[[z]]], method=method))
      dimnames(sd1) <- dimnames(sd2) <- dimnames(sd3) <- list(nms, nms)
      list(sd1, sd2, sd3)
    })
    
    # calculate weighted string distances
    sdist_final <- lapply(seq_along(sdist), function(z) {
      d1 <- weighting[1] * sdist[[z]][[1]]
      d2 <- weighting[2] * sdist[[z]][[2]]
      d3 <- weighting[3] * sdist[[z]][[3]]
      res <- matrix(mapply(sum,d1,d2,d3, MoreArgs=list(na.rm=T)),ncol=ncol(d1))/3
      diag(res) <- NA
      res[upper.tri(res)] <- NA
      dimnames(res) <- dimnames(d1)
      res
    })
    
    # setup duplicates matrix
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
    # string distances
    # stringdists <- unlist(sapply(seq_along(sdist_final), function(z) as.numeric(as.dist(sdist_final[[z]]))))
    # names(stringdists) <- 1:length(stringdists)
    dups <- subset(dups, stringdist < cutoff)
    # stringdists <- stringdists[stringdists < cutoff]
  }
  
  # Method to look at comments/reference to other datasets
  if (opt == "comments") {
    chars <- paste(oldcharnames, oldstatenames)
    # Search for author-year references (e.g., Eliason et al. 2015 character 13)
    m <- str_match_all(chars, "([A-Z][a-z]+\\s(?:et\\sal\\.\\s|(?:\\&|and)\\s[A-Z][a-z]+\\s)?)[\\(]?(\\d{4})[\\)]?.{1,5}[Ch]ar(\\.|acter)?\\s(\\d+)")
    names(m) <- seq_along(chars)
    # convert to matrix with only non-missing cases
    dups <- suppressWarnings(lapply(seq_along(m), function(z) {cbind(df_name = z, m[[z]])}))
    dups <- as.data.frame(do.call('rbind', dups))
    dups <- dups[, c(1:4, 6)]
    names(dups) <- c('charnum', 'text', 'author', 'year', 'charmatch')
    dups$id <- as.numeric(as.character(dups$charnum))
    dups$infile <- x$file[dups$id]
    dups$charnum <- x$charnum[dups$id]
    dups$matchfile <- paste0(tolower(str_extract(dups$author, "[A-Z][a-z]+")), "_", dups$year)
    head(dups)
    dups <- dups[, c("infile", "charnum", "author", "year", "charmatch", "matchfile")]
    # match 1
    id1 <- sapply(1:nrow(dups), function(i) {
      res <- x$file == dups$infile[i] & x$charnums==dups$charnum[i]
      ifelse(any(res), which(res), NA)
      })
    # match 2
    id2 <- sapply(1:nrow(dups), function(i) {
      res <- x$file==dups$matchfile[i] & x$charnums==dups$charmatch[i]
      ifelse(any(res), which(res), NA)
      })
    dups <- cbind(id1, id2)
    dups <- subset(dups, complete.cases(dups))
    dups <- cbind(dups, stringdists=NA)
  }

  # Method looking at trait correlations
  if (opt == "traitcor") {
    mat <- x$data
    mat <- gsub("\\?|\\-", NA, mat)  # convert '?' or '-' scorings to NA:
    pairids <- combn(1:ncol(mat), m=2)  # get all pairs of characters
    if (method == "polycor") {
      dmat <- pbmclapply(1:ncol(pairids), function(i) {
        id1 <- pairids[1, i]
        id2 <- pairids[2, i]
        char1 <- mat[, id1]
        char2 <- mat[, id2]
        M <- cbind(char1, char2)
        M <- na.omit(M)
        M <- apply(M, 2, as.numeric)
        suppressWarnings(tryCatch(polychor(M), error=function(e) NA))
      }, mc.cores = cores)
      cors <- unlist(dmat)
      weights <- rep(1, length(cors))
    }
    if (method == "hamming") {
      dmat <- pbmclapply(1:ncol(pairids), function(i) {
        id1 <- pairids[1, i]
        id2 <- pairids[2, i]
        char1 <- mat[, id1]
        char2 <- mat[, id2]
        M <- na.omit(cbind(char1,char2))
        wt <- nrow(M)
        M <- apply(M, 2, paste0, collapse="")
        sdist <- as.numeric(stringdistmatrix(M, method="hamming"))
        sdist <- sdist/str_length(M[1])  # standardize by number of characters
        sdist <- abs(sdist-0.5)
        sdist <- sdist/0.5
        c(sdist=sdist, weight=wt)
        }, mc.cores = cores)
      cors <- sapply(dmat, "[[", "sdist")
      weights <- sapply(dmat, "[[", "weight")  
    }
    stringdists <- 1-abs(cors)
    dups <- cbind(char1 = pairids[1,], char2 = pairids[2, ], cor = cors, weight = weights)
  }

  res <- x  # resultant NEXUS file for outputting later

  # create final duplicates matrix
  if (nrow(dups) == 0) {
    warning("No duplicates found, try a different method")
    res$dups <- NULL
  } else {
    res$dups <- data.frame("char1" = dups[, 1], "charnum1" = x$charnum[dups[, 1]], "char2" = dups[, 2], "charnum2" = x$charnum[dups[, 2]], stringdist = dups[, 3])
    # res$dups <- res$dups[order(res$dups$stringdist), ]
  }

  # get clusters
  if (opt == "terms") {
    res$clusters <- grps
    res$g <- g
    res$cl <- cl
  }
  
  if (opt == "fuzzy" & nrow(dups) != 0) {
    g <- graph_from_edgelist(t(apply(dups[,1:2,drop=F], 1, as.character)), directed=FALSE)
    # if (cutoff == Inf) {
      grps <- communities(components(g))
    # } else {
      # E(g)$weight <- 1/dups[,'stringdist']
      # cl <- clustfun(g, E(g)$weight)
      # grps <- communities(cl)
    # }    
    res$clusters <- grps
  }

  if (opt == "comments" & nrow(dups) != 0) {
    g <- graph_from_edgelist(t(apply(dups[,1:2,drop=F], 1, as.character)), directed=FALSE)
    cl <- clustfun(g)
    grps <- communities(cl)  # get clusters
    res$clusters <- grps
  }
  
  res$cutoff <- cutoff
  
  res

}


## sub-function to get groups of duplicated characters based on input map (char1 in column 1, char2 in column 2)
findgroups <- function(x) {
  g <- graph_from_data_frame(x, directed=FALSE)
  grps <- max_cliques(g)
  grps <- sapply(1:length(grps), function(x) {as.numeric(names(grps[[x]]))})
  names(grps) <- paste('dup', seq_along(grps),sep='')
  grps
}

# Small function to update a duplicated nex object
# changing cutoff will affect the characters identified as duplicates for later
# visualization (e.g., with `printout` or `duptree` functions)
update.nex <- function(x, cutoff=Inf) {
  dups <- x$dups
  newdups <- dups[dups$stringdist < cutoff, ]
  if (nrow(newdups)==0) {
    stop("No duplicates below the set cutoff")
  }
  g <- graph_from_edgelist(t(apply(newdups[,c('char1','char2'),drop=F], 1, as.character)), directed=FALSE)
  grps <- communities(components(g))
  x$clusters <- grps
  x$cutoff <- cutoff
  return(x)
}
