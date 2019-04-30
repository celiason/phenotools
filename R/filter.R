#' Filter duplicated characters in a nexus file
#' 
#' Given a list of specified duplicates, drop and/or merge characters
#' 
#' @param x (required) a nexus input file
#' @param dups an optional matrix or list specifying which character are duplicates (e.g., dups = list('1733' = c(1741,1745,1755)))
#' 
#' @examples \dontrun{
#' x <- read.nex(system.file("extdata", "clarke_2006.nex", package = "phenotools"))
#' y <- read.nex(system.file("extdata", "nesbitt_2015.nex", package = "phenotools"))
#' xy <- concat(list(x,y))
#' dups <- duplicated(xy, opt="terms")
#' # a pair of duplicate characters
#' xy[,c(144,345)]$charlab
#' # drop from dataset
#' xy2 <- filter(xy, cbind(144,345))
#' xy
#' xy2
#' }
#' 
#' @export
#' 
filter <- function(x, dups = NULL) {
    if (class(dups) == "list") {      	
      	# dups <- matrix(unlist(dups), ncol=length(dups))
      	dups <- cbind(as.numeric(rep(names(dups), sapply(dups, length))), as.numeric(unlist(dups)))
    }
	if (is.null(dups)) {
		dups <- x$dups
	}
	dups <- as.matrix(dups)
	res <- x
	# keep <- rep(TRUE, ncol(res$data))
	drops <- rep(NA, nrow(dups))
	# drops <- numeric(length = nrow(dups))
	for (i in seq_along(drops)) {
	  # i=1
	  id <- as.numeric(dups[i, ])
	  # get character scorings
	  scores <- x$data[, match(id, x$charnum)]
	  scores1 <- scores[,1]
	  scores2 <- scores[,2]
	  # create reverse scores (in the case where someone scores 0 & 1 and someone
	  # else 1 & 0 for the same structure)
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
	    # if both traits have scored the same taxa, but scores differ either
	    # keep (if just reversed scores) or leave alone
	    if (all(scores1[cc] == scores2rev[cc], na.rm=TRUE)) {
	      newscores <- scores1
	      drops[i] <- id[1]
	      message('Scores differ systematically (reversed values); dropping trait ',
	      	id[1], '; check state labels to confirm')
	    }
	    if (!all(scores1[cc] == scores2[cc]) & !all(scores1[cc] == scores2rev[cc])) {
	      warning('Traits differ in their scorings. Keeping both.')
	    }
	  }
	  # Case 2: Different number of taxa scores, some overlapping scores
	  if (diff(scorings) != 0 & n.overlap > 0) {
	    # if trait 1 scores more taxa, overlapping scores same - keep character
	    # with more scorings
	    if (all(scores1[cc] == scores2[cc], na.rm=TRUE)) {
	      drops[i] <- id[which.min(scorings)]  # drop trait with fewer scorings
	      message('Score overlap identical; keeping character with more scorings
	      	and dropping character ', id[which.min(scorings)])
	    }
	    # if trait 1 scores more taxa, overlapping scores different - if reversed keep:
	    if (all(scores1[cc] == scores2rev[cc], na.rm=TRUE)) {
	      drops[i] <- id[which.min(scorings)]
	      warning('Assuming characters are equivalent (reversed scorings) and
	      	dropping character with fewer scorings (', id[2], '); check state
	      	labels to confirm')
	    }
	    # otherwise do nothing:
	    if (!all(scores[cc,1] == scores[cc,2]) & !all(scores[cc,1] == rev(scores[cc,2]))) {
	      warning('Traits differ in their scorings; keeping both')
	    }
	  }
	  # Case 3: trait 1 scores non-overlapping with trait 2 scores; merge characters
	  if (n.overlap == 0) {
	    newscores <- pmin(scores[,1], scores[,2], na.rm=TRUE)
	    if (any(apply(scores, 2, allNA))) {
	    	keepchar <- which(!apply(scores, 2, allNA))
	    	dropchar <- which(apply(scores, 2, allNA))
	    }
	    res$data[, match(id[keepchar], x$charnum)] <- newscores  # replace values
	    # res$data[, id[2]] <- newscores  
	    drops[i] <- id[dropchar]  # drop trait with fewer scorings
	    warning('Merging non-overlapping character scorings; dropping character ',
	    	id[1], '; check state labels to confirm')
	  }  
	}

	keep <- !1:ncol(res$data) %in% drops

	# Output, drop duplicated characters, labels, etc.:
	res$data <- res$data[, keep]
	res$charlabels <- res$charlabels[keep]
	res$statelabels <- res$statelabels[keep]
	res$charset <- res$charset[keep]
	res$charnums <- res$charnums[keep]
	res$charpartition <- res$charpartition[keep]
	return(res)
}
