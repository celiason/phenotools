#' Filter duplicated characters in a nexus file
#' 
#' Given a list of specified duplicates, drop and/or merge characters
#' 
#' @param x (required) a nexus input file
#' @param dups an optional matrix or list specifying which character are duplicates (e.g., dups = list('1733' = c(1741,1745,1755)))
#' 
#' @examples \dontrun{
#' dat <- read.nex("/Users/chadeliason/Dropbox/phenome dataset/data/2015-09-02/original/final_reordered.nex")
#' drops <- grep("livezey_2006", dat$charlab)
#' tmp <- capture_comments(file = "~/Dropbox/phenome dataset/data/2015-09-02/modified/final_reorderedJAC+CME.txt")
#' ss <- tmp$markup
#' x <- dat[, ss]
#' dups <- tmp$dups
#' xx <- filter(x, dups=tmp$dups)
#' sum(is.na(x$data))/length(x$data)  # 86.5% missing data
#' sum(is.na(xx$data))/length(xx$data)  # 86.3% missing data
#' }
#' 
filter.nex <- function(x, dups=NULL) {
    if (class(dups)=="list") {      	
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
	  id <- as.numeric(dups[i, ])
	  # get character scorings
	  scores <- x$data[, match(id, x$charnum)]
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
	    if (any(apply(scores, 2, allNA))) {
	    	keepchar <- which(!apply(scores, 2, allNA))
	    	dropchar <- which(apply(scores, 2, allNA))
	    }
	    res$data[, match(id[keepchar], x$charnum)] <- newscores  # replace values
	    # res$data[, id[2]] <- newscores  
	    drops[i] <- id[dropchar]  # drop trait with fewer scorings
	    warning('Merging non-overlapping character scorings; dropping character ', id[1], '; check state labels to confirm')
	  }  
	}
	keep <- !res$charnum %in% drops
	################################################################################
	# Output, drop duplicated characters, labels, etc.:
	################################################################################
	res$data <- res$data[, keep]
	res$charlabels <- res$charlabels[keep]
	res$statelabels <- res$statelabels[keep]
	res$charset <- res$charset[keep]
	res$charnums <- res$charnums[keep]
	res$charpartition <- res$charpartition[keep]
	return(res)
}

# small function to count NAs
allNA <- function(x) {
	ifelse(all(is.na(x)), T, F)
}
