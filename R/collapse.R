#' Write a nexus data object to a file
#'
#' Function to collapse taxa and merge data in taxon 1 with other taxa
#'
#' @param x (required) a `nex` object
#' @param map a list specifying equivalent taxa to map from and to (e.g., `c('Tinamus' = 'Tinamus_major')` will map all characters from Tinamus to Tinamus_major)
#' @param method whether to convert characters scored in both taxa to polymorphisms (`merge`) or retain original characters (default)
#' 
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' 
#' @author Chad Eliason \email{chad_eliason@@utexas.edu}
#' 
#' @export
#' 
collapse <- function(x, map, method = c('retain', 'merge')) {

# TODO maybe use mutate? or merge?

  res <- x

  method <- match.arg(method)

# add new characters to the matrix first

  added <- sort(unique(as.character(unlist(map))))

  kept <- intersect(added, x$taxlabels)

  # added <- setdiff(sort(unique(as.character(unlist(map)))), x$taxlabels)

  df <- data.frame(from = rep(names(map), lapply(map, length)),
                   to = unlist(map),
                   index = rep(seq_along(map), lapply(map, length)),
                   stringsAsFactors=FALSE)

  # if (any(df$from %in% df$to)) {
  #   stop("Taxa in 'from' are also in 'to' taxa. This is redundant -- please fix.")
  # }

  dropped <- names(map)[!names(map) %in% added]

  # this is tricky
  # if you are saying to map two species to a new terminal, then you need to run the above iteratively to
  # make sure you aren't overwriting the character data..

  # find new taxa that aren't already in the dataset
  if (any(!added %in% x$taxlabels)) {
    newid <- which(!added %in% x$taxlabels)
    res$taxlabels <- c(x$taxlabels, added[newid])
    newdata <- matrix(NA, nrow=length(newid), ncol=ncol(x$data))
    res$data <- rbind(x$data, newdata)
  }

  #
  if (!any(dropped %in% x$taxlabels)) {
    cat('No taxa to drop\n')
    return(x)
  }

  dropid <- which(x$tax %in% dropped)

  # loop over all original taxa - i.e. names of 'map' list

  for (i in 1:nrow(df)) {

      # extract character scorings
      id1 <- which(res$tax == df$from[i])
      id2 <- which(res$tax == df$to[i])
      scores1 <- res$data[id1, ]
      scores2 <- res$data[id2, ]

      # find overlaps in scored characters
      # scores in original, NAs in target:
      id.scored1 <- which(!is.na(scores1) & is.na(scores2))
      # NAs in original, scores in target:
      id.scored2 <- which(is.na(scores1) & !is.na(scores2))
      # different scores in target, original:
      id.diff <- which(scores1 != scores2)
      # same scores in target, orginal:
      id.overlap <- which(scores1 == scores2)

      # replace NAs in target/new with original/old taxon
      res$data[id2, id.scored1] <- scores1[id.scored1]
      
      # merge method
      if (method=='merge') {  
      
        # if data values are not the same, convert to polymorphisms
        if (any(id.diff)) {
            warning('Different character scorings [characters ', paste0(id.diff, sep=' '), '] between taxa ', df$from[i], ' and ', df$to[i], '; converting to polymorphisms')
            # splits into digits
            char1 <- stringr::str_extract_all(res$data[id1, id.diff], '[\\d\\-]', simplify=T)
            char2 <- stringr::str_extract_all(res$data[id2, id.diff], '[\\d\\-]', simplify=T)
            chars <- cbind(char1, char2)
            # convert to polymorphic characters
            newchars <- sapply(seq_along(id.diff), function(z) paste0(sort(unique(chars[z, ])), collapse=""))
            # add parentheses
            newchars <- ifelse(stringr::str_length(newchars) > 1, paste0('(',newchars,')'), newchars)
            # check if any gaps coded with normal scorings
            id.gap <- grep('-', newchars)
            if (any(id.gap)) {
              warning('Error in character', id.diff[id.gap], ': polymorphism containing a gap `-` in character 
                      and scored character; replacing with gap scoring')
              newchars[id.gap] <- gsub("[0-9()]", "", newchars[id.gap])
            }
            res$data[id2, id.diff] <- newchars  # replace data
        }
        # merge data if no overlap in scorings or if character scorings are the same
        if (any(id.overlap)) {
          warning('Redundant character scorings between taxa ', df$from[i], ' and ', df$to[i],
           '[characters ', paste0(id.overlap, sep=" "), ']; combining taxa')
        }
      }
      if (method=='retain') {
        # if data values are not the same, retain original taxa
        if (any(id.diff)) {
          warning('Different character scorings [characters ', paste0(id.diff, sep=' '), '] between taxa ', df$from[i], ' and ', df$to[i], '; retaining both taxa')
              id <- id[!id == id1] # don't drop
              dropped <- dropped[-i]
          } else if (any(id.overlap)) {
              warning('Redundant character scorings between taxa ', df$from[i], ' and ', df$to[i], '[characters ', paste0(id.overlap, sep=" "), ']; combining taxa')   
          }
      }
  }

  # drop taxa
  
  res$data <- res$data[-dropid, ]
  
  res$taxlabels <- res$taxlabels[-dropid]
  
  cat('Dropped taxa:', setdiff(dropped, kept), 'Added taxa:', setdiff(added, kept), sep='\n')
  
  # return results
  
  res

}
