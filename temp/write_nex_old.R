# old version!!

# Mon Aug 10 16:04:53 2015

#' Write a nexus data object to a file
#'
#' Function to collapse taxa and merge data in taxon 1 with other taxa
#'
#' @param x (required) a `nex` object
#' @param file file patch for exported file
#' @param map a list specifying equivalent taxa to map from and to (e.g., `c('Tinamus' = 'Tinamus_major')` will map all characters from Tinamus to Tinamus_major)
#' @param method whether to convert characters scored in both taxa to polymorphisms (`merge`) or retain original characters (default)
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' @examples \dontrun{
#'
#' x <- read.nex(file='example/toy1.nex')
#' write.nex(x2, file='output/testnexus.nex')
#'
#' @author Chad Eliason \email{chad_eliason@@utexas.edu}
#'
x <- read.nex(files[4])
map <- list('Ancestor' = c('Tinamus_solitarius', 'Crypturellus_undulatus', 'Rhynchotus_rufescens'), 'Dromaius' = 'Rhynchotus_rufescens')
# map
# maybe use mutate? or merge?
collapse.nex <- function(x, map, method = c('retain', 'merge')) {

  res <- x

  method <- match.arg(method)


# add new characters to the matrix first

# added <- setdiff(unlist(map), x$taxlabels)

added <- setdiff(sort(unique(as.character(unlist(map)))), x$taxlabels)

df <- data.frame(from = rep(names(map), lapply(map, length)), to = unlist(map),
  index = rep(seq_along(map), lapply(map, length)), stringsAsFactors=FALSE)

# library(igraph)

g <- graph_from_data_frame(df, directed=T)

if (any(df$from %in% df$to)) {
  stop("Taxa in 'from' are also in 'to' taxa. This is redundant -- please fix.")
}

dropped <- names(map)

# this is tricky
# if you are saying to map two species to a new terminal, then you need to run the above iteratively to
# make sure you aren't overwriting the character data..

  # find new taxa that aren't already in the dataset
  if (any(!added %in% res$taxlabels)) {
    newid <- which(!added %in% res$taxlabels)
    res$taxlabels <- c(res$taxlabels, added[newid])
    newdata <- matrix(NA, nrow=length(newid), ncol=ncol(res$data))
    res$data <- rbind(res$data, newdata)
    # fill in data (do this later I think...)
    # for (i in newid) {
    #   id1 <- match(df$to[i], res$taxlabels)
    #   # id2 <- match(df$from[i], res$taxlabels)
    #   # res$data[id1, ] <- res$data[id2, ]
    # }
  }

  dropid <- which(x$tax %in% names(map))

  # loop over all original taxa - i.e. names of 'map' list
  
  for (i in 1:nrow(df)) {
    
  }

  for (i in seq_along(names(map))) {
    # and taxa merging with...
    for (j in seq_along(map[[i]])) {
      id1 <- which(x$tax == names(map)[[i]])				
      id2 <- which(x$tax == map[[i]][j])
      # extract character scorings
      scores1 <- x$data[id1, ]
      scores2 <- x$data[id2, ]
      # find overlaps in scored characters
      id.scored1 <- which(!is.na(scores1) & is.na(scores2))  # score in original, NA in target
      id.scored2 <- which(is.na(scores1) & !is.na(scores2))  # NA in original, score in target
      id.diff <- which(scores1 != scores2)  # different scores in target, original
      id.overlap <- which(scores1 == scores2)  # same scores in target, orginal
      # replace NA in target with original
      res$data[id2, id.scored1] <- scores1[id.scored1]
      if (method=='merge') {
        # if data values are not the same, convert to polymorphisms
        if (any(id.diff)) {
            warning('Different character scorings [characters ', paste0(id.diff, sep=' '), '] between taxa ', names(map)[i], ' and ', map[[i]][j], '; converting to polymorphisms')
            # splits into digits
            char1 <- str_extract_all(x$data[id1, id.diff], '[\\d\\-]', simplify=T)
            char2 <- str_extract_all(x$data[id2, id.diff], '[\\d\\-]', simplify=T)
            chars <- cbind(char1, char2)
            # convert to polymorphic characters
            newchars <- sapply(seq_along(id.diff), function(z) paste0(sort(unique(chars[z, ])), collapse=""))
            # add parens
            newchars <- ifelse(str_length(newchars) > 1, paste0('(',newchars,')'), newchars)
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
          warning('Redundant character scorings between taxa ', names(map)[i], ' and ', map[[i]][j],
           '[characters ', paste0(id.overlap, sep=" "), ']; combining taxa')
        }
      }
      if (method=='retain') {
        # if data values are not the same, retain original taxa
        if (any(id.diff)) {
          warning('Different character scorings [characters ', paste0(id.diff, sep=' '), '] between taxa ', names(map)[i], ' and ', map[[i]][j], '; retaining both taxa')
              id <- id[!id == id1] # don't drop
              dropped <- dropped[-i]
          } else if (any(id.overlap)) {
              warning('Redundant character scorings between taxa ', names(map)[i], ' and ', map[[i]][j], '[characters ', paste0(id.overlap, sep=" "), ']; combining taxa')   
          }
      }
    }
  }

  # drop taxa
  res$data <- res$data[-dropid, ]
  res$taxlabels <- res$taxlabels[-dropid]
  cat('Dropped taxa:', dropped, 'Added taxa:', added, sep='\n')
  # return results
  res  
}
