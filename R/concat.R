#' Merge a list of nexus files
#'
#' A function that combines numerous nexus datasets
#'
#' @param x (required) a list of `nex` objects to merge
#' @param taxa a character vector of names of taxa to include in concatenated nexus file (optional)
#' 
#' @return an object of class `nex` for use in further `phenotools` functions
#' 
#' @examples \dontrun{
#' x <- read.nex(system.file("extdata", "clarke_2006.nex", package = "phenotools"))
#' y <- read.nex(system.file("extdata", "nesbitt_2015.nex", package = "phenotools"))
#' xy <- concat(list(x, y), taxa = c('Anas_platyrhynchos', 'Gallus_gallus', 'Crypturellus_undulatus'))
#' xy
#' }
#' 
#' @author Chad Eliason \email{celiason@@fieldmuseum.org}
#' 
#' @export
#' 
concat <- function(x, taxa=NULL) {
  mat <- lapply(x, '[[', 'data')
  alltaxlabels <- unique(unlist(lapply(x, '[[', 'taxlabels')))
  # get taxon labels for each input dataset
  taxbydataset <- lapply(x, '[[', 'taxlabels')
  # find overlapping species
  if (!is.null(taxa)) {
    ids <- lapply(1:length(taxbydataset), function(z) {match(taxa, taxbydataset[[z]])})
    taxlabels <- taxa
  } else {
    ids <- lapply(1:length(taxbydataset), function(z) {match(alltaxlabels, taxbydataset[[z]])})
    taxlabels <- alltaxlabels
  }
  # combine data
  dat <- do.call(cbind, lapply(1:length(mat), function(z) {mat[[z]][ids[[z]], ]}))
  rownames(dat) <- 1:length(taxlabels)
  names(dat) <- 1:ncol(dat)
  nchar <- ncol(dat)
  ntax <- nrow(dat)
  charlabels <- unlist(lapply(x, '[[', 'charlabels'))
  charnums <- unlist(lapply(x, '[[', 'charnums'))
  charsets <- unlist(lapply(x, '[[', 'charset'))
  file <- unlist(lapply(x, '[[', 'file'))
  charpartition <- unlist(lapply(x, '[[', 'charpartition'))
  statelabels <- unlist(lapply(x, '[[', 'statelabels'))
  missing <- unlist(lapply(x, '[[', 'missing'))
  gap <- unlist(lapply(x, '[[', 'gap'))
  symbols <- paste(unique(unlist(strsplit(gsub('[\\(\\)\\??\\-]', '', sort(unique(as.vector(dat)))), ""))),collapse="")
  # symbols <- paste(sort(unique(as.vector(mat))), collapse="")
  symbols <- gsub('[\\?\\-]','',symbols)
  if (length(unique(missing)) != 1) {
    stop('Datasets differ in characters used for missing data, fix before merging')
  } else { missing <- unique(missing) }
  if (length(unique(gap)) != 1) {
    stop('Datasets differ in characters used for gap data, fix before merging')
  } else { gap <- unique(gap) }
  res <- list(taxlabels = taxlabels, data = dat, charlabels = charlabels, statelabels = statelabels, charnums = charnums, charset = charsets, missing = missing, gap = gap, charpartition = charpartition, symbols = symbols, file = file)
  class(res) <- c('nex', 'list')
  res
}
