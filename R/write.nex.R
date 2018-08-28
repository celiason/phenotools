#' Write a nexus data object to a file
#'
#' Function to write nexus file to data for, e.g., analysis in MESQUITE, etc.
#'
#' @param x (required) nexus file (`nex`) object
#' @param file file patch for exported file
#' @param gap character representing incomparable data
#' @param missing character representing missing data
#' @return an object of class \code{nex} for use in further \code{nexustools} functions
#' @examples \dontrun{
#' x <- read.nex(file='example/toy1.nex')
#' write.nex(x2, file='output/testnexus.nex')
#' }
#' @author Chad Eliason \email{chad_eliason@@utexas.edu}
#'
write.nex <- function(x, file, missing=NULL, gap=NULL, mrbayes=FALSE, ngen=NULL,
  phy=NULL, run=FALSE, format=c("nexus", "tnt")) {

    format <- match.arg(format)

    if (is.null(missing)) {
      missing <- x$missing
    }

    if (is.null(gap)) {
      gap <- x$gap
    }
    dat <- as.matrix(x$data)
    dat <- ifelse(is.na(dat), missing, dat)
    zz <- file(file, 'w')
    ntax <- nrow(x$data)
    nchar <- ncol(x$data)
    
    # file concat function from read.nexus.data in the ape package
    fcat <- function(..., file = zz, sep = '') cat(..., file = file, sep = sep, append = TRUE)
  
  # Write character partitions
  if (!all(x$charpart=="''")) {
      test <- factor(x$charpartition)
      # Function for finding starts and lengths of string/number sequences
      seqle <- function(x,incr=1) { 
        if(!is.numeric(x)) x <- as.numeric(x) 
        n <- length(x)  
        y <- x[-1L] != x[-n] + incr 
        i <- c(which(y|is.na(y)),n) 
        list(lengths = diff(c(0L,i)),
             values = x[head(c(0L,i)+1L,-1L)]) 
      } 
      # Find sequence locations and lengths
      xx <- lapply(levels(test), function(x) {
          seqle(which(test==x))
        }
      )
      # Paste starts and ends of sequences
      res <- sapply(seq_along(xx), function(i) {
        starts <- xx[[i]]$values
        len <- xx[[i]]$lengths
        res <- ifelse(len>1, paste0(starts, "-", starts+(len-1)), starts)
        paste0(res, collapse=" ")
      })
      names(res) <- levels(test)
      charpart <- res
    }

  # Write character partition for file names
  # x <- nexs[[1]]
  if (!all(x$file=="''")) {
    filepart <- factor(x$file)
      # Find sequence locations and lengths
      xx <- lapply(levels(filepart), function(x) {
          seqle(which(filepart==x))
        })
      # Paste starts and ends of sequences
      res <- sapply(seq_along(xx), function(i) {
        starts <- xx[[i]]$values
        len <- xx[[i]]$lengths
        res <- ifelse(len>1, paste0(starts, "-", starts+(len-1)), starts)
        paste0(res, collapse=" ")
      })
      names(res) <- levels(filepart)
      filepart <- res
    # filepart
  }

  if (format=="nexus") {

    if (mrbayes) {
      fcat('#NEXUS\n[Data written by write.nex.R, ', date(), "]\n")
      fcat('begin data;\n')
      fcat('dimensions ntax=', ntax, ' nchar=', nchar, ';\n')
      fcat('format datatype=standard gap=- missing=',missing,' symbols="',x$symbols,'";\n')
      fcat('matrix\n')
      # write data matrix
      sapply(1:ntax, function(z) { fcat(x$taxlabels[z], '\t\t', dat[z,], '\n') })
      fcat(';\nend;\n')
      fcat('begin mrbayes;\n')
      if (!is.null(phy)) {
        ntips <- length(phy$tip.label)
        nnodes <- phy$Nnode
        ids <- Descendants(phy, (ntips+1):(ntips + nnodes), type = 'tips')
        ids <- lapply(1:length(ids), function(x) { phy$tip[ids[[x]]] })
        fcat(paste('constraint node', 1:phy$Nnode, ' = ', sapply(ids, paste0, collapse=" "), ';', sep=''), sep='\n')
        fcat('prset topologypr=constraints(', paste('node',1:phy$Nnode,sep='',collapse=','), ');\n', sep='')
      }
      fcat('prset ratepr=variable;\n')
      fcat('lset coding=all rates=gamma; [morphology, using all characters not just variable ones coding=variable for the latter]\n')
      # fcat('prset nodeagepr=calibrated;\n\n')
      fcat('[setup and run]', paste0('mcmc ngen=', as.integer(format(ngen, scientific=F)), ';'), 'sump burninfrac=0.4;', 'sumt burninfrac=0.4;', 'set nowarn=yes;', paste0('execute ', file), 'end;', sep='\n')
      close(zz)
      # taken from the ips package (CITE!)
      if (run) {
          if (.Platform$OS.type == "unix") {
              system(paste("mb", file))
          }
          else {
              system(paste("mrbayes ", file, ".bayes", sep = ""))
          }
          tr <- read.nexus(paste(file, ".con.tre", sep = ""))
          tr
      }
    } else {
        fcat('#NEXUS\n[Data written by write.nex.R, ', date(), "]\n")
        fcat('\tBEGIN TAXA;\n')
        fcat('\t\tDIMENSIONS NTAX=', ntax, ';\n')
        fcat('\t\tTAXLABELS\n')
        fcat(paste0("\t\t\t", "'", x$taxlabels, "'"), sep="\n")
        fcat('\t\t;\nEND;\n')
        fcat('\tBEGIN CHARACTERS;\n')
        fcat('\t\tDIMENSIONS ', ' NTAX=', ntax, ' NCHAR=', nchar, ';\n')
        fcat('\t\tFORMAT DATATYPE=STANDARD GAP=- MISSING=',missing,' SYMBOLS="',x$symbols,'";\n')
        if (!is.null(x$charlabels)) {
          fcat('\t\tCHARLABELS\n')
          fcat(paste('\t\t\t[', x$charnums, '] ', x$charlabels, sep=""), sep='\n')
          # fcat(paste("\t\t\t[", 1:ncol(dat), "] '", x$charlabels, "'", sep=""), sep="\n")
          fcat("\t\t;\n")
        }
        if (!is.null(x$statelabels)) {
          fcat('\t\tSTATELABELS\n')
          fcat(paste('\t\t\t', 1:nchar, ' ', x$statelabels, ',', sep=""), sep='\n')
          fcat('\t\t;\n')
        }
        fcat('\t\tMATRIX\n')
        # write data matrix
        sapply(1:ntax, function(z) { fcat("\t\t\t", "'", x$taxlabels[z], "'", "\t\t", dat[z,], "\n") })
        fcat('\t;\nEND;\n')

        # write character partitions
        fcat('BEGIN SETS;\n')
        fcat('\tCHARPARTITION set1=', paste0(names(charpart), ":", charpart, collapse=","), ';\n')
        if (!all(x$file=="''")) {
          fcat('\tCHARPARTITION file=', paste0(names(filepart), ":", filepart, collapse=","), ';\n')
        }
        fcat('END;\n')
        close(zz)
      }
  }

  if (format == "tnt") {

  dat <- gsub("\\(", "\\[", dat)
  dat <- gsub("\\)", "\\]", dat)
  fcat("xread\n'Data written by write.nex.R\n", date(), "'\n", sep="")
  fcat(nchar, " ", ntax, "\n")
  # write data matrix
  sapply(1:ntax, function(z) {fcat(x$taxlabels[z], " ", dat[z,], "\n")})
  fcat(";\n")
  close(zz)

  }

}
