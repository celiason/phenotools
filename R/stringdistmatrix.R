# modified function from stringdist package
stringdistmatrix <- function (a, b, method = c("osa", "lv", "dl", "hamming", "lcs", 
    "qgram", "cosine", "jaccard", "jw", "soundex"), useBytes = FALSE, 
    weight = c(d = 1, i = 1, s = 1, t = 1), maxDist = Inf, q = 1, 
    p = 0, bt = 0, useNames = c("none", "strings", "names"), 
    ncores = 1, cluster = NULL, nthread = getOption("sd_num_thread")) 
{
    if (maxDist < Inf) 
        warning("Argument 'maxDist' is deprecated for function 'stringdistmatrix'. This argument will be removed in the future.")
    if (ncores > 1) {
        warning("Argument 'ncores' is deprecated as stringdist now uses multithreading by default. This argument is currently ignored and will be removed in the future.")
        ncores <- 1
    }
    if (!is.null(cluster)) {
        message("Argument 'cluster' is deprecaterd as stringdust now uses multithreading by default. The argument is currently ignored and will be removed in the future")
    }
    if (is.list(a) || (!missing(b) && is.list(b))) {
        warning(listwarning("stringdistmatrix", "seqdistmatrix"))
    }
    if (identical(useNames, FALSE)) 
        useNames <- "none"
    if (identical(useNames, TRUE)) 
        useNames <- "strings"
    useNames <- match.arg(useNames)
    method <- match.arg(method)
    nthread <- as.integer(nthread)
    stopifnot(all(is.finite(weight)), all(weight > 0), all(weight <= 
        1), q >= 0, p <= 0.25, p >= 0, is.logical(useBytes), 
        ifelse(method %in% c("osa", "dl"), length(weight) >= 
            4, TRUE), ifelse(method %in% c("lv", "jw"), length(weight) >= 
            3, TRUE), ncores > 0, nthread > 0)
    if (method == "jw") 
        weight <- weight[c(2, 1, 3)]
    if (missing(b)) {
        if (useNames == "names") {
            a <- setNames(as.character(a), names(a))
        }
        else {
            a <- as.character(a)
        }
        return(stringdist:::lower_tri(a, method = method, useBytes = useBytes, 
            weight = weight, q = q, p = p, bt = bt, useNames = useNames, 
            nthread = nthread))
    }
    if (useNames == "names") {
        rowns <- names(a)
        colns <- names(b)
    }
    a <- as.character(a)
    b <- as.character(b)
    if (useNames == "strings") {
        rowns <- a
        colns <- b
    }
    if (!useBytes) {
        a <- enc2utf8(a)
        b <- enc2utf8(b)
    }
    if (length(a) == 0 || length(b) == 0) {
        return(matrix(numeric(0)))
    }
    x <- pbapply::pbsapply(b, stringdist:::do_dist, USE.NAMES = FALSE, a=a, method=method, weight=weight,
        maxDist=maxDist, q=q, p=p, bt=bt, useBytes=useBytes, nthread=nthread)
    if (useNames %in% c("strings", "names")) {
        structure(matrix(x, nrow = length(a), ncol = length(b), 
            dimnames = list(rowns, colns)))
    }
    else {
        matrix(x, nrow = length(a), ncol = length(b))
    }
}