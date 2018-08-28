get.biWeightedProjection <- function(A,vertex=FALSE,mode='shared-neighbours'){
    if (class(A) != "matrix"){
        stop("The function must apply to matrix object.\n")
        
    }
    g <- graph.incidence(A,mode = 'all')
    
    if (!bipartite.mapping(g)$res){
        stop("The function applies to bipartite graphs.\n")
    }

    if (is.bipartite(g)){
        V(g)$type <- bipartite.mapping(g)$type
        if (vertex) proj <- bipartite.projection(g)[[2]]
        else proj <- bipartite.projection(g)[[1]]
        V(proj)$name <- V(g)[type==vertex]$name
        ## number of shared vertices as edge weights 
        if (mode=='shared-neighbours'){
            incidenceMatrix <- get.incidence(g)	## TODO: might need some more testing if type and vType match
            if (!vertex){
                ## off-diagonal elements represent edge weights, diagonal elements vertex weights (currently unused)
                m <- incidenceMatrix %*% t(incidenceMatrix)
            }
            else {
                m <- t(incidenceMatrix) %*% incidenceMatrix
            }
            ## assign edge weights from matrix m
            if (length(E(proj))>0){
                for (i in 1:(length(E(proj)))){
                    E(proj)[i]$weight <- m[V(proj)[get.edge(proj,i)[[1]]]$name, V(proj)[get.edge(proj,i)[[2]]]$name]
                }
            }
            return(proj)
        }
        ## Newman2001 Formula
        if (mode=='newman'){
            incidenceMatrix <- get.incidence(g)
            if (!vertex){
                incidenceMatrix <- t(incidenceMatrix)
                if (length(c(which(rowSums(incidenceMatrix)==1)))!=0){
                    m2 <- incidenceMatrix[c(which(rowSums(incidenceMatrix)==1))*(-1),]
                }
                else{
                    m2 <- incidenceMatrix
                }
                m <- t(m2) %*% (m2*(1/(rowSums(m2)-1)))
            }
            else {	
                if (length(c(which(rowSums(incidenceMatrix)==1)))!=0){
                    m2 <- incidenceMatrix[c(which(rowSums(incidenceMatrix)==1))*(-1),]
                }
                else{
                    m2 <- incidenceMatrix
                }
                m <- t(m2) %*% (m2*(1/(rowSums(m2)-1)))
            }
            if (length(E(proj))>0){
                for (i in 1:(length(E(proj)))){
                    E(proj)[i]$weight <- m[V(proj)[get.edge(proj,i)[[1]]]$name, V(proj)[get.edge(proj,i)[[2]]]$name]
                }
            }
            return(proj)
        }
    }
}