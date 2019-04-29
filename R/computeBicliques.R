# graph = igraph object
# k = number of shared terms (increasing this takes exponentially more time)
# l = minimum number of characters in a module
computeBicliques <- function(graph, k, l) {
  vMode1 <- c()
  if (!is.null(V(graph)$type)) {
    vMode1 <- which(!V(graph)$type)
    vMode1 <- intersect(vMode1, which(igraph::degree(graph) >= l))
  }
  nb <- get.adjlist(graph)
  bicliques <- list()
  if (length(vMode1) >= k) {
    comb <- combn(vMode1, k)
    i <- 1
    sapply(1:ncol(comb), function(c) {
      commonNeighbours <- c()
      isFirst <- TRUE
      sapply(comb[,c], function(n) {
        if (isFirst) {
          isFirst <<- FALSE
          commonNeighbours <<- nb[[n]]
        } else {
          commonNeighbours <<- intersect(commonNeighbours, nb[[n]])
        }
      })
      if (length(commonNeighbours) >= l) {
        bicliques[[i]] <<- list(m1=comb[,c], m2=commonNeighbours)
      }
      i <<- i + 1
    })
  }
  bicliques
}


# computeBicliques(g, k=2, l=1)

# M <- (rbind(c(1,0,1,1), c(1,1,1,1), c(1,0,1,0)))
# rownames(M) <- c('a', 'b', 'c')
# colnames(M) <- c('1', '2', '3', '4')
# V(g)$name
# g <- graph_from_incidence_matrix(M)
# plot(g, layout=layout_as_bipartite)
# system.time(bc <- computeBicliques(g, k=2, l=5)) # takes ~38s for 1227 chars x 1658 terms
# bc <- bc[!sapply(bc, is.null)]
# bc # m1 = terms in a module, m2 = characters in module

# table(V(g)$type)

# system.time(clust <- computeBicliques(g, k=2, l=2))

# clust <- clust[!sapply(clust, is.null)]

# length(clust)

# V(g)$name[clust[[6]]$m1]
# V(g)$name[clust[[6]]$m2]

# twig$charlab[c(1,854)]

# not weighted