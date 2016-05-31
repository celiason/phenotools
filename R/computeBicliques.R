computeBicliques <- function(graph, k, l) {

  vMode1 <- c()
  if (!is.null(V(graph)$type)) {

    vMode1 <- which(!V(graph)$type)
    vMode1 <- intersect(vMode1, which(degree(graph) >= l))
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
# g <- graph_from_incidence_matrix(M)

# plot(g, layout=layout_as_bipartite)

# bc <- computeBicliques(g, k=3, l=1)

# bc

# wc <- cluster_fast_greedy(g)

# wc

# plot(wc, g, layout=layout_as_bipartite)

# wc[[3]]

# plot(wc[1])

# membership(wc)

# plot(induced_subgraph(g, membership(wc)==3))
