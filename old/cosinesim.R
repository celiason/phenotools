# testing

library(textmineR)

dtm <- TermDocumentMatrix(termlist, control = list(wordLengths = c(2, Inf), weighting = function(x) weightTfIdf(x, normalize = TRUE), stemming=TRUE))  

dtm <- TermDocumentMatrix(termlist, control = list(wordLengths = c(2, Inf)))

dtm <- as.matrix(dtm)
head(dtm)

# dtm <- as.matrix(TermDocumentMatrix(tmp, control = list(weighting = function(x) weightTfIdf(x, normalize = T))))
# dtm <- as.matrix(TermDocumentMatrix(tmp, control = list(weighting = function(x) weightTfIdf(x, normalize = F))))

tf_mat <- TermDocFreq(t(dtm))

# TF-IDF and cosine similarity
tfidf <- t(t(dtm)[ , tf_mat$term ]) * tf_mat$idf

tfidf <- t(tfidf)

# The next step is to calculate cosine similarity and change it to a distance. We’re going to use some linear algebra to do this. The dot product of two positive-valued, unit-length vectors is the cosine similarity between the two vectors. For a deeper explanation of the math and logic, read this article.

csim <- tfidf / sqrt(rowSums(tfidf * tfidf))

csim <- csim %*% t(csim)

csim

head(csim)

g2 <- graph_from_adjacency_matrix(csim, weighted=TRUE)

plot(g2)

communities(cluster_walktrap(g2))


# R’s various clustering functions work with distances, not similarities. We convert cosine similarity to cosine distance by subtracting it from 1
# This works because cosine similarity is bound between 0 and 1
# While we are at it, we’ll convert the matrix to a dist object.

cdist <- as.dist(1 - csim)

range(csim)



# The last step is clustering. There are many clustering algorithms out there. My preference is agglomerative hierarchical clustering using Ward’s method as the merge rule. Compared to other methods, such as k-means, hierarchical clustering is computationally inexpensive.

# In the example below, I choose to cut the tree at 10
# 10
#  clusters. This is a somewhat arbitrary choice. I often prefer to use the silhouette coefficient. You can read about this method here. Performing this is an exercise I’ll leave to the reader.

hc <- hclust(cdist, "ward.D")

K <- 100

clust <- cutree(hc, K)

# print character pairs
strwrap(x$charlab[which(clust==1)], simplify=F)
strwrap(x$charlab[which(clust==2)], simplify=F)
strwrap(x$charlab[which(clust==3)], simplify=F)

clustering <- cutree(hc, K)

plot(hc, main = "Hierarchical clustering of 100 NIH grant abstracts",
     ylab = "", xlab = "", yaxt = "n")
rect.hclust(hc, K, border = "red")

# v <- sort(rowSums(m), decreasing=TRUE)
# d <- data.frame(word = names(v), freq = v)

g <- graph_from_incidence_matrix(dtm, weighted = TRUE)

g

cl <- cluster_fast_greedy(g, weights = E(g)$weight)
grps <- communities(cl)

plot(cl, g)

plot(g, edge.width = 2*E(g)$weight)


# m[c('carpometacarp','intermetacarp','pisiform','spac','supratrochlear'), c(389,390,391)]
# methclust?
# plot(jitter(as.numeric(m)))
# m[1:5, 1:5]




# meth clust

library(methClust)

dtm <- DocumentTermMatrix(termlist, control = list(weighting = weightBin, wordLengths = c(2, Inf), stemming=TRUE))

head(as.matrix(dtm))

reg_data <- as.matrix(dtm)

dim(reg_data)

res <- meth_topics(meth=reg_data, unmeth=1-reg_data, K=5, tol=10, use_squarem=FALSE)

str(res)
head(res$freq)
head(res$omega)
plot(res$omega[5, ], type='o', ylim=c(0 ,1))


strwrap(x$charlab[c(1,5)], simplify=F)
