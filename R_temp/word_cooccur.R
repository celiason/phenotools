# word co-occurrence network for phenomes
# source: http://f.briatte.org/r/turning-keywords-into-a-co-occurrence-network

library(dplyr)
library(ggnetwork)
library(ggplot2)
library(readr)
library(stringr)
library(tnet)
library(network) # keep after tnet

x <- twig1$charlab
x <- gsub("0|1", "", x)
x <- cleantext(x)
x <- gsub("\\d*", "", x)
x <- gsub("[[:punct:]]", "", x)
tocut <- c("to", "the", "and", "an", "as", "a", "or", "of", "absent", "present", "form", "process", "state", "view", "margin", "shape", "placed", "recess", "(UN)?ORDERED", "(un)?ordered", "lateral", "distal", "ventral", "posterior", "anterior", "medial", "dors", "external")
tocut <- paste0("\\b", tocut, "\\b", collapse="|")
x <- gsub(tocut, "", x)
x <- gsub("½|¼", "", x)
x <- str_split(x, "\\s")
x <- lapply(x, unique)
x <- lapply(seq_along(x), function(i) { x[[i]][x[[i]]!=""] })



xx <- x[c(1, 500, 100)]


# weighted edge list
system.time(e <- xx %>%
  lapply(function(x) {
    expand.grid(x, x, w = 1 / length(x), stringsAsFactors = FALSE)
  }) %>%
  bind_rows)

# rate-limiting step:
# sorting..
system.time(e <- apply(e[, -3], 1, str_sort) %>%
  t %>%
  data.frame(stringsAsFactors = FALSE) %>%
  mutate(w = e$w))

# newman fowler weights
system.time(e <- group_by(e, X1, X2) %>%
  summarise(w = sum(w)) %>%
  filter(X1 != X2))

head(e)
tail(e)

# network object
n <- network(e[, -3], directed = FALSE)  # undirected network

stopifnot(nrow(e) == network.edgecount(n))
set.edge.attribute(n, "weight", e$w)

# weighted degree at alpha = 1
t <- as.edgelist(n, attrname = "weight") %>%
  symmetrise_w %>%
  as.tnet %>%
  degree_w

stopifnot(nrow(t) == network.size(n))
set.vertex.attribute(n, "degree_w", t[, "output" ])

# show only keywords at or above median weighted degree
l <- n %v% "degree_w"
l <- ifelse(l >= median(l), network.vertex.names(n), NA)

stopifnot(length(l) == network.size(n))
set.vertex.attribute(n, "label", l)



ggplot(n, aes(x, y, xend = xend, yend = yend)) +
  geom_edges(aes(color = weight)) +
  geom_nodes(color = "grey50") +
  geom_nodelabel(aes(size = degree_w, label = label),
                 color = "grey20", label.size = NA) +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_gradient2(low = "grey25", midpoint = 0.75, high = "black") +
  guides(size = FALSE, color = FALSE) + 
  theme_blank()


