# user-generated ontology

test <- subset(allnex2, allnex2$file == "mayr_2011")

test$charlab

x <- str_replace_all(test$charlab, "[\\s]*[\\[\\(].+[\\]\\)][\\s]*", "")  # remove comments in square brackets

x

# split up into words
# x <- sapply(x, str_split, "\\s")


onto <- readLines("baumel ontology.txt")
onto

?read.table
onto <- read.table("baumel ontology.txt", sep="\t", header=FALSE, fill=TRUE, na.strings="", strip.white=TRUE, stringsAsFactors=FALSE, col.names=1:8)
head(onto)

# tail(sapply(onto, na.locf, na.rm=FALSE))



onto <- gsub("\\t", "", onto)
onto <- schinke(onto)
onto <- sapply(onto, str_split, "\\s")

grep(onto[[2]][1], tolower(x), value=TRUE)


# 


str_detect(x[[50]], onto[[5]])

x[50]

onto[5]



matches <- lapply(terms$"search term", grep, x = tomatch, ignore.case = TRUE, perl = TRUE)



pairs <- expand.grid(tolower(x), tolower(onto))

system.time(res <- mapply(stringdist, pairs[,1], pairs[,2], method="lcs"))

res2 <- matrix(res, nrow = length(x), ncol = length(onto))

image(res2)

rownames(res2) <- paste0("char", seq_along(x))

plot(hclust(dist(res2)), cex=.5)

which.min(res2[2, ])

onto[276]
x[2]

# yuck - I think I need to figure out this latin stemming thing.
# this matched "Os parietale" to "Ossa palatina completely fused along midline"


res[3, 208]


g <- graph_from_incidence_matrix(res[1:10, 1:100], weighted = TRUE)
g
pdf(file = "~/Desktop/test.pdf")
plot(g, edge.width = 5/E(g)$weight, vertex.size = 0, layout = layout_as_bipartite)
dev.off()


sapply(x, stringdist, onto[3], method="lcs")

plot(unlist(lapply(x, stringdist, onto[3], method="lcs")))



class(x)

length(x)




ontology <- str_match_all(x, "(.*?)(\\,|\\:|\\.|$)")


# figure out if word is 


# nested terms
ontology <- lapply(ontology, "[", , 2)


# remove blanks
ontology <- lapply(seq_along(ontology), function(x) {
	ontology[[x]][!ontology[[x]]==""]
})


# remove spaces at start/end of words
ontology <- lapply(seq_along(ontology), function(x) {gsub("^ | $", "", ontology[[x]])})

ontology <- lapply(ontology, tolower)

ontology

res <- c(NA, NA)

for (i in seq_along(ontology)) {
	ss <- ontology[[i]]
	l <- length(ss)
	for (i in (l-1):1) {
		res <- rbind(res, c(ss[i+1], ss[i]))
	}
}

res <- res[-1, ]

res <- res[!grepl("[Ff]orma", res[,1]), ]

res <- res[!grepl("[Ss]tatus", res[,1]), ]

res <- res[complete.cases(res), ]

# id <- paste(res[,1], res[,2])

# res <- res[!duplicated(id), ]

res <- gsub("\\(.*?\\)", "", res)

res

g <- graph_from_data_frame(res)

g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

g

plot(g)


# should we use a previously generated ontology to search on?

id <- grep("[Ff]em", V(g)$name)

# ego(g, 4, mode = "in")

id

subcomponent(g, id[1], mode="out")


plot(induced_subgraph(g, id))


unique(V(g)$name)

# now, for the graph/ontology, the idea is to get keywords to search on
# then wherever they're found, use the ontology to link them to higher-up traits

# maybe use string distances?



# find location terms (distal, proximal, etc.)

chars <- allnex2$charlabel

distal <- grep("distal", chars)
proximal <- grep("proxim", chars)
anterior <- grep("cran|anteri", chars)
lateral <- grep("latero|lateral", chars)
posterior <- grep("caud|poster", chars)

table(allnex2$charpart[distal])
table(allnex2$charpart[proximal])
table(allnex2$charpart[anterior])
table(allnex2$charpart[posterior])
table(allnex2$charpart[lateral])
