#' plot data matrix and phylogeny
#' plot a nexus file
#' @param bw whether to plot in black (scoring present) and white (no scoring)
#' @param colour any factor present in the nexus file (statelabels, file, charpartition, charset, etc.)
#'
#' x <- allnex.genome
#' phy <- tree

# TODO compute changes along different branches, set as branch thickness
# TODO create a shiny app for exploring phenomic data?

plot.nex <- function(x, phy = NULL, legend.pos = 'none', bw = FALSE, na.value = 'lightgray', fsize = 8, fill = c('statelabels', 'charpartition', 'charset', 'file')) {

	require(grid)
	require(ape)
	require(gridExtra)
	require(RColorBrewer)

  	source('/Users/chadeliason/R/ggplotphylo.R')

  	fill <- match.arg(fill)

	xx <- 1:ncol(x$data)
	yy <- x$taxlabels
	df <- expand.grid(x=xx,y=yy)
	df$z <- zz <- as.character(t(x$data))
	df$z <- ifelse(df$z %in% c(x$missing, x$gap), NA, df$z)

	if (fill == 'statelabels') {
		df$fill <- as.character(t(x$data))
		df$fill <- factor(ifelse(df$fill %in% c(x$missing, x$gap), NA, df$fill))
	}

	# need to think about this. doesn't make sense now to have whole block colored one way. maybe alpha for character states?

	if (fill == 'charpartition') {
		df$fill <- ifelse(is.na(df$z), 0, 1)
		df$charpartition <- x$charpartition[df$x]
		df$fill <- ifelse(is.na(df$z), NA, df$charpartition)
	}

	if (fill == 'file') {
		df$fill <- ifelse(is.na(df$z), 0, 1)
		df$file <- x$file[df$x]
		df$fill <- ifelse(is.na(df$z), NA, df$file)
	}

	# if (fill == 'charset') {
	# 	df$fill <- as.factor(x$charset)
	# }

	# df$fill <- as.factor(t(x$data))

	if (length(levels(df$fill)) < 12) {
		# pal <- scale_fill_brewer(palette = 'Set3', na.value = na.value)
		colourCount = length(unique(df$fill))
		pal <- scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(colourCount))
	} else {
		pal <- scale_fill_discrete(na.value = na.value)
	}

	if (bw) {
		df$fill <- ifelse(df$fill=='-', 'gray', df$fill)
		df$fill <- ifelse(is.na(df$fill), 'white', df$fill)
		df$fill <- ifelse(!df$fill %in% c('white','gray'), 'black', df$fill)
		# df$fill <- ifelse(is.na(df$fill), 'white', 'black')
		pal <- scale_fill_identity(na.value = na.value)
	}

	if (!is.null(phy)) {

		# only keep data in phylogeny, and match tips
		if (any(phy$tip %in% x$tax)) {
			df <- df[df$y %in% phy$tip, ]
			# x <- x[x$tax %in% phy$tip, ]
			phy <- drop.tip(phy, which(!phy$tip %in% df$y))
			tip.order <- phy$tip.label[phy$edge[phy$edge[, 2] <= Ntip(phy), 2]]
			matches <- match(df$y, tip.order)
		} else {
			stop('Tip labels not found in data')
		}

		# create ggplot phylogeny object
		p1 <- ggplot(data = phy) +
				geom_segment(aes(x = x, y = y, xend = xend, yend = yend), colour = "black", alpha = 1, lineend='square') +
				coord_cartesian(y = c(.75, Ntip(phy) + .75)) +
				theme(
					legend.position = "none",
					axis.text = element_text(colour = NA),
					axis.text.x = element_blank(),
					axis.text.y = element_blank(),
					axis.ticks = element_line(colour = NA),
					axis.title.x = element_text(colour = NA),
					axis.title.y = element_blank(),
					axis.line = element_line(colour = NA),
					plot.margin = unit(c(0,0,0,0), "cm"),
					plot.background = element_rect(fill = NA),
					panel.background = element_rect(fill= NA),
					panel.border = element_rect(fill = NA, colour = NA),
					panel.grid = element_blank()
				)
		
		df$matches <- match(df$y, tip.order)
		
		p2 <- ggplot(df, aes(x, reorder(y, matches), fill = fill)) + geom_tile() +
			  	xlab("Trait") + 
				coord_cartesian(y = c(.75, Ntip(phy) + .75)) +
				pal +
			  	scale_x_continuous(breaks = seq(0, ncol(x$data), by=10)) +
			  		theme(
						legend.position = "none",
						axis.text.x = element_blank(),
						axis.text.y = element_text(hjust = 0, size = fsize),
						axis.title.x = element_text(colour = NA),
						axis.title.y = element_blank(),
						axis.line.y = element_line(colour=NA), axis.text.y=element_blank(),
						axis.ticks.y = element_line(colour = NA), axis.title.y = element_blank(),
						plot.margin = unit(c(0,0,0,0), "cm"),
						panel.background = element_blank(), #element_rect(fill = NA, colour = NA),
						panel.border = element_blank(), #element_rect(fill=NA, colour = NA),
						panel.grid = element_blank()
						)

		grid.arrange(p1, p2, nrow = 1, widths = c(.2, 1))
		
		return(list(phy = phy, data = x))

	}

	if (is.null(phy)) {
		
		ggplot(df, aes(x, y, fill = abbreviate(fill, minlength=10))) + geom_tile() +
		  	xlab("Trait") + 
		  	ylab("Species") +
		  	ylim(rev(levels(df$y))) +
		  	theme(
		  		legend.position = legend.pos,
		  		axis.text.x = element_blank(),
				axis.text.y = element_text(hjust = 0, size = fsize),
				axis.line.y = element_line(colour=NA),
				axis.ticks.y = element_line(colour = NA), axis.title.y = element_blank(),
				plot.margin = unit(c(0,0,0,0), "cm"),
				panel.background = element_blank(),
				panel.border = element_blank(),
				panel.grid = element_blank(),
				legend.text = element_text(size = fsize * .75),
				legend.title = element_blank(),
				legend.key.size = unit(8, 'points')
				) +
		  	scale_x_continuous(breaks = seq(0, ncol(x$data), by=10)) +
		  	pal +
		  	guides(fill=guide_legend(nrow=2, byrow=TRUE))
	}

}
