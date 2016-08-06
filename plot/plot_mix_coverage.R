#! /usr/bin/env Rscript
#
# This script generates a plot of the read coverage across the reference
# genome for the complete mixture sample as well as the assembled contributor
# haplotypes and any reads that were not assigned to a contributor. The count
# of the majority base at each position is also plotted along with coverage
# so that variant position can be identified by eye.
#
# Sam Vohr (svohr@soe.ucsc.edu)


library(ggplot2)


args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("usage: plot_mix_coverage stat_file_prefix")
}

stats.prefix <- args[1]


obs.tab.fn <- paste(stats.prefix, "obs.tab", sep='.')
obs.tab <- read.table(obs.tab.fn)

obs <- data.frame(Contributor=paste(obs.tab[, 1], obs.tab[, 2], sep=':'),
                  position=obs.tab[, 3],
                  coverage=obs.tab[, 8])
obs$agreement <- apply(obs.tab[, 4:7], 1, max)

# Remove coverage and agreement values at positions with no coverage.
obs[obs$coverage == 0, c(3, 4)] <- NA
# Change the positions to 1-based for plotting.
obs$position <- obs$position + 1


cov.plot <- ggplot(obs) +
  geom_line(aes(x=position, y=agreement), size=0.3, color="lightblue") +
  geom_line(aes(x=position, y=coverage), size=0.6) +
  facet_wrap(~Contributor, ncol=1) +
  theme_minimal() +
  expand_limits(y = 0) +
  scale_x_continuous(breaks=c(1, 5000, 10000, 15000, 16569),
                     minor_breaks=seq(500, 16000, 500)) +
  labs(y="Coverage", x="RSRS Reference Position")

plot.fn <- paste(stats.prefix, "mix_coverage.png", sep='.')
n.facets <- length(levels(obs$Contributor))

ggsave(filename=plot.fn, plot=cov.plot, width=10, height=n.facets + 1, dpi=100)

quit(save='no', status=0)
