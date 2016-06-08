#! /usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

stats.prefix <- args[1]

obs.tab.fn <- paste(stats.prefix, "obs.tab", sep='.')

obs.tab <- read.table(obs.tab.fn)

obs <- data.frame(Contributor=paste(obs.tab[,1], obs.tab[,2], sep=':'),
                  position=obs.tab[,3],
                  coverage=obs.tab[,8])
obs$agreement <- apply(obs.tab[,4:7], 1, max)

# Remove coverage and agreement values at positions with no coverage.
obs[obs$coverage == 0, c(3,4)] <- NA

# Change the positions to 1-based for plotting.
obs$position <- obs$position + 1


cov.plot <- ggplot(obs) +
  geom_line(aes(x=position, y=agreement), size=0.3, color="lightblue") +
  geom_line(aes(x=position, y=coverage), size=0.6) +
  facet_wrap(~Contributor, ncol=1) +
  theme_minimal() +
  labs(y="Coverage", x="RSRS Reference Position")

plot.fn <- paste(stats.prefix, "hap_coverage.png", sep='.')
n.facets <- length(levels(obs$Contributor))

ggsave(filename=plot.fn, plot=cov.plot, width=10, height=n.facets * 1, dpi=100)
