#! /usr/bin/env Rscript
#
# This script plots the coverage of each contributor + unassigned reads in a
# single plot of the mtDNA reference.
#
# Sam Vohr (svohr@soe.ucsc.edu)


library(ggplot2)


set_zero_to_empty <- function(x) {
  # sets the 0 values of vector X to NA unless it is the first or last in a
  # run of 0's.
  v <- x
  for (i in 1:length(x)) {
    if (v[i] == 0) {
      if (i > 1 && !is.na(v[i - 1]) && v[i-1] > 0) {
        v[i] <- 0
      }
      else if (i < length(v) && !is.na(v[i + 1]) && v[i + 1] > 0) {
        v[i] <- 0
      }
      else {
        v[i] <- NA
      }
    }
  }
  return(v)
}


# Check the command-line argument
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("usage: plot_mix_coverage stat_file_prefix")
}

stats.prefix <- args[1]


obs.tab.fn <- paste(stats.prefix, "obs.tab", sep='.')
orig.tab <- read.table(obs.tab.fn)
obs.tab <- orig.tab[orig.tab$V1 != 'all' & orig.tab$V1 != 'unassigned', ]

obs <- data.frame(Contributor=paste(obs.tab[, 1], obs.tab[, 2], sep=':'),
                  position=obs.tab[, 3],
                  coverage=obs.tab[, 8])
obs$agreement <- apply(obs.tab[, 4:7], 1, max)

unasn.tab <- orig.tab[orig.tab$V1 == 'unassigned', ]
unasn <- data.frame(Contributor=unasn.tab[, 1],
                    position=unasn.tab[, 3],
                    coverage=unasn.tab[, 8])

# Remove coverage and agreement values at positions with no coverage.
obs[obs$coverage == 0, 4] <- NA
obs$coverage <- set_zero_to_empty(obs$coverage)

# Change the positions to 1-based for plotting.
obs$position <- obs$position + 1
unasn$position <- unasn$position + 1

cov.plot <- ggplot()

if (nrow(unasn) > 0) {
  cov.plot <- cov.plot +
    geom_area(data=unasn, aes(x=position, y=coverage), alpha=0.1)
}

cov.plot <- cov.plot +
  geom_line(data=obs,
            aes(x=position, y=coverage, color=Contributor), size=0.6) +
  theme_minimal() +
  scale_x_continuous(breaks=c(1, 5000, 10000, 15000, 16569),
                     minor_breaks=seq(500, 16000, 500)) +
  labs(y="Coverage", x="RSRS Reference Position")

plot.fn <- paste(stats.prefix, "hap_coverage.png", sep='.')

ggsave(filename=plot.fn, plot=cov.plot, width=15, height=4, dpi=100)

quit(save='no', status=0)

