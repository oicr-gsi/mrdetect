library(data.table)
library(ggplot2)

#setwd('/Volumes/cgi/scratch/fbeaudry/plasmaWG/TGL49_0143/')

plasma_coverage_file <- 'detectionsPerSite.txt'

plasma_coverage <- read.table(plasma_coverage_file)

ggplot(plasma_coverage,aes(x=V1)) + 
  geom_histogram( bins = max(plasma_coverage$V1)) +
  # geom_density() + 
  theme_classic() + labs(x="Reads per site",y="Number of Sites") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_blank())