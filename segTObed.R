args <- commandArgs(trailingOnly = TRUE)
seg <- read.table(args[1],header=TRUE)
seg$type <- 'NEU'
seg$type[seg$seg.mean < -0.3] <- 'DEL'
seg$type[seg$seg.mean > 0.3] <- 'DUP'
seg <- seg[,c('chrom','loc.start','loc.end','type','seg.mean')]
names(seg) <- c('#chr', 'start', 'end', 'type', 'log2')
write.table(
		  seg,
		  file = paste0(args[2],".seg.bed"),
		  append = F, quote = FALSE, sep = "\t", 
		  eol = "\n", na = "NA",dec = ".", row.names = FALSE, 
		  col.names = TRUE
)