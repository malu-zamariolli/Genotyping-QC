het_qc = snakemake@input[["het_qc"]]
het_out = snakemake@output[["het_out"]]

dat <- read.table(het_qc, header=T) 
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], het_out, quote=F, row.names=F) # print FID and IID for valid samples
q() # exit R


