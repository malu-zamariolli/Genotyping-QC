# Make a list of individuals after excluding sex mismatchs

#Snakemake input

valid_samples = snakemake@input[["valid"]]
sexcheck = snakemake@input[["sexcheck"]]
plot = snakemake@output[["plot"]]

gender <- read.table(sexcheck, header=T,as.is=T)

pdf(plot)
hist(gender[,6],main="Gender", xlab="F")

male=subset(gender, gender$PEDSEX==1)
hist(male[,6],main="Men",xlab="F")

female=subset(gender, gender$PEDSEX==2)
hist(female[,6],main="Women",xlab="F")
dev.off()

# Create list of individuals that pass check-sex
# Read in file
valid <- read.table(valid_samples, header=T)
dat <- read.table(sexcheck, header=T)
valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], snakemake@output[["filtered_sex"]], row.names=F, col.names=F, sep="\t", quote=F) 
q() # exit R
