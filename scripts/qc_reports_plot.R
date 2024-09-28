
#Plots for QC report values:
library(ggplot2)

# Load input files:
lmiss <- read.table(snakemake@input[["lmiss"]], header = T)
imiss <- read.table(snakemake@input[["imiss"]], header = T)
maf <- read.table(snakemake@input[["maf"]], header = T)
hwe <- read.table(snakemake@input[["hwe"]], header = T)
het <- read.table(snakemake@input[["het"]], header = T)

#Define output
plots <- snakemake@output[["plots"]]

# Theme
my_theme <-  theme_classic() + 
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10, margin = margin(t = 12, r = 0, b = 0, l = 0), face = "bold"),
    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 12, b = 0, l = 0), face = "bold"),
    axis.title.y.right = element_text(size = 17, margin = margin(t = 0, r = 0, b = 0, l = 17)),
    plot.title = element_text(size = 12, color = "gray30", face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"))

# Open pdf
pdf(plots)

# MAF:
maf_05 <- nrow(maf[maf$MAF< 0.05, ])
maf_01 <- nrow(maf[maf$MAF< 0.01, ])
total <- nrow(maf)

ggplot(maf, aes(x=MAF)) +
  labs(title = "MAF distribution", y = "# of SNPs") +
  geom_histogram(fill = "#1e28b7", color = "black") +
  annotate("text", x = 0.2, y = 100000, label = paste0("MAF < 0.05: ",maf_05, " variants")) +
  annotate("text", x = 0.2, y = 95000, label = paste0("MAF < 0.01: ",maf_01, " variants")) +
  annotate("text", x = 0.2, y = 110000, label = paste0("Total: ",total, " variants")) + 
  my_theme

hist(maf$MAF, n=1000, main="MAF<0.05", xlim=c(0,0.05), xlab= "MAF", ylab= "# of SNPs", col= "#1e28b7") 

#Frequência de individual missingness
nimiss <- nrow(imiss[imiss$F_MISS > 0.01, ])

ggplot(imiss, aes(x=F_MISS)) +
  labs(title = "Individual missingness", y = "# of individuals", x = "Missing genotype rate" ) +
  geom_histogram(fill = "darkgreen", color = "black") +
  annotate("text", x = 0.1, y = 300, label = paste0("# Individuals (missing genotype rate > 0.01): ",nimiss)) +
  my_theme

hist(imiss$F_MISS, n=1000, main="Individual missingness <0.02", xlim=c(0,0.02), col= "darkgreen",  ylab= "# of individuals", xlab= "Missing genotype rate")


#Frequência de locus missingness
nlmiss <- nrow(lmiss[lmiss$F_MISS > 0.01, ])

ggplot(lmiss, aes(x=F_MISS)) +
  labs(title = "Histogram SNP missingness", y = "# of SNPs", x = "Missing individual rate" ) +
  geom_histogram(fill = "purple", color = "black") +
  annotate("text", x = 0.6, y = 200000, label = paste0("# SNPs (missing rate > 0.01): ",nlmiss)) +
  my_theme

hist(lmiss$F_MISS, n=100, main="Histogram SNP missingness <0.1", xlim=c(0,0.1), col= "purple", xlab= "Missing individual rate", ylab= "# of SNPs")

#P-valor associado ao HWE
hist(-log10(hwe$P), n=1000, main="Hardy-Weinberg Equilibrium p>10^-5", xlim=c(0,5), col= "lightpink4", xlab = "-log10(P)", ylab="# of SNPs")

# Heterozygosity
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
dev.off()
