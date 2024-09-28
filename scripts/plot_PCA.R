library(ggplot2)
library(data.table)

# Input and output
data <- fread(snakemake@input[["eigenvec"]], header = T)
plot_pca <- snakemake@output[["plot_pca"]]

# My theme
my_theme <-  theme_classic() + 
  theme(text = element_text(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_text(size = 16, margin = margin(t = 12, r = 0, b = 0, l = 0), face = "bold"),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0), face = "bold", hjust = 0.75),
        plot.title = element_text(size = 16, color = "gray30", hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"))

# Plot
p <- ggplot(data, aes(x = PC1, y = PC2)) +  
  geom_point(aes(), lwd = 0.8) +
  my_theme +
  theme(legend.justification = c(0, 0), legend.position = "bottom", 
        legend.box.background = element_rect(color="grey", size=0.5), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) 

ggsave(plot_pca, p, width = 17, height= 10)

