setwd('C:/Users/Лена/Documents/R_projects/araport')
agi <- read.table('C:/Users/Лена/Documents/флешка/gene_tf.txt', stringsAsFactors = F, header = F)
agi <- agi[!duplicated(agi$V1), ]
agi <- agi[order(agi$V1), ]
files <- dir(pattern = '*.csv')
##check that your files have AGI codes
files <- files[grep('AT[1-5]G*', files)]
level_heatmap <- c('light-grown seedling', 'dark-grown seedling', 'leaf', 'root',
                   'shoot apical meristem', 'root apical meristem')
##choose tissue to mean counts and for heatmap
mean_data <- function(file_name, level) {
  df <- read.csv(file_name)
  gene_mean <- unlist(strsplit(file_name, '_'))[1]
  i <- 2
  for (tissue in level) {
    #a <- mean(df[df$experiment_tissue == tissue])
    gene_mean[i] <- mean(df[df$experiment_tissue == tissue, 2])
    i <- i+1
  }
  return(gene_mean)
}
#define function to average counts in specified tissue
#this not beutiful, but creating dataset
temp_df <- mean_data(files[1], level_heatmap)
i <- 2
while (i <= length(files)) {
  temp_df <- rbind(temp_df, mean_data(files[i], level_heatmap))
  i <- i+1
}
genes_heat <- as.data.frame(temp_df, row.names = temp_df[, 1])
colnames(genes_heat) <- c('gene_id', level_heatmap)
rownames(genes_heat) == agi[, 1]
rownames(genes_heat) <- agi[, 2]

genes_heat <- genes_heat[order(rownames(genes_heat)), ]
write.csv(genes_heat, 'tissue_specific_expr.csv', row.names = T, quote = F)
library(gplots)
genes_heat <- genes_heat[, -1]
genes_heat[, 1] <- as.character(genes_heat[, 1])
genes_heat[, 2] <- as.character(genes_heat[, 2])
genes_heat[, 3] <- as.character(genes_heat[, 3])
genes_heat[, 4] <- as.character(genes_heat[, 4])
genes_heat[, 5] <- as.character(genes_heat[, 5])
genes_heat[, 6] <- as.character(genes_heat[, 6])
class(genes_heat[, 6])
genes_heat <- as.matrix(genes_heat)
my_palette <- colorRampPalette(c("black", "red"))(n = 299)
pl3 <- heatmap.2((genes_heat), Colv = 'none', trace='none', col=my_palette, scale="none", distfun = dist, na.rm = T, sepwidth=c(0.05,0.05), 
                 notecex=1.0, notecol="cyan", na.color=par("bg"), density.info = 'histogram')
