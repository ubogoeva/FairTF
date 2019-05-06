setwd('~/TranscInterTF/data/araport/')
aco <- read.table('~/TranscInterTF/data/araport/aco1.txt', stringsAsFactors = F, header = T)
if (nrow(aco) != nrow(unique(aco))) {
  aco <- aco[!duplicated(aco[, 1]), ]
}
aco <- aco[order(aco$gene_id), ]
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

aco_file <- vector()
i <- 1
for (gene in aco$gene_id) {
  aco_file <- c(aco_file, files[grep(aco$gene_id[i], files)])
  i <- i+1
}
#choose only target files for aco
#this not beutiful, but creating dataset
temp_df <- mean_data(aco_file[1], level_heatmap)
i <- 2
while (i <= length(aco_file)) {
  temp_df <- rbind((temp_df), mean_data(aco_file[i], level_heatmap))
  i <- i+1
}
genes_heat <- as.data.frame(temp_df, row.names = temp_df[, 1], stringsAsFactors = F)
colnames(genes_heat) <- c('gene_id', level_heatmap)

(rownames(genes_heat) == aco[, 1])
rownames(genes_heat) <- aco[, 2]

genes_heat <- genes_heat[order(rownames(genes_heat)), ]
write.csv(genes_heat, 'tissue_specific_expr_aco1.csv', row.names = T, quote = F)
library(gplots)

library(pheatmap)

genes_heat <- genes_heat[, -1]

class(genes_heat[, 1])
#genes_heat <- as.matrix(genes_heat)
genes_heat[, 1] <- as.numeric(as.character(genes_heat[, 1]))
genes_heat[, 2] <- as.numeric(as.character(genes_heat[, 2]))
genes_heat[, 3] <- as.numeric(as.character(genes_heat[, 3]))
genes_heat[, 4] <- as.numeric(as.character(genes_heat[, 4]))
genes_heat[, 5] <- as.numeric(as.character(genes_heat[, 5]))
genes_heat[, 6] <- as.numeric(as.character(genes_heat[, 6]))

my_palette <- colorRampPalette(c("black", "red"))(n = 299)
my_array <- as.matrix(genes_heat)
pheatmap(my_array, color = my_palette, cluster_rows = F, cluster_cols = F)
#heatmap.2(my_array, Rowv = F, Colv = F, trace='none', col=my_palette, scale="none", distfun = dist, na.rm = T, sepwidth=c(0.05,0.05), 
#          notecex=1.0, notecol = 'red4', na.color=par("bg"), density.info = 'histogram', dendrogram = 'none')
