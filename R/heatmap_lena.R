genes <- read.table('~/TranscInterTF/data/heatmap_lena.txt', header = T, stringsAsFactors = F)
whole_table <- read.csv('~/TranscInterTF/data/degs_annotated.csv', header = T, stringsAsFactors = F)
library(pheatmap)
sub_genes <- whole_table[match(genes$ID, whole_table$GeneID), ]
sub_genes[2, 3] <- log2(18.17/0.01)
sub_genes[16, 3] <- 0

rownames(sub_genes) <- genes$trivial
sub_genes$GeneID <- NULL
sub_genes <- as.matrix(sub_genes)
my_palette <- colorRampPalette(c("green", 'black', "red"))(n = 299)
pheatmap(sub_genes, color = my_palette, cluster_rows = F, cluster_cols = F)

heatmap.2(sub_genes, Rowv = F, Colv = F, trace='none', col=my_palette, scale="none", distfun = dist, na.rm = T, sepwidth=c(0.05,0.05), 
                 notecex=1.0, notecol="cyan", na.color=par("bg"), density.info = 'histogram', dendrogram = 'none')
##now lets count heatmaps for TF 
load('~/TranscInterTF/data/sign_datasets.RData')
tf_list <- read.table('~/TranscInterTF/data/heatmap_tf.txt', header = T, stringsAsFactors = F)
##only unique tf we should use
tf_list <- tf_list[!duplicated(tf_list$trivial), ]
tf_genes <- rnaseq_1h_sign[match(tf_list$ID, rnaseq_1h_sign$GeneID), ]
tf_genes$RNAseq_logFC_6h <- rnaseq_6h_sign[match(tf_list$ID, rnaseq_6h_sign$GeneID), 2]
tf_genes[, 4:10] <- lewis_sign[match(tf_list$ID, lewis_sign$GeneID), 2:8]
rownames(tf_genes) <- tf_list$trivial
write.csv(tf_genes, file = '~/TranscInterTF/data/TF_folds.csv', row.names = T, quote = F)
tf_genes$GeneID <- NULL
tf_genes <- as.matrix(tf_genes)

heatmap.2(tf_genes, Rowv = F, Colv = F, trace='none', col=my_palette, scale="none", distfun = dist, na.rm = T, sepwidth=c(0.05,0.05), 
          notecex=1.0, notecol="cyan", na.color=par("bg"), density.info = 'histogram', dendrogram = 'none')
