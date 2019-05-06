tf_from_dapseq <- read.table('~/TranscInterTF/data/dapseq/name_of_TF.txt', header = F, stringsAsFactors = F)
load('~/TranscInterTF/data/raw_datasets.RData')
combined_data <- read.csv('~/TranscInterTF/data/degs_annotated.csv', header = T, stringsAsFactors = F)
source('~/TranscInterTF/R/make_deg.R')
#rnaseq_1h <- make_deg(rnaseq_1h[, c(1, 3, 7)], yname = 'RNAseq_logFC_1h')
make_deg <- function(my_data, yname='new_log') {
  #my_data <- my_data[, c(1, 2, 3)]
  
  colnames(my_data) <- c('GeneID', 'logFC', 'FDR')
  new_gens <- as.data.frame(my_data$GeneID)
  for (i in 1:nrow(my_data)) {
    if (is.na(my_data[i, 3])) { next
      cat(n)
    }
    if (my_data[i, 3] < 0.05 & (my_data$logFC[i] > log2(1.5) | my_data$logFC[i] < -log2(1.5))) {
      new_gens[i, 2] <- my_data$logFC[i]
    }
    else {
      new_gens[i, 2] <- 0
    }
  }
  colnames(new_gens) <- c('GeneID', yname)
  return(new_gens)
}

df_1h_subset <- make_deg(rnaseq_1h[match(tf_from_dapseq[, 1], rnaseq_1h[, 1]), ], yname = 'RNAseq_logFC_1h')
df_6h_subset <- make_deg(rnaseq_6h[match(tf_from_dapseq[, 1], rnaseq_6h[, 1]), ], yname = 'RNAseq_logFC_6h')
df_lewis_subset <- make_deg(lewis_data[match(tf_from_dapseq[, 1], lewis_data[, 1]), ], yname = 'logFC_05h')

sum_table <- cbind(df_1h_subset$X, df_6h_subset$ID, df_lewis_subset$ID, tf_from_dapseq[, 1])

df_1h_subset[is.na(df_1h_subset$X), ] <- 0
df_1h_subset[is.numeric(df_1h_subset$X), ]
base::setdiff(gene_set, df_1h_subset$X)

subset_1h <- make_deg(df_1h_subset[, c(1, 3, 7)], yname = 'RNAseq_logFC_1h')