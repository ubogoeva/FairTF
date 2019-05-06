tf_from_dapseq <- read.table('~/TranscInterTF/data/dapseq/name_of_TF.txt', header = F, stringsAsFactors = F)
load('~/TranscInterTF/data/sign_datasets.RData')


#df_1h_subset <- make_deg(rnaseq_1h[match(tf_from_dapseq[, 1], rnaseq_1h[, 1]), ], colname = 'RNAseq_logFC_1h')


intersect_table <- cbind(tf_from_dapseq, rnaseq_1h_sign[match(tf_from_dapseq[, 1], rnaseq_1h_sign[, 1]), -1], 
                   rnaseq_6h_sign[match(tf_from_dapseq[, 1], rnaseq_6h_sign[, 1]), -1], 
                   lewis_sign[match(tf_from_dapseq[, 1], lewis_sign[, 1]), -1])
colnames(intersect_table)[1:4] <- c('GeneID', 'trivial_name', "RNAseq_logFC_1h", "RNAseq_logFC_6h")
head(intersect_table)
# combine dataset and get auxin sensitivity
for (i in 1:nrow(intersect_table)){
  intersect_table$aux_sens[i] <- any(intersect_table[i, -(1:2)] != 0, na.rm = T)
}
intersect_table <- intersect_table[, c(1:2, 12, 3:11)] #reorder columns

#get common table for overlapping TF vs whole 132 (130) genes
source('~/TranscInterTF/R/read_data_overlap.R')
# all_TF_old <- read_data_overlap(working_dir = '~/R_projects/rnaseq_chipseq/count_IAA_degs/ethylene_overlap/whole_chr')
all_TF_inter <- read_data_overlap(working_dir = '~/TranscInterTF/data/dapseq/dapseq_pos/whole_chr')
library(stringr)
#parse our TF files to get name from filename 
all_TF_new <- str_split_fixed(all_TF_inter$V1, '_', 3)
all_TF_inter$V1 <- all_TF_new[, 3]
all_TF_inter$V1 <- gsub('.overlap15', '', all_TF_inter$V1)
head(all_TF_inter)
rm(all_TF_new)
#load enriched genes from functional annotation
gene_set <- read.table('~/TranscInterTF/data/gens_enriched.txt', header = F, stringsAsFactors = F)[, 1]

nrow(intersect_table)
##find intersection for my TF
for (i in 1:nrow(intersect_table)) {
  intersect_table$intersection[i] <- paste(all_TF_inter[grep(intersect_table[i, 2], all_TF_inter$V1), 4], collapse = '/')
}
##возможно это можно сделать без цикла, использовав функцию match()
intersect_table <- intersect_table[, c(1:3, 13, 4:12)]
write.csv(intersect_table, file = '~/TranscInterTF/data/whole_aux_inter_TF.csv', row.names = F, quote = F)
write.csv(intersect_table[, 1:4], file = '~/TranscInterTF/data/short_aux_inter_TF.csv', row.names = F, quote = F)
