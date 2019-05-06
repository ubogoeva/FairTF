#gene set we should get from sunctional annotation with combining gens and filter duplicated
#but I wait for annotation gaf and load it from my old gene set

gene_set <- read.table('~/TranscInterTF/data/gens_enriched.txt', header = F, stringsAsFactors = F)[, 1]
#then we should get output as table with significant values of logExpression from each dataset
source('~/TranscInterTF/R/select_vector.R')
source('~/TranscInterTF/R/make_deg.R')
load('~/TranscInterTF/data/raw_datasets.RData')
#source('~/TranscInterTF/R/select_degs_up_down.R')
df_1h_subset <- rnaseq_1h[match(gene_set, rnaseq_1h[, 1]), ]
df_1h_subset[is.na(df_1h_subset$X), ] <- 0
df_1h_subset[is.numeric(df_1h_subset$X), ]
  base::setdiff(gene_set, df_1h_subset$X)

subset_1h <- make_deg(df_1h_subset[, c(1, 3, 7)], yname = 'RNAseq_logFC_1h')
df_6h_subset <- select_vector(rnaseq_6h, gene_set) 
diff_6h <- base::setdiff(gene_set$V1, df_6h_subset$ID) #"AT2G47520" "AT5G55620" "AT5G07310" "AT1G28160" "AT3G59060" "AT1G71520" "AT2G25820" "AT1G22810"
#its for inputation for bind DFs
subset_6h <- make_deg(df_6h_subset[, c(1, 5, 7)], yname = 'RNAseq_logFC_6h')
##GET DEGS from LEWIS dataset

lewis_all <- cbind(make_deg(select_vector(lewis_data[, c(1,2,3)], gene_set$V1), yname = 'logFC_05h'), 
                   make_deg(select_vector(lewis_data[, c(1,4,5)], gene_set$V1), yname = 'logFC_1h')[, 2], 
                   make_deg(select_vector(lewis_data[, c(1,6,7)], gene_set$V1), yname = 'logFC_2h')[, 2],
                   make_deg(select_vector(lewis_data[, c(1,8,9)], gene_set$V1), yname = 'logFC_4h')[, 2],
                   make_deg(select_vector(lewis_data[, c(1,10,11)], gene_set$V1), yname = 'logFC_8h')[, 2],
                   make_deg(select_vector(lewis_data[, c(1,12,13)], gene_set$V1), yname = 'logFC_12h')[, 2],
                   make_deg(select_vector(lewis_data[, c(1,14,15)], gene_set$V1), yname = 'logFC_24h')[, 2])
#subset lewid dataframe with only significant values of logFC (FDR < 0.05)


head(lewis_all)
colnames(lewis_all) <- c('AGI', 'logFC_05h', 'logFC_1h', 'logFC_2h', 'logFC_4h', 'logFC_8h', 'logFC_12h', 'logFC_24h')
head(lewis_all)
diff_lewis <- setdiff(gene_set$V1, lewis_all$AGI)


##its okay, now we have to combine three datasets to one
#we cant make cbind because we have different nrow


subset_1h$GeneID <- as.character(subset_1h$GeneID)
subset_1h[132, 1] <- "AT1G22810"
tail(subset_1h)

#realize appending to 1h_RNaseq
old_row <- nrow(subset_6h)
subset_6h$GeneID <- as.character(subset_6h$GeneID)
subset_6h[(old_row+1):(old_row+length(diff_6h)), 1] <- as.character(diff_6h)
lewis_all$AGI <- as.character(lewis_all$AGI)
lewis_all[131:132, 1] <- diff_lewis
#now lets cBIND our datasets
combined_data <- cbind(subset_1h[order(subset_1h$GeneID), ], subset_6h[order(subset_6h$GeneID), 2], lewis_all[order(lewis_all$AGI), -1])
colnames(combined_data)[3] <- 'RNAseq_logFC_6h'
head(combined_data)
write.csv(combined_data, file = '~/TranscInterTF/data/degs_annotated.csv', quote = F, row.names = F)
