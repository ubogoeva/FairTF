lewis <- read.table('~/R_projects/Lewis_all/lewis_combine.txt', header = T, stringsAsFactors = F)
my_tf <- read.table('~/R_projects/dapseq_new/all_TF_568.txt', header = F, stringsAsFactors = F)[, 1]

tf_lewis <- lewis[match(my_tf, lewis$ID), ]
tf_lewis <- na.omit(tf_lewis)
write.csv(tf_lewis, '~/R_projects/rnaseq_chipseq/lewis_tf.csv', row.names = F)

read_data  <- function(working_dir){
  old_dir <- getwd()
  setwd(working_dir)
  df  <- data.frame()
  number  <- 0
  for (i in dir(pattern = "*.overlap15")){
    temp_df  <- read.table(i, stringsAsFactors = F)
    temp_df$V1 <- i
    df  <- rbind(df, temp_df)
    number <- number + 1
  }
  print(paste(as.character(number), "files were combined"))
  setwd(old_dir)
  return(df)
}
all_TF <- read_data(working_dir = '~/R_projects/rnaseq_chipseq/count_IAA_degs/ethylene_overlap/whole_chr')
a <- read.table("acc_m1500p1_ZFHD_ATHB24_col_a.overlap15", stringsAsFactors = F)

write.csv(all_TF, 'combine_TF.csv', row.names = F)
library(stringr)
all_TF_new <- str_split_fixed(all_TF$V1, '_', 3)
all_TF$V1 <- all_TF_new[, 3]
all_TF_new <- str_split_fixed(all_TF$V1, '.', 1)
all_TF$V1 <- gsub('.overlap15', '', all_TF$V1)
head(all_TF)
write.csv(all_TF, 'combine_TF_upd.csv', row.names = F)
gene_set <- read.table('~/TranscInterTF/data/gens_enriched.txt', header = F, stringsAsFactors = F)[, 1]
my_tf <- read.table('~/R_projects/dapseq_new/all_TF_568_triv.txt', header = F, stringsAsFactors = F)[, 2]
gens_tf <- data.frame(matrix(ncol = length(agi), nrow = length(my_tf)))
row.names(gens_tf) <- my_tf
colnames(gens_tf) <- agi
#row.names = dir(pattern = '*.overlap')
head(gens_tf)
# 'AT1G19220' %in% all_TF[all_TF$V1 == 'ZFHD_ATHB34_colamp_a', 4]


# i <- 1
# j <- 1
# while (i <= 2) {
#   while (j <= 2)
#     if (colnames(gens_tf)[as.numeric(j)] %in% all_TF[all_TF$V1 == rownames(gens_tf)[as.numeric(i)], 4]) {
#       gens_tf[as.numeric(i), as.numeric(j)] <- T
#     }
#     else {
#       gens_tf[as.numeric(i), as.numeric(j)] <- F
#     }
#   i <- i+1
#   j <- j+1
# }
# i <- 1
# j <- 1
#try to intersect
#gens_tf <- rownames(gens_tf)
#gens_tf <- as.data.frame(gens_tf)
i <- 1
while (i <= nrow(gens_tf)) {
  j <- 1
  while (j <= ncol(gens_tf)) {
    tf_subset <- all_TF[all_TF$V1 == rownames(gens_tf)[i], 4]
    if (colnames(gens_tf)[j] %in% tf_subset) {
      gens_tf[i, j] <- colnames(gens_tf)[j]
    }
    else {
      gens_tf[i, j] <- NA
    }
    
    j <- j+1
  }

  i <- i+1
}
gens_tf$summary <- gens_tf[]
length(gens_tf[2, is.na(gens_tf[2, ])])
head(gens_tf)
write.csv(gens_tf, 'tf_final.csv', quote = F)






