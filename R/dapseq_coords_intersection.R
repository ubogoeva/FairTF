all_peaks <- read.table('~/R_projects/dapseq_new/all_m1500p1.bed', header = F, stringsAsFactors = F)
my_gens <- (read.table('~/TranscInterTF/data/dapseq/TF_vas.txt', header = F, stringsAsFactors = F))[, 1]
#files was separated previously, overlap15 files is in directory
#
#my_tf <- read.table('~/R_projects/dapseq_new/all_TF_568.txt', header = F, stringsAsFactors = F)[, 1] #old file
#forget about it

source('~/TranscInterTF/R/read_data_overlap.R')

all_TF <- read_data_overlap(working_dir = '~/TranscInterTF/data/dapseq/dapseq_pos/vas')
#274 files was combined
colnames(all_TF)[1] <- 'TF'
write.csv(all_TF, '~/TranscInterTF/data/dapseq/TF_coords_for_vasilina.csv', quote = F, row.names = F)
library(stringr)

TF_name <- str_split_fixed(all_TF$TF, '_', 3)
all_TF$TF <- TF_name[, 3]
all_TF$TF <- gsub('.overlap15', '', all_TF$TF)
my_tf <- read.table('~/R_projects/dapseq_new/all_TF_568_triv.txt', header = F, stringsAsFactors = F)
gens_tf <- data.frame(matrix(ncol = length(my_gens), nrow = length(my_tf[, 2])))
row.names(gens_tf) <- my_tf[, 2]
colnames(gens_tf) <- my_gens
i <- 1
j <- 1
while (i <= nrow(gens_tf)) {
  j <- 1
  while (j <= ncol(gens_tf)) {
    tf_subset <- all_TF[all_TF$TF == rownames(gens_tf)[i], 4]
    if (colnames(gens_tf)[j] %in% tf_subset) {
      gens_tf[i, j] <- 'yes'
    }
    else {
      gens_tf[i, j] <- ""
    }
    
    j <- j+1
  }
  
  i <- i+1
}
sum(rownames(gens_tf) == my_tf[, 2])
gens_tf$gene_id <- my_tf[, 1]
gens_tf <- gens_tf[, c(6, 1:5)]
write.csv(gens_tf, '~/TranscInterTF/data/dapseq/TF_intersect_Vasilina.csv', quote = F)

peak_for <- my_peaks[my_peaks$V4 == '+', ]
peak_rev <- my_peaks[my_peaks$V4 == '-', ]
head(all_TF)
intersect_for <- function(bed, all_TF) {
  df <- data.frame()
  for (gene in bed[, 5]) {
    sub_TF <- all_TF[all_TF[, 4] == gene, ]
    sub_TF$gene <- gene
    sub_TF$begin <- sub_TF[, 2] - bed[bed[, 5] == gene, 3]
    sub_TF$end <- sub_TF[, 3] - bed[bed[, 5] == gene, 3]
    df <- rbind(df, sub_TF)
  }
  return(df)
}

tf_for <- intersect_for(peak_for, all_TF)

intersect_rev <- function(bed, all_TF) {
  df <- data.frame()
  for (gene in bed[, 5]) {
    sub_TF <- all_TF[all_TF[, 4] == gene, ]
    sub_TF$gene <- gene
    sub_TF$begin <- -(sub_TF[, 3] - bed[bed[, 5] == gene, 2])
    sub_TF$end <-  -(sub_TF[, 2] - bed[bed[, 5] == gene, 2])
    df <- rbind(df, sub_TF)
  }
  return(df)
}
tf_rev <- intersect_rev(peak_rev, all_TF)
tf_raw <- rbind(tf_for, tf_rev)
write.csv(tf_raw, '~/TranscInterTF/data/dapseq/tf_raw_Vas.csv', row.names = F, quote = F)
tf_raw <- tf_raw[, c(1, 5:7)]
write.csv(tf_raw, '~/TranscInterTF/data/dapseq/TF_coords_Vasilina.csv', row.names = F, quote = F)
