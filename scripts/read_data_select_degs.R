##read data
# setwd('~/TranscInterTF/') 
#not require
rnaseq_1h <- read.csv('~/TranscInterTF/data/1h_RNAseq_star_deseq2.csv', header = T, stringsAsFactors = F)
rnaseq_6h <- read.table('~/TranscInterTF/data/6h_RNAseq.txt',header =T, stringsAsFactors = F)
lewis_data <- read.table('~/TranscInterTF/data/lewis.txt', header = T, stringsAsFactors = F)

save(rnaseq_1h, rnaseq_6h, lewis_data, file = '~/TranscInterTF/data/raw_datasets.RData')
##make degs for saving and downstream work
source('~/TranscInterTF/R/make_deg.R')
lewis_sign <- cbind(make_deg(lewis_data[, c(1,2,3)], colname = 'logFC_05h'), make_deg(lewis_data[, c(1,4,5)], colname = 'logFC_1h')[, 2],
                     make_deg(lewis_data[, c(1,6,7)], colname = 'logFC_2h')[, 2], make_deg(lewis_data[, c(1,8,9)], colname = 'logFC_4h')[, 2],
                     make_deg(lewis_data[, c(1,10,11)], colname = 'logFC_8h')[, 2], make_deg(lewis_data[, c(1,12,13)], colname = 'logFC_12h')[, 2],
                     make_deg(lewis_data[, c(1,14,15)], colname = 'logFC_24h')[, 2])
colnames(lewis_sign) <- c('GeneID', 'logFC_05h', 'logFC1h', 'logFC2h', 'logFC4h', 'logFC8h', 'logFC12h', 'logFC24h')
rnaseq_1h_sign <- make_deg(rnaseq_1h[, c(1, 3, 7)], 'RNAseq_logFC_1h', lg = 0)
rnaseq_6h_sign <- make_deg(rnaseq_6h[, c(1, 5, 7)], 'RNAseq_logFC_6h', lg = 0)
save(rnaseq_1h_sign, rnaseq_6h_sign, lewis_sign, file = '~/TranscInterTF/data/sign_datasets.RData')
# source('~/TranscInterTF/scripts/select_degs.R')
degs_1h <- select_degs(rnaseq_1h[, c(1, 3, 6, 7)], lg=log2(1.5), FDR_thresh = 0.05)
degs_6h <- select_degs(rnaseq_6h[, c(1, 5, 6, 7)], lg=log2(1.5), FDR_thresh = 0.05)
#check our 6h dataset
nrow(degs_6h[degs_6h$logFC == 'Inf', ]) #21
#we can find value Inf in column logFC 
#filter these values
degs_6h <- degs_6h[degs_6h$logFC != 'Inf', ]
##get lewis degs

##background genes for functional annotation
bg <- rnaseq_1h$X

##make functional annotation
#source("~/TranscInterTF/scripts/func_annot_functions.r")

#onto <- read.obo2("~/TranscInterTF/data/go.obo")
#annot <- read.gaf2("~/TranscInterTF/data/tair.gaf")
#pre_annot <- prepare_annotation(onto, annot$data, bg)
source('~/R_projects/my_functions/annot_ethylene.R')
gaf <- read.table("~/R_projects/FSGOR/tair.gaf")
