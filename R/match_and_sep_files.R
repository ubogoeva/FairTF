all_peaks <- read.table('~/TranscInterTF/data/dapseq/all_m1500p1.bed', header = F, stringsAsFactors = F)
gene_set <- read.table('~/TranscInterTF/data/dapseq/TF_vas.txt', header = F, stringsAsFactors = F)[ ,1]
my_peaks <- all_peaks[match(gene_set, all_peaks$V5), ]
##match file with my coords
my_peaks <- na.omit(my_peaks)
write.table(my_peaks, file = '~/TranscInterTF/data/dapseq/vas_m1500p1_all.bed', sep = '\t', quote = F, col.names = F, row.names = F)

sep_files <- function(df, path, pref, chr) {
  chr_sub <- my_peaks[my_peaks$V1 == chr, ]
  write.table(chr_sub, file = paste0(path, pref, chr, '.bed'), sep = '\t', quote = F, col.names = F, row.names = F)
  return(nrow(chr_sub))
}
path <- "~/TranscInterTF/data/dapseq/"
prefix <- 'vas_m1500p1_'
sep_files(my_peaks, path, prefix, 'chr1')
sep_files(my_peaks, path, prefix, 'chr2')
sep_files(my_peaks, path, prefix, 'chr3')
sep_files(my_peaks, path, prefix, 'chr4')
sep_files(my_peaks, path, prefix, 'chr5')

