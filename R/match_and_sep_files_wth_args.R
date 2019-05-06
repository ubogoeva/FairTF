file_prefix <- 'degs_m1500p1_'
##
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  file_prefix <- args
} else {
  stop("one argument must be supplied (input file).n", call.=FALSE)
}

all_peaks <- read.table('~/TranscInterTF/data/dapseq/all_m1500p1.bed', header = F, stringsAsFactors = F)
gene_set <- read.table('~/TranscInterTF/data/gens_enriched.txt', header = F, stringsAsFactors = F)[ ,1]
my_peaks <- all_peaks[match(gene_set, all_peaks$V5), ]
##match file with my coords

my_peaks <- na.omit(my_peaks)
write.table(my_peaks, file = paste0('~/TranscInterTF/data/dapseq/', file_prefix, 'all.bed'), sep = '\t', quote = F, col.names = F, row.names = F)

sep_files <- function(df, path, pref, chr) {
  chr_sub <- my_peaks[my_peaks$V1 == chr, ]
  write.table(chr_sub, file = paste0(path, pref, chr, '.bed'), sep = '\t', quote = F, col.names = F, row.names = F)
  return(nrow(chr_sub))
}
path <- "~/TranscInterTF/data/dapseq/"
sep_files(my_peaks, path, file_prefix, 'chr1')
sep_files(my_peaks, path, file_prefix, 'chr2')
sep_files(my_peaks, path, file_prefix, 'chr3')
sep_files(my_peaks, path, file_prefix, 'chr4')
sep_files(my_peaks, path, file_prefix, 'chr5')
write(paste('call ajo_chr_n_ds_2.bat', file_prefix, sep = ' '), file = paste0(path, 'run.bat'))
