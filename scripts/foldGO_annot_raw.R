library(FoldGO)
library(org.At.tair.db)
library(topGO)
library(mgsa)

data("example")
gaf <- readGAF("~/R_projects/FSGOR/tair.gaf")
setwd('~/R_projects/FSGOR')
star <- read.csv('~/R_projects/RNA-seq/deseq_result_upd_sample.csv', header = T, stringsAsFactors = F)
gaf_path <- system.file("extdata", "~/R_projects/FSGOR/tair.gaf", package = "FoldGO")
gaf <- GAFReader(file = "~/R_projects/FSGOR/tair.gaf",  geneid_col = 10)
source('~/R_projects/my_functions/select_degs.R')
star <- star[, c(1, 3, 6, 7)]
colnames(star) <- c('GeneID',	'FC',	'pval',	'qval')
degs_star <- select_degs(star, log2(1.5))
degs_up <- degs_1h[degs_1h$log2FoldChange > 0, ]
degs_down <- degs_1h[degs_1h$log2FoldChange < 0, ]
colnames(degs_up) <- c('GeneID',	'FC',	'pval',	'qval')
colnames(degs_down) <- c('GeneID',	'FC',	'pval',	'qval')
up_groups <- GeneGroups(degs_up, 6)

down_groups <- GeneGroups(degs_down, 6)
bbgenes <- bg

up_annotobj <- FuncAnnotGroupsTopGO(genegroups = up_groups, namespace = "CC", 
                                    customAnnot = gaf, annot = topGO::annFUN.GO2genes, bggenes = bggenes)
down_annotobj <- FuncAnnotGroupsTopGO(genegroups = down_groups, namespace = "BP", 
                                    customAnnot = gaf, annot = topGO::annFUN.GO2genes, bggenes = bggenes)

res <- getResultList(up_annotobj)
res_down <- getResultList(down_annotobj)
res_go <- res$`1-2`
down_table <- res_down$`1-2`
ethylene <- res_go[grep('\\bethylene\\b', res_go$name), ]
a <- down_table[grep('\\bethylene\\b', down_table$name), ]
ethylene <- rbind(ethylene, a)
annot <- getAnnotation(gaf)
vector_go <- c()
for (go in (ethylene$GO_id)) {
  vector_go <- c(vector_go, annot[[go]])
}
vector_go <- unique(vector_go)
gens_enrich <- intersect(vector_go, degs_star$AGI)
