
source("~/R_projects/FSGOR/func_annot_functions.r")
onto <- read.obo2("~/R_projects/FSGOR/go.obo")
annot <- read.gaf2("~/R_projects/FSGOR/gene_association.tair")
pre_annot <- prepare_annotation(onto, annot$data, bg)
annot_ethylene <- function(df, onto, annot){
  df <- df
  bg <- as.character(df[, 1]) 
  source("~/R_projects/my_functions/select_degs.r")
  degs <- select_degs(df, log2(1.5))
  print(nrow(degs)) ##Столько дэг получилось"
  res <-  FuncAnnotTest(pre_annot, degs$ID, onto$real_names, 
                            onto$namespaces)
  eth_all <<- res[grep('\\bethylene\\b', res$name), ]
  qits <- unlist(strsplit(eth_all$qit_genes, "/", fixed = T))
  qit_uniq <- unique(qits) ##отбор уникальных дэг
  print(length(qits)) ##количество не уникальных ДЭГ
  return(qit_uniq)
}

