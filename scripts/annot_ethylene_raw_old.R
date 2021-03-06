
source("~/TranscInterTF/scripts/func_annot_functions.r")

annot_ethylene <- function(df, onto, annot){
  df <- df
  bg <- as.character(df[, 1]) 
  source("~/R_projects/my_functions/select_degs.r")
  degs <- select_degs(df, log2(1.5))
  print(nrow(degs)) ##������� ��� ����������"
  res <-  FuncAnnotTest(pre_annot, degs$ID, onto$real_names, 
                            onto$namespaces)
  eth_all <<- res[grep('\\bethylene\\b', res$name), ]
  qits <- unlist(strsplit(eth_all$qit_genes, "/", fixed = T))
  qit_uniq <- unique(qits) ##����� ���������� ���
  print(length(qits)) ##���������� �� ���������� ���
  return(qit_uniq)
}

