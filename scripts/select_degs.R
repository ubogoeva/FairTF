#function read table, with logFC as second column and FDR as fourth
print('function read table, with logFC as second column and FDR as fourth')
print('you can choose value of log and FDR for selecting ')
select_degs <- function(dataset, lg=log2(1.5), FDR_thresh=0.05){
  logFC <- dataset[, 2]
  FDR <- dataset[, 4]
  deg_table <- subset(dataset, FDR < 0.05 & (logFC > lg | 
                                               logFC < -lg))
  return(deg_table)
}
