print('This function is for subset up-degs or down-degs for FoldGO')
select_degs_up_down <- function(dataset, lg=log2(1.5), direction=c('up', 'down', 'both')) {
  logFC <- dataset[, 2]
  FDR <- dataset[, 4]
  if (direction == 'up') {
    deg_table <- subset(dataset, dataset[, 4] < 0.05 & (dataset[, 2] > lg))
    return(deg_table)
  }
  if (direction == 'down') {
    deg_table <- subset(dataset, FDR < 0.05 & (logFC < -lg))
    return(deg_table)
  }
  else {
      deg_table <- subset(dataset, FDR < 0.05 & (logFC > lg | logFC < -lg))
      return(deg_table)
  }
}

