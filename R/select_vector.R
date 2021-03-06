##������ �������, �� ������� ����� �������� ��������
#df <- read.table('~/R_projects/microarrays/peri_bargmann.txt', header = T, stringsAsFactors = F) 
#gens <- df$ID[110:220] ##��� ������ ��������, ������� �� ����� �������� � ������� ������� df
#res <- data.frame()
##������� ��� ������ �������� �� �������
select_vector <- function(df, vect) {
  vect <- as.character(vect)
  res <- data.frame()
  for (gene in (vect)) {
    res <- rbind(res, df[grep(gene, df[, 1]), ])
  }
  return(res)
}
##������� ��� ������ �������� ���������������
select_vector_invert <- function(df, vect) {
  res <- data.frame()
  for (i in 1:length(vect)) {
    b <- vect[i]
    res <- rbind(res, df[grep(b, df$ID, invert = T), ])
  }
  return(res)
}
