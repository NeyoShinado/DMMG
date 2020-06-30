eva <- function(path="result/pars/", pattern, col=1){
  files = list.files(path, pattern=pattern)
  nmi = matrix(0, 3, 3)
  colnames(nmi) = c("l2_1", "l2_10", "l2_100")
  rownames(nmi) = c("l3_1", "l3_10", "l3_100")
  for(i in 1:length(files)){
    res = readRDS(paste0(path, files[i]))
    nmi[i] = res$nmi[dim(res$nmi)[1], col]
  }
  print(t(nmi))
}