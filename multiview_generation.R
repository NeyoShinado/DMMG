multiview_generation <- function(X, droprate, clust, sigmac, sigmag, cell_thre=20, gene_thre=50){
  
  # init 
  N = dim(droprate)[1]
  P = dim(droprate)[2]
  map = unique(clust)
  K = length(map)
  split_data = list()
  Dg = list()
  Dc = list()
  view_id = list()
  index = 1    # index of split_data
  
  meandrop = apply(droprate, 2, mean)
  cat("### constructing multiview_data.../n")
  for(i in c(1:K)){
    id = which(clust == map[i])
    if(length(id) < cell_thre){
      next
    }else{
      clustdrop = droprate[id, ]
      clust_meandrop = apply(clustdrop, 2, mean)    # gene meandrop of select clust
      
      # lowdrop data generate
      lowdrop = which(clust_meandrop < 0.1)
      if(length(lowdrop) > gene_thre){
        split_data[[index]] = X[, lowdrop]
        #* need resive by group gene's total drop
        Dg[[index]] = Matrix(diag(exp(-sigmag * meandrop[lowdrop])), sparse=TRUE)
        Dc[[index]] = Matrix(diag(exp(-sigmac * apply(droprate[, lowdrop], 1, mean))), sparse=TRUE)
        view_id[[index]] = lowdrop
        index = index + 1
      }else{
        next
      }
      #* meandrop > 0.8 且频数高于其余90%基因的基因集，即间隔最短的十分位
      # highdrop data generate
      dropintev = quantile(clust_meandrop[setdiff(1:P, lowdrop)], seq(0, 10)*0.1)    #* 十等分
      dropgap = sapply(1:(length(dropintev)-1), function(i){
        dropintev[i+1]-dropintev[i]
      })
      hdropl = dropintev[which.min(dropgap)]
      hdropr = dropintev[which.min(dropgap)+1]
      highdrop = which(clust_meandrop >= max(0.8, hdropl) & clust_meandrop < hdropr)
      if(length(highdrop) > gene_thre){
        split_data[[index]] = droprate[, highdrop]
        #* drop weight of highdrop
        Dg[[index]] = Matrix(diag(exp(sigmag * (meandrop[highdrop]-1))), sparse=TRUE)
        Dc[[index]] = Matrix(diag(exp(sigmac * (apply(droprate[, highdrop], 1, mean)-1))), sparse=TRUE)
        view_id[[index]] = highdrop
        index = index + 1
      }else{
        next
      }
    }
  }
  
  res = list(split_data = split_data, Dg = Dg, Dc = Dc, view_id = view_id)
}