# clustering vars updating
library(Matrix)
library(bignmf)


var_update <- function(lg_X, K, npc, S, neighbors, guide_cluster, gt, mu_g, sd_g, rep=1, sigmac=0.5, sigmag=0.5, lambda1=0.1, lambda2=0.1, 
                       thre_J=1e-4, drop_thre=0.5, iteration=1, 
                       clust_iteration=100, imp_iteration=100, 
                       output="result/", res_save=TRUE, data_name){
  s = Sys.time()
  # lg_X without pseudo count
  N = dim(lg_X)[1]
  P = dim(lg_X)[2]
  
  
  imp_iter = 1
  droprate = 0
  imp_X = 0
  
  
  ### imputing_update
  cat(paste0("## carrying ", imp_iter, "th cluster imputation...\n\n"))
  #*clust = clust_set[[i]]
  imp_res = imputing_update(log10(10^(lg_X) + 0.01), guide_cluster, neighbors, 
                            imp_iteration=imp_iteration, 
                            out_file=paste0(output, "imp_res/", imp_iter, "th_impres.rds"), 
                            res_save=res_save)
  
  #imp_X = log10(10^(imp_X) - 0.01)   # log-tran 1.01 for imputing, 1.0 for clustering
  
  ## arrange imp_output
  #local_sim = imp_res$local_sim
  imp_X = imp_res$imp_X
  droprate = imp_res$droprate
  meandrop = apply(droprate, 2, mean)   # 0 - 0.9  gene meandrop
  meangene = apply(lg_X, 2, mean)
  
  # fea separation
  meandrop_id = order(meandrop, decreasing=FALSE)
  dropmean = cbind(meandrop, meangene)
  v = svd(dropmean)$v
  gene_group = dropmean %*% v[, 1]
  
  # locate center
  gene_class = kmeans(gene_group, centers=sample(gene_group[meandrop_id[1:100]], K), nstart = 20, iter.max = 1000)$cluster
  
  
  if(FALSE){
  meandrop = apply(droprate, 1, mean)
  meangene = apply(lg_X, 1, mean)
  id = order(meandrop, decreasing=TRUE)
  png(paste0("../MGGE/result/figure/cellmeandrop verus meangene/cell\'s meandrop ", 
             data_name, ".png"), width=600, height=600)
  plot(meandrop[id])
  dev.off()
  png(paste0("../MGGE/result/figure/cellmeandrop verus meangene/cell\'s meangene ", 
             data_name, ".png"), width=600, height=600)
  plot(meangene[id])
  dev.off()
  png(paste0("../MGGE/result/figure/cellmeandrop verus meangene/cell\'s grouped meangene of", 
             data_name, ".png"), width=600, height=600)
  data = data.frame(id = 1:N, gene = meangene[id], class=as.character(gt[id]))
  p = ggplot(data, aes(x=id, y=gene, color=class)) + geom_point()
  plot(p)
  dev.off()
  }
  
  if(FALSE){
  # 16\gene normalize
  split_data = lapply(split_data, function(i){
    gene_norm = sqrt(colSums(i^2))
    i = i / t(matrix(gene_norm, dim(i)[2], N))
    return(i)
  })
  }
    
  
  # divide gene
  #####
  map = unique(gene_class)
  view_id = lapply(1:length(map), function(i){
    id = which(gene_class == map[i])
    if(length(id) < K){
      
    }else{
      return(id)
    }
  })
  view_id = view_id[sapply(view_id, function(i){all(!is.null(i))})]
  
  meandrop = sapply(view_id, function(i){
    res = sum(meandrop[i]) / length(i)
    return(res)
  })
  
  #* 8\ log_X or imp_X
  split_data = lapply(view_id, function(i){
    res = droprate[, i]  
    return(res)
  })
  
  
  #* 12\weighted cell
  Dc = rowSums(droprate)
  Dg = colSums(droprate)
  Dg = lapply(view_id, function(i){
    res = Dg[i]
    return(res)
  })
  Dc = Matrix(diag(exp(-sigmac * Dc / P)), sparse=TRUE)
  Dg = lapply(1:length(view_id), function(i){
    Matrix(diag(exp(-2 * meandrop[i] * Dg[[i]] / N)), sparse=TRUE)    # scale factor 2 for small dataset and 3 for the big one
  })
  
  # init vars
  set.seed(rep)   #*15\ produce reproductive result
  #H = matrix(runif(dim(lg_X)[1]*K), N, K)
  H = eigen(diag(colSums(S)) - S)$vector
  H = H[, (N-K):(N-1)]
  H[H<=0] = 1e-10
  W = matrix(runif(N * npc), N, npc)
  Wv = lapply(1:length(view_id), function(i){
    pi = dim(split_data[[i]])[2]
    Wi = matrix(runif(N * npc), N, npc)
    return(Wi)
  })
  V = lapply(1:length(view_id), function(i){
    pi = dim(split_data[[i]])[2]
    Vi = matrix(runif(npc * pi), npc, pi)
    return(Vi)
  })
  

  ### Clustering update
  cat(paste0("## running ", imp_iter, "th vars update via clustering and imputation...\n"))
  clust_res = clustering_update(split_data, K, npc, lambda1, lambda2, W=W, V=V, Wv=Wv, 
                                H=H, Dc=Dc, Dg = Dg, S=S, iteration=clust_iteration, meandrop=meandrop,
                                out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), 
                                res_save=res_save)
  #clust_res = clustering_update(imp_X, K, npc, lambda1, lambda2, W=W, V=V, H=H, Dc=Dc, L=L, iteration=clust_iteration,
  #                             out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), res_save=res_save)

  H = clust_res$H
  V = clust_res$V
  W = clust_res$W
  Wv = clust_res$Wv
  dw = clust_res$dw
  alpha = clust_res$alpha
  cluster = clust_res$cluster
  nmi = clust_res$nmi
  res_J = clust_res$J
  lambda1 = clust_res$lambda1
  lambda2 = clust_res$lambda2
  J_DR = clust_res$J_DR
  J_HE = clust_res$J_HE
  J_LE = clust_res$J_LE
  sigma = clust_res$sigma
  #neighbors = clust_res$neighbors
  
  cat("# MGGE iteration complete!\n")
  time = Sys.time() - s
  message("## Consume", time, "seconds.\n")
  
  res = list(cluster=cluster, imp_X = imp_X, guide_cluster = guide_cluster, lambda1 = lambda1, lambda2 = lambda2,
             droprate=droprate, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, J = res_J, nmi = nmi, sigma = sigma, beta = beta, 
             W = W, V = V, Wv=Wv, H = H, Dc = Dc, Dg = Dg, drop_id = view_id, S = S, dw = dw, weight = alpha, time = time)  

  return(res)

}

imputing_update <- function(lg_X, cluster, neighbors, drop_thre=0.5, out_file="result/imp_res/res.rds", ncores=1, imp_iteration=100, res_save=TRUE){
  # return log-tran X
  # lg 1.01
  # not parallel version on imputation
  imp_res = imputation_wlabel_model(lg_X, cluster, neighbors=neighbors, point=log10(1.01),
                                    drop_thre=drop_thre, ncores=ncores, imp_iteration=imp_iteration)
  imp_X = imp_res$count_imp
  droprate = imp_res$droprate
  local_sim = imp_res$local_sim
  if(res_save){
    saveRDS(imp_res, out_file)
  }
  rm(imp_res)
  
  return(list(imp_X = imp_X, local_sim = local_sim, droprate = droprate))
}

