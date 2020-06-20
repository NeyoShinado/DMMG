clustering_update <- function(X, K, npc, lambda1, lambda2, W=NULL, V=NULL, Wv=NULL,
                              H=NULL, Dc=NULL, Dg=NULL, S=NULL, iteration=50, meandrop=NULL,
                              thre=1e-4, out_file="result/clust_res/res.rds", 
                              res_save=TRUE){
  cat("## clustering vars updating...\n")
  
  # initialization
  M = length(X)
  N = dim(W)[1]
  J_set = NULL
  J_DR = matrix(0, iteration, M)
  J_MG = matrix(0, iteration, M)
  J_HE = NULL
  J_LE = NULL
  par1 = NULL
  par2 = NULL
  par_sigma = NULL
  par_beta = NULL
  nmi = c()
  #nmi = matrix(0, iteration, 2)
  # eigen val
  Xe = c(sqrt(eigen(X[[1]] %*% t(X[[1]]))$values)[1], sqrt(eigen(X[[2]] %*% t(X[[2]]))$values)[1], 
         sqrt(eigen(X[[3]] %*% t(X[[3]]))$values)[1])
  #Xn1 = sapply(1:3, function(i){
  #  res = norm(X[[i]], '1')
  #  return(res)
  #})
  
  Se = eigen(S)$values[1]
  # var init
  #* 6\
  #E = H
  #theta = matrix(0, N, K)
  #sigma = 4
  #beta = 1    #* par of sparse norm
  Ds = diag(colSums(S))
  L = Ds - S
  alpha = matrix(1/M, iteration, M)
  #* Dc/Dg init
  
  
  for(iter in 1:iteration){
    cat(paste0("### Updating ", iter, "th round of ", iteration, "\n"))
    cat(paste0(strrep("#", round(30*iter/iteration)), " ", 
               100*signif(iter/iteration, digits=2), "%\n"))
    cat("### updating Wv...\n")
    SH = H %*% t(H)
    DH = diag(colSums(SH))
    #pr = alpha[iter, 1] * V[[1]] %*% t(V[[1]]) + alpha[iter, 3] * V[[3]] %*% t(V[[3]])
    #sec_ord_g = lambda1 * LH %*% W + W %*% pr + alpha[iter, 2] * Dc^2 %*% W %*% V[[2]] %*% Dg^2 %*% t(V[[2]])
    #nW = alpha[iter, 1] * X[[1]] %*% t(V[[1]]) + alpha[iter, 2] * Dc^2 %*% X[[2]] %*% Dg^2 %*% t(V[[2]]) + 
    #  alpha[iter, 3] * X[[3]] %*% t(V[[3]])
    pr = lapply(1:M, function(i){
      res = Dc^2 %*% Wv[[i]] %*% V[[i]] %*% Dg[[i]]^2 %*% t(V[[i]])
      return(res)
    })

    nW = lapply(1:M, function(i){
      res = Dc^2 %*% X[[i]] %*% Dg[[i]]^2 %*% t(V[[i]])
      return(res)
    })
    
    Wv = lapply(1:M, function(i){
      res = Wv[[i]] * ((nW[[i]] + alpha[iter, i] * W) / (pr[[i]] + alpha[iter, i] * Wv[[i]]))
      return(res)
    })
    #W = W * ((nW + sigma*E + theta) / (sec_ord_g + sigma*W))
    rm(pr, nW)
    
    
    # 5\ unitize
    #norms = rowSums(W)
    #norms[norms == 0] = 1e-10
    #W = W / matrix(norms, N, dim(W)[2])
    
    
    cat("### updating V...\n")
    V = lapply(1:M, function(i){
      res = V[[i]] * ((t(W) %*% Dc^2 %*% X[[i]] %*% Dg[[i]]^2) / (t(W) %*% Dc^2 %*% W %*% V[[i]] %*% Dg[[i]]^2+1e-10))
    })
    
    #* 5\ can not normalize W & V    #####
    #norms = NormalizeUV(Wv, V)
    #W = norms$U
    #V = norms$V
    #rm(norms)
    
    cat("### updating W...\n")
    mvpr = 0
    for(i in 1:M){
      mvpr = mvpr + alpha[iter, i] * Wv[[i]]
    }
    mvnw = sum(alpha[iter, ]) * W
    W = W * ((SH %*% W + mvpr) / (mvnw + SH %*% W))
    rm(mvpr, mvnw)
    
    
    cat("### updating H...\n")
    d = rowSums(W^2)
    d = matrix(d, nrow=N, ncol=N)
    d = d + t(d) - 2*W %*% t(W)
    
    #* 5\ ortho
    #H = H * ((lambda2 * S %*% H + 2 * H) / 
    #           ((lambda2 * Ds + lambda1 / 2 * d) %*% H + 2 *H %*% t(H) %*% H + 1/sigma * H %*% (theta + t(theta))))
    #H = eigen(0.5*lambda1 * d + lambda2 * L)$vectors[, (N-K+1):N]  #* eigen
    H = H * ((lambda2 * S %*% H) / ((lambda2 * Ds + lambda1 / 2 * d) %*% H))  #* normal
    
    #* 6\ sparse
    #H = H * ((lambda2 * S %*% H + sigma * E + theta) / ((lambda2 * Ds + lambda1 / 2 * d + sigma * diag(1, N)) %*% H))  #* normal
    norms = rowSums(H)    #* 5\row based H normalize
    norms[norms == 0] = 1e-10
    H = H / matrix(norms, N, K)
    
    #* 5\ Orthogonal constraint
    #sigma = 1.618*sigma
    #theta = theta + sigma * (t(H) %*% H - diag(1, K))
    
    ##*6\ sparse regularize
    #cat("### updating E and theta...\n")
    #E = soft(H - theta/sigma, beta/sigma)
    #sigma = 1.1 * sigma
    #theta = theta + sigma * (E - H)
    
    
    cat('### updating alpha...\n')
    J_MG[iter, ] = t(sapply(1:M, function(i){
      cost = norm(W - Wv[[i]], 'F')^2
      return(cost)
    }))
    
    # veiw check
    alpha[iter, ] = t(sapply(1:M, function(i){
      #if(i %in% ignore | i == 3){
      #  return(0)
      #}else{
        res = 0.5 / sqrt(J_MG[iter, i])
        res = res / meandrop[i]    #* 9\ weight revise by meandrop
        return(res)
      #}
    }))
    # 4\unitize weight
    #alpha[iter, ] = alpha[iter, ] / sum(alpha[iter, ])
    
    
    #** 7\update parameters    option one
    if(iter>1){# %% 10 == 0){
      #He = sqrt(eigen(H %*% t(H))$values)[1]
      #Ve = c(sqrt((eigen(V[[1]] %*% t(V[[1]]))$values)[1]), sqrt((eigen(V[[2]] %*% t(V[[2]]))$values)[1]),
      #       sqrt((eigen(V[[3]] %*% t(V[[3]]))$values)[1]))
      #lambda1 = sum(alpha[iter, ] * Xe) / sum(He + alpha[iter, ] * Ve)
      
      #de = eigen(d)$values[1]
      #lambda2 = (lambda1 * de) / (2 * Se)
      gamma = 2
      #lambda1 = (1 / J_HE[iter-1])^(1/(gamma - 1)) / sum((1 / J_HE[iter-1]) + (1 / J_LE[iter-1]))^(1/(gamma - 1))
      #lambda2 = (1 / J_LE[iter-1])^(1/(gamma - 1)) / sum((1 / J_HE[iter-1]) + (1 / J_LE[iter-1]))^(1/(gamma - 1))
      lambda1 = 2
      lambda2 = 2
      par1[iter] = lambda1
      par2[iter] = lambda2
      
      #* sparse norm
      #beta = sum(alpha[iter, ] * Xn1) / sum(alpha[iter, ] * Ve)
      #par_beta[iter] = beta
      #par_sigma[iter] = sigma
      
      #par1[(iter-9):iter] = lambda1
      #par2[(iter-9):iter] = lambda2
    }
    
    
    # cost calculation
    # if verbose
    J_DR[iter, ] = t(sapply(1:M, function(i){
      cost = norm(Dc^2 %*% (X[[i]] - Wv[[i]] %*% V[[i]]) %*% Dg[[i]]^2, 'F')^2
      return(cost)
    }))
    J_LE[iter] = sum(diag(t(H) %*% L %*% H))
    J_HE[iter] = sum(diag(t(W) %*% (DH-SH) %*% W))
    J_set[iter] = sum(J_DR[iter, ]) + sum(alpha[iter, ] * J_MG[iter, ]) + lambda1 * J_HE[iter] + lambda2 * J_LE[iter]
    cat("### Current cost:", J_set[iter], "\n")
    
    
    if(iter > 1){
      var_J = abs(J_set[iter] - J_set[iter-1])
      
      # convergence check
      if(var_J <= thre){
        cluster = sapply(1:N, function(i){
          which.max(H[i, ])
        })
        
        res = list(W = W, V = V, H = H, cluster = cluster, Dc = diag(Dc), Dg = Dg,
                   J = J_set, J_DR = J_DR, J_MG = J_MG, J_HE = J_HE, J_LE = J_LE, nmi = nmi,
                   lambda1 = par1, lambda2 = par2, alpha = alpha, dw = d)   #* lambda save
        if(res_save){
          saveRDS(res, out_file)
        }
        return(res)
      }
    }
    
    
    # recording nmi
    if(FALSE){
      cluster = specc(as.matrix(d), K)@.Data
      nmiW_sp = NMI(cluster, gt)
      
      F = eigen(d)$vectors
      F = F[, (N-K):(N-1)]
      clust = kmeans(F, K)$cluster
      nmiW_L = NMI(clust, gt)
      
      cluster = sapply(1:N, function(i){
        which.max(H[i, ])
      })
      nmiH = NMI(cluster, gt)
    }
    

    cluster = sapply(1:N, function(i){
      which.max(H[i, ])
    })
    nmiH_max = NMI(cluster, gt)
    nmi = c(nmi, nmiH_max)
    tryCatch(
      {
        cluster_k = kmeans(H, K, nstart=20, iter.max = 200)$cluster
        nmiH_kmean = NMI(cluster_k, gt)
        cat(paste0("## NMI:   nmiH_max: ", nmiH_max, "    nmiH_kmean: ", nmiH_kmean, "\n"))
        #nmi[iter, ] = c(nmiH_max, nmiH_kmean)
      }, error=function(e){
        cat(paste0("## NMI:   nmiH_max: ", nmiH_max, "\n"))
      }
    )
  }
  
  res = list(W = W, V = V, H = H, cluster = cluster, Dc = diag(Dc), Dg = Dg,
             J = J_set, J_DR = J_DR, J_MG = J_MG, J_HE = J_HE, J_LE = J_LE, nmi = nmi,
             lambda1 = par1, lambda2 = par2, beta = par_beta, sigma = par_sigma, alpha = alpha, dw = d)   #* lambda save
  if(res_save){
    saveRDS(res, out_file)
  }
  return(res)
}


soft <- function(x, theta){
  if(max(theta == 0)){
    res = x
    return(x)
  }else{
    res = abs(x) - theta
    res[which(as.matrix(res<0))] = 0
    res = sign(x) * res
    return(res)
  }
}