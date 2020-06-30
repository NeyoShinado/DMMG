clustering_update <- function(X, K, NN, lambda=1, Sv=NULL, S=NULL, H=NULL, Dc=NULL, Dg=NULL, 
                              iteration=50, meandrop=NULL, thre=1e-4, out_file="result/clust_res/res.rds", 
                              res_save=TRUE){
###* be careful of the diagonal of both S,Sv and dH of multiply opt
  
  cat("## clustering vars updating...\n")
  # initialization
  M = length(X)
  N = dim(X[[1]])[1]
  NN = min(N-2, NN)
  J = NULL
  J_DE = matrix(0, iteration, M)
  J_FC = matrix(0, iteration, M)
  J_GE = matrix(0, iteration, M)
  J_GC = NULL
  alpha = matrix(0, iteration, M)
  wv = matrix(1/M, iteration, M)
  par = c()
  ari = matrix(1/M, iteration, M)
  #nmi = c()
  nmi = matrix(0, iteration, 2)
  
  #* with Drop weight
  dx = lapply(X, function(i){
    dv = rowSums(i^2)
    dv = matrix(dv, nrow=N, ncol=N)
    dv = dv + t(dv) - 2 * i %*% t(i)
    dv = dv - diag(diag(dv))
    return(dv)
  })

  #* initialize alpha & Sv
  res = lapply(1:M, function(i){
    res = update_alpha(dx[[i]], wv[1, i], S, NN, init=TRUE)
    return(res)
  })
  #####
  alphas = sapply(res, function(i){
    alpha = i$alpha
    return(alpha)
  })
  for(i in 1:M){
    alpha[, i] = alphas[i]
  }
  
  Sv = lapply(res, function(i){
    Sv = i$Sv
    return(Sv)
  })
  rm(res)
  
## optimize version 1 of matrix multiplication
  for(iter in 1:iteration){

    cat(paste0("### Updating ", iter, "th round of ", iteration, "\n"))
    cat(paste0(strrep("#", round(30*iter/iteration)), " ", 
               100*signif(iter/iteration, digits=2), "%\n"))
    
    cat("### updating Sv...\n")
    Sv = lapply(1:M, function(i){
      res = Sv[[i]] * (2*wv[iter, i]*S) / (2*wv[iter, i]*Sv[[i]] + 2*alpha[iter, i]*Sv[[i]]+dx[[i]])
      return(res)
    })
    
    #* 1\norm
    Sv = lapply(Sv, function(i){
      N = dim(i)[1]
      norm =  rowSums(i)
      norm = matrix(norm, N, N)
      i = i / (norm + 1e-10)
      i = (i + t(i)) / 2
      #i = i - diag(diag(i)) + diag(1, N)
      return(i)
    })
    

    cat("### updating S...\n")
    dH = rowSums(H^2)
    dH = matrix(dH, nrow=N, ncol=N)
    dH = dH + t(dH) - 2*H %*% t(H)
    dH = dH - diag(diag(dH))
    
    temp = lapply(1:M, function(i){
      res = wv[iter, i] * Sv[[i]]
      return(res)
    })
    S = S * (2*Reduce(f="+", temp)) / (2*sum(wv[iter, ]) * S + lambda * dH)
    
    #* 1\norm
    norm = rowSums(S)
    norm = matrix(norm, N, N)
    S = S / (norm + 1e-10)
    S = (S + t(S)) / 2
    #S = S - diag(diag(S)) + diag(1, N)
    rm(dH, temp)

    
    cat("### updating wv...\n")
    J_GE[iter, ] = sapply(1:M, function(i){
      res = norm(S - Sv[[i]], 'F')^2
      return(res)
    })
    
    wv[iter, ] = sapply(1:M, function(i){
      res = 0.5 / sqrt(J_GE[iter, i])
      return(res)
    })
    
    browser()
    cat("### updating alpha...\n")#####
    #alpha[iter, ] = sapply(1:M, function(i){
    #  res = update_alpha(dx[[i]], wv[iter, i], S, NN)
    #  return(res)
    #})
    
    
    cat('### updating H...\n')
    D = diag(colSums(S))
    H = H * (S %*% H) / (D %*% H)
    
    #* 5\ H unitize
    norms = rowSums(H)    #* 5\row based H normalize
    norms[norms == 0] = 1e-10
    H = H / matrix(norms, N, K)
    
    
    #** 7\update parameters
    if(iter>1){
      par[iter] = lambda
    }
    
    
    # cost calculation
    temp_DE = lapply(1:M, function(i){
      res = dx[[i]] * Sv[[i]]
      return(res)
    })
    J_DE[iter, ] = sapply(temp_DE, sum)
    rm(temp_DE)
    J_FC[iter, ] = sapply(1:M, function(i){
      res = alpha[iter, i] * norm(Sv[[i]], "F")^2
      return(res)
    })
    J_GE[iter, ] = sapply(1:M, function(i){
      res = wv[iter, i] * norm(S - Sv[[i]], "F")^2
      return(res)
    })
    J_GC[iter] = sum(diag(t(H) %*% S %*% H))
    J[iter] = sum(J_DE[iter, ]) + sum(J_FC[iter, ]) + sum(J_GE[iter, ]) + lambda * J_GC[iter]
    cat("### Current cost:", J[iter], "\n")
    
    
    if(FALSE){#iter > 1){
      var_J = abs(J[iter] - J[iter-1])
      
      # convergence check
      if(var_J <= thre){
        cluster = sapply(1:N, function(i){
          which.max(H[i, ])
        })
        
        res = list(S = S, Sv = Sv, H = H, cluster = cluster, Dc = Dc, Dg = Dg,
                   J = J, J_DE = J_DE, J_FC = J_FC, J_GE = J_GE, J_GC = J_GC, nmi = nmi, ari = ari, 
                   lambda = par, alpha = alpha, wv = wv)
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
    ari[iter, 1] = ARI(cluster, gt)
    nmi[iter, 1] = nmiH_max
    #nmi = c(nmi, nmiH_max)
    tryCatch(
      {
        cluster_k = kmeans(H, K, nstart=20, iter.max = 200)$cluster
        nmiH_kmean = NMI(cluster_k, gt)
        cat(paste0("## NMI:   nmiH_max: ", nmiH_max, "    nmiH_kmean: ", nmiH_kmean, "\n"))
        nmi[iter, 1:2] = c(nmiH_max, nmiH_kmean)
        ari[iter, 2] = ARI(cluster_k, gt)
      }, error=function(e){
        cat(paste0("## NMI:   nmiH_max: ", nmiH_max, "\n"))
      }
    )
  }
  
  res = list(S = S, Sv = Sv, H = H, cluster = cluster, Dc = Dc, Dg = Dg,
             J = J, J_DE = J_DE, J_FC = J_FC, J_GE = J_GE, J_GC = J_GC, nmi = nmi, ari = ari, 
             lambda = par, alpha = alpha, wv = wv)
  if(res_save){
    saveRDS(res, out_file)
  }

  return(res)
}


# sub_function
update_alpha <- function(dx, w, S, NN, init=FALSE){
  N = nrow(S)

  dist = dx - 2 * w * S
  id = t(apply(dist, 1, function(i){
    res = order(i, decreasing=FALSE)
    return(res)
  }))
  dist = t(apply(dist, 1, function(i){
    res = sort(i, decreasing=FALSE)
    return(res)
  }))
  
  
  alphas = sapply(1:N, function(i){
    di = dist[i,2:(NN+2)]
    res = 0.5*(NN*di[NN+1]-sum(di[1:NN])) - w
    return(res)
  })

  alpha = mean(alphas)
  
  if(init){
    # initialize S
    Sv = t(sapply(1:N, function(i){
      res = matrix(1e-10, 1, N)
      di = dist[i,2:(NN+2)]
      res[id[i, 2:(NN+2)]] = (di[NN+1]-di)/(NN*di[NN+1] - sum(di[1:NN]) + 1e-10)
      return(res)
    }))
    
    Sv = (Sv + t(Sv)) / 2
    return(list(alpha=alpha ,Sv=Sv))
  }else{
    return(alpha)
  }
}