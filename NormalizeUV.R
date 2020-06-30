# normalize updating vars

NormalizeUV <- function(U, V, Norm=2){
  # Object term ||X - UV||,  which X ~ N*P, noramlize U by default
  # Norm -- 2: Fro Norm
  #      -- 1: 1 Norm
  if(TRUE){
    M = length(U)
    N = dim(U[[1]])[1]
    k = dim(U[[1]])[2]
    
    if(Norm != 1 && Norm != 2){
      stop("## Incorrect nrom pars for norm type identification!")
    }
    
    
    if(Norm == 2){
      U = lapply(1:M, function(i){
        norms = sqrt(colSums(U[[i]]^2))
        norms[norms == 0] = 1e-10
        res = U[[i]] / t(matrix(norms, k, N))
        V[[i]] = V[[i]] * matrix(norms, k, dim(V[[i]])[2])
        return(res)
      })
    }else{
      U = lapply(1:M, function(i){
        norms = sqrt(colSums(abs(U[[i]])))
        norms[norms == 0] = 1e-10
        res = U[[i]] / t(matrix(norms, k, N))
        V[[i]] = V[[i]] * matrix(norms, k, dim(V[[i]])[2])
        return(res)
      })
    }
    
    
    res = list(U = U, V = V)
    return(res)
  }else{
    
    ## row normalize
    M = length(U)
    N = dim(U[[1]])[1]
    k = dim(U[[1]])[2]
    U = lapply(1:M, function(i){
      norms = rowSums(U[[i]])
      norms[norms == 0] = 1e-10
      res = U[[i]] / matrix(norms, N, k)
      return(res)
    })
  }
  res = list(U = U, V = V)   
}