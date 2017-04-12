matrix_depth_by_descent <- function(Z, n, p, K = 2000, ntry = 2000, structure.directions, step.check = 20, print = FALSE){
    
  XX_T <- apply(Z, 1, function(x){
    return(x%*%t(x))
  })
  
  ######### generate directions ##########
  require(MASS)
  
  uniform.directions <- mvrnorm(K, rep(0, p), diag(p))
  uniform.directions <- apply(uniform.directions, 1, function(x){
    x/sqrt(sum(x^2))
  })
  
  if(!is.null(structure.directions)){
    all.directions <- data.frame(uniform.directions, structure.directions)  
  } else {
    all.directions <- uniform.directions
  }
  
  kendall.cor <- cor(Z, method = "kendall")
  cor.X <- sin(kendall.cor * pi /2) # this is the correlation matrix
  scale.matrix <- diag(apply(Z, 2, function(x){
    sqrt(median(x^2))
  }))
  scaled.kendall <- scale.matrix %*% cor.X %*% scale.matrix
    
  project_to <- projection_directions(XX_T, as.matrix(all.directions))
  M.initial <- as.vector(scaled.kendall)
  
  deepest <- descent_algorithm_matrix_depth(XX_T, as.matrix(M.initial), as.matrix(all.directions), ntry, step.check, print = print) 

  res <- list(deepest$best_depth, matrix(deepest$best, p), scaled.kendall)
  return(res)
}