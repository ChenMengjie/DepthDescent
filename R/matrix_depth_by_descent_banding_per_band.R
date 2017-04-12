matrix_depth_by_descent_banding_per_band <- function(Z, n, p, per_band = 1000, ntry = 2000, structure.directions, step.check = 20, band){
  
  XX_T <- apply(Z, 1, function(x){
    return(x%*%t(x))
  })
  
  ######### generate directions ##########
  require(MASS)
  
  directions <- matrix(0, ncol = per_band*(p - 2*band+1), nrow = p)
  
  uniform.directions <- mvrnorm(per_band*(p - 2*band+1), rep(0, 2*band), diag(2*band))
  
  uniform.directions <- apply(uniform.directions, 1, function(x){
    x/sqrt(sum(x^2))
  })
  
  for(i in 1:(p - 2*band+1)){
    directions[i:(i+2*band-1), c(1:per_band) + per_band*(i-1)] <- uniform.directions[, c(1:per_band) + per_band*(i-1)]
  }  
  
  if(!is.null(structure.directions)){
    all.directions <- data.frame(directions, structure.directions)  
  } else {
    all.directions <- directions
  }
  
  kendall.cor <- cor(Z, method = "kendall")
  cor.X <- sin(kendall.cor * pi /2) # this is the correlation matrix
  scale.matrix <- diag(apply(Z, 2, function(x){
    sqrt(median(x^2))
  }))
  
  scaled.kendall <- scale.matrix %*% cor.X %*% scale.matrix
  
  banded.kendall <- banding_a_matrix(scaled.kendall, band)
  
  project_to <- projection_directions(XX_T, as.matrix(all.directions))
  M.initial <- as.vector(banded.kendall)
  
  deepest <- descent_algorithm_matrix_depth_bandable(XX_T, as.matrix(M.initial), as.matrix(all.directions), ntry, step.check, band) 
  
  res <- list(deepest$best_depth, matrix(deepest$best, p), banded.kendall)
  return(res)
}