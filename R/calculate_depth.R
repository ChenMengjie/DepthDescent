calculate_depth <- function(data, estimator, directions){
  XX_T <- apply(Z, 1, function(x){
    return(x%*%t(x))
  })
  project_to <- projection_directions(cbind(XX_T, as.vector(estimator)), directions)
  project.M.rank <- apply(data.frame(project_to), 1, rank)
  M.rank <- project.M.rank[nrow(project.M.rank), ]
  M.depth <- ifelse(M.rank > (nrow(data) + 1)/2, nrow(data) + 1 - M.rank, M.rank)
  depth <- min(M.depth)
  return(depth)
}