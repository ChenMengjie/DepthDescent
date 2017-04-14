###  DepthDescent


### Introduction

**DepthDescent** implements an algorithm to compute the deepest matrix estimator proposed in https://arxiv.org/abs/1506.00691.  The core idea of the algorithm is to start with some robust scatter matrix estimator that is easy to compute, and then iteratively evaluate its depth on some direction subset of finite cardinality to approximate its matrix depth function and move towards the directions that can improve the depth. 

### Installation

**DepthDescent** relies on the following R packages: **Rcpp** and **RcppArmadillo**. These package need to be pre-installed. 

  ```R
  install.packages ("Rcpp")
  install.packages ("RcppArmadillo")
  ```
  
**DepthDescent** can be installed from github directly as follows:

  ```R
  install.packages ("devtools")
  library(devtools)
  install_github("ChenMengjie/DepthDescent")
  ```
  
### Tutorial

  ```R
p <- 10
Sigma1 <- matrix(0, nrow = p, ncol = p)
for(i in 1:p){
	for(j in 1:p){
		Sigma1[i, j] <- 4*(0.5)^abs(i-j)	
	}
}
require(MASS)
per <- 0.1 # Percentage of contamination
indicators <- sample(c(0, 1), n, replace = TRUE, prob = c(per, 1 - per))
outlier <- length(indicators[indicators == 0])
Z1 <- mvrnorm(n-outlier, rep(0, p), Sigma1)
if(per!=0){
	Z2 <- matrix(10, ncol = p, nrow = outlier)
} else {
	Z <- Z1
}
depth <- matrix_depth_by_descent(Z, n, p, K = 500, ntry = 500)
  ```
 **K** is the number of directions to use. *ntry* is the number of iterations to run.
  
### Author

**Mengjie Chen** (UChicago)

Bug report, comments or questions please send to mengjiechen@uchicago.edu.
