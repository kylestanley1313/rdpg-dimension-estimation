## =============================================================================
## PROJECT: STAT 597 (NET) -- Final Project
## PURPOSE: Simulation Code
## DATE: Fall 2022
## AUTHOR: Kyle Stanley
## =============================================================================

library(pbmcapply)


## Utilities ===================================================================


gen_config <- function(id, reps, N, D, B) {
  list(
    id = id, 
    reps = reps, 
    N = N, 
    D = D, 
    B = B, 
    support = list(min = rep(0, D), max = rep(1/sqrt(D), D))
  )
}


sim_A_from_Z <- function(Z) {
  
  N <- nrow(Z)
  
  ## simulate adjacency matrix
  A <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      if (i < j) {
        p <- Z[i,] %*% Z[j,]
        p <- min(1, max(0, p))
        A[i,j] <- rbinom(n = 1, size = 1, prob = p)
      }
    }
  }
  A <- A + t(A)  ## symmetricize
  
  return(A)
  
}


sim_A_from_config <- function(config) {
  
  ## unpack config
  N <- config$N
  D <- config$D
  support <- config$support
  
  ## simulate latent positions
  Z <- matrix(nrow = N, ncol = D)
  for (n in 1:N) {
    for (d in 1:D) {
      min <- support$min[d]
      max <- support$max[d]
      Z[n,d] <- runif(1, min = min, max = max)
    }
  }
  
  ## simulate adjacency matrix
  A <- sim_A_from_Z(Z)
  
  return(list(
    Z = Z,
    A = A
  ))
  
}


estimate_Z <- function(A.vecs, A.vals, dim) {
  if (dim == 1) {
    Z.hat <- (matrix(A.vecs[,dim], ncol = 1)) * sqrt(A.vals[dim])
  }
  else {
    Z.hat <- A.vecs[,1:dim] %*% diag(sqrt(A.vals[1:dim]))
  }
  return(Z.hat)
}


test_dim <- function(B, A.vecs, A.vals, D.null) {
  
  lambdas.star <- rep(NA, B)
  Z.hat <- estimate_Z(A.vecs, A.vals, D.null)
  
  for (b in 1:B) {
    A.star <- sim_A_from_Z(Z.hat)
    lambdas.star[b] <- svd(A.star)$d[D.null+1]
  }
  
  p <- sum(lambdas.star > A.vals[D.null+1]) / B
  return(p)
  
}


test_dims <- function(B, A, D) {
  
  A.svd <- svd(A)
  cols <- c('D.null', 'p')
  res <- data.frame(matrix(nrow = 0, ncol = length(cols)))
  colnames(res) <- cols
  
  for (d in 1:(D+2)) {
    p <- test_dim(B, A.svd$u, A.svd$d, d)
    res[nrow(res)+1,] <- c(d, p)
  }
  
  return(res)
  
}

process_config <- function(config) {
  A <- sim_A_from_config(config)$A
  
  cols <- c('id', 'rep', 'N', 'D', 'D.null', 'p')
  res <- data.frame(matrix(nrow = 0, ncol = length(cols)))
  
  for (r in 1:config$reps) {
    temp <- test_dims(config$B, A, config$D)
    temp$id <- config$id
    temp$rep <- r
    temp$N <- config$N
    temp$D <- config$D
    res <- rbind(res, temp)
  }
  return(res)
}


## Execution ===================================================================


## Generate configs
Ns <- c(100, 250, 500, 1000)
Ds <- c(2, 3, 4, 5)
reps <- 20
B <- 20
configs <- list()
id <- 1
for (N in Ns) {
  for (D in Ds) {
    configs[[id]] <- gen_config(id, reps, N, D, B) 
    id <- id + 1
  }
}

configs <- list(configs[[1]], configs[[2]])  ## DEBUGGING


## Process configs in parallel
print("----- START PROCESSING -----")
set.seed(1)
out <- pbmclapply(configs, process_config, ignore.interactive = TRUE)
print("----- END PROCESSING -----")

## Process results
res <- out[[1]]
for (i in 2:length(out)) {
  res <- rbind(res, out[[i]])
}
res <- res[,c('id', 'rep', 'N', 'D', 'D.null', 'p')]

## Write results
write.table(
  res, 
  file.path('results', 'results01.csv'), 
  sep = ',',
  col.names = TRUE,
  row.names = FALSE,
  append = FALSE,
  quote = FALSE
)




