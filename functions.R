library(hitandrun)
library(igraph)
library(Matrix)

##### Sampling functions #####

sample_SAD <- function(n, # number of species
                       constant = FALSE, # return equal abundances (1/n) for all species?
                       lambda = 1 # rate parameter for exponential distribution
                       ){
  # sample a random vector of abundances for n species
  # current implementation uses iid exponential samples with rate lambda, but any distribution could be substituted
  # normalize abundances to sum to one wlog
  
  x <- if(constant) rep(1, n) else rexp(n, lambda)
  x <- x / sum(x)
  return(x)
}


sample_A <- function(n, # number of species
                     aij_dist = rexp, # distribution for matrix elements (exp with lambda = 1 by default)
                     norm_rows = FALSE # should rows be normalized to sum to one?
                     ){
  # sample a random matrix with iid elements (followed by optional row normalization)
  
  A <- matrix(aij_dist(n^2), n, n)
  if(norm_rows) A <- diag(1 / rowSums(A)) %*% A
  return(A)
}


sample_r <- function(n, # number of species
                     r_dist = rexp, # distribution for individual growth rates (exp with lambda = 1 by default)
                     constant = FALSE, # return equal growth rates for all species?
                     normalize = FALSE # should growth rates be normalized to sum to one?
                     ){
  # sample a random vector of growth rates 
  # iid from input distribution (exp by default)
  
  r <- if(constant) rep(1, n) else r_dist(n)
  if(normalize) r <- r / sum(r)
  return(r)
}


sample_A_r_pair <- function(x, # vector of species abundances
                            require_stable = TRUE, # keep sampling until a stable equilibrium is found?
                            method, # sampling method (one of "sample_A_solve_r", "sample_r_renormalize_A", "hit_and_run", "assemble_from_pool")
                            params = c(1,1) # method-specific parameters
                            ){
  # sample an interaction matrix (A) and growth rate vector (r) pair for an input SAD (x) such that A x = r
  # optional rejection sampling is used to find A, r pairs corresponding to stable equilibrium
  # sampling methods include:
  ##  "sample_A_solve_r": sample A iid and then solve for r = A x
  ##  "sample_r_renormalize_A": sample A' and r iid, then define A = D(r) D(1/rowSums(A)) A D(1/x)
  ##  "hit_and_run": use hit-and-run sampling algorithm to sample uniformly given constraints A x = r and 0< r < r_max and/or 0 < a_ij < a_max
  ##  "assemble_from_pool": sample A' and r' iid, then integrate GLV dynamics to find persistent set of species. Return A and r for the persistent set (n_final <= n_initial)
  
  n <- length(x)
  stability_flag <- require_stable
  
  repeat{
    
    if(method == "sample_A_solve_r"){
      
      A <- sample_A(n)
      r <- A %*% x
    }
    
    if(method == "sample_r_renormalize_A"){
      
      r <- sample_r(n, constant = TRUE)
      A <- diag(as.vector(r)) %*% sample_A(n, norm_rows = TRUE) %*% diag(1 / x)
    }
    
    if(method == "hit_and_run"){ # NOTE: this method slows down rapidly with increasing n
      
      r_max <- params[1]
      a_max <- params[2]
      
      # set up system of constraints for hitandrun sampler
      
      # matrix of constraints
      constr <- rbind(
        cbind(kronecker(diag(n), t(x)), -diag(n)), # constraints of the form A_(i) x - r_i = 0
        diag(n * (n + 1)), # constraints of the form a_ij < a_max and/or r_i < r_max
        -diag(n * (n + 1)) # constraints of the form a_ij > 0 and r_i > 0
      )
      
      # type of constraint (equality or inequality)
      dir <- c(rep("=", n), rep("<=", n * (n + 1)), rep("<=", n * (n + 1)))
      
      # right-hand side of constraint
      rhs <- c(rep(0, n), rep(a_max, n^2), rep(r_max, n), rep(0, n * (n + 1)))
      
      out <- hitandrun(constr = list(constr = constr, dir = dir, rhs = rhs), 
                       n.samples = 1)
      
      A <- matrix(out[, 1:(n^2)], n, n, byrow = TRUE)
      r <- out[, -(1:(n^2))]
    }
    
    if(method == "assemble_from_pool"){
      
      # still in development
    }
    
    # if require_stable, check for stability
    # break if stable or resample and try again
    if(!require_stable | all(Re(eigen(-diag(x) %*% A)$val) < 0)) break
  }
  
  r <- as.vector(r) # ensure r is a vector (not matrix)
  
  return(list(A = A, r = r))
}

##### Analysis functions #####

build_pairwise_outcome_matrix <- function(A, r){
  
  # compute an n x n matrix of pairwise competitive outcomes
  # Codes for outcomes:
  ##  M_ij = 1 : i excludes j
  ##  M_ij = -1 : j excludes i
  ##  M_ij = 0 : coexistence
  ##  M_ij = NA : priority effect (no coexistence, outcomes depends on initial conditions)
  # Note diagonals are set to NA
  
  M <- diag(1 / r) %*% A # create matrix of a_ij / r_i (which determines pairwise outcomes)
  
  Q <- M * 0 
  for(i in 1:nrow(M)){
    for(j in 1:ncol(M)){
      if((M[i, i] > M[j, i]) & (M[i, j] > M[j, j])) Q[i,j] <-  -1 # j beats i
      if((M[i, i] < M[j, i]) & (M[i, j] < M[j, j])) Q[i,j] <-   1 # i beats j
      if((M[i, i] > M[j, i]) & (M[i, j] < M[j, j])) Q[i,j] <-   0 # coexistence
      if((M[i, i] < M[j, i]) & (M[i, j] > M[j, j])) Q[i,j] <-   NA # priority effect
    }
  }
  diag(Q) <- NA
  return(Q)
}



compute_hierarchy_score <- function(Q){
  
  # compute fraction of pairwise exclusions that go against hierarchy, as determined following Chang et al. 2023
  # output is random due to random tie breaking
  
  Q <- replace(Q, is.na(Q), 0) # replace NA by 0
  n <- nrow(Q)
  
  outcome_diff <- rowSums(Q) # for each species, number of pairwise competitions won - lost
  comp_rank <- rank(outcome_diff, ties = "random") # competitive rank based on pairwise outcome differential
  
  Q <- Q[order(comp_rank), order(comp_rank)] # reorder Q by competitive ranking 
  against_rank <- sum(Q[lower.tri(Q)] == -1) # count exclusions that go "against rank"
  all_exclusions <- sum(Q[lower.tri(Q)] != 0) # total number of pairwise exclusions
  
  return(against_rank / all_exclusions) # return the ratio
}

count_intransitive_triangles <- function(Q){
  
  # count the number of directed triangles (RPS outcomes) in the network of pairwise outcomes
  
  Q <- replace(Q, (is.na(Q) | Q == -1), 0) # replace NA and -1 by 0 (remaining 1s are directed edges)
  g <- graph_from_adjacency_matrix(Q) # convert matrix to directed graph
  t <- triad_census(g)[10] # count directed triangles
  
  return(t)
}

count_fraction_excluded <- function(Q){
  
  # compute the fraction of pairwise outcomes that are competitive exclusion
  
  Q <- replace(Q, (is.na(Q)), 0) # replace NA by 0
  n <- nrow(Q)
  e <- sum(Q == -1)
  
  return(e / (n * (n - 1) / 2)) # return fraction of all pairwise outcomes that are exclusions
}

count_fraction_coexist <- function(Q){
  
  # compute the fraction of pairwise outcomes that are coexistence
  
  Q <- replace(Q, (is.na(Q)), Inf) # replace NA by Inf (since we will count 0s)
  n <- nrow(Q) 
  e <- sum(Q == 0)
  
  return(e / (n * (n - 1))) # return fraction of all pairwise outcomes that are coexistence (NOTE: 0s are counted twice for each pairwise coexistence)
}

test_for_loops <- function(Q){
  
  # check if there are any directed (intransitive) loops of pairwise exclusions 
  
  Q <- replace(Q, (is.na(Q) | Q == -1), 0) # replace NA and -1 by 0
  expQ <- expm(Q > 0) # compute matrix exponential -- loops will cause diagonal > 1
  return(all(diag(expQ) != 1)) # if all diagonal == 1, there are no directed loops
}

count_fraction_positive_net_effects <- function(A){
  
  # count the fractin of net effects (elements of A^-1) that are positive
  
  n <- nrow(A)
  invA <- solve(A)
  
  return(sum(invA < 0) / (n * (n - 1)))
}