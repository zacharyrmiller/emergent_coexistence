# Simulation functions to accompany "Multispecies coexistence emerges from pairwise competitive exclusions in communities with competitive hierarchy" by Z.R. Miller and D. Max
# Questions about code? Contact: zachary.miller@yale.edu

library(igraph)
library(Matrix)
library(deSolve)
library(tidyverse)

### Functions used for numerical simulations, analysis of simulation results, and visualizations ###

##### Sampling functions #####

sample_SAD <- function(n, # number of species
                       constant = FALSE, # return equal abundances (1/n) for all species?
                       lambda = 1 # rate parameter for exponential distribution
){
  # sample a random vector of abundances for n species
  # current implementation uses iid exponential samples with rate lambda, but any distribution could be substituted
  
  x <- if(constant) rep(1, n) else rexp(n, lambda)
  x <- x / sum(x) # normalize abundances to sum to one wlog
  
  return(x)
}


sample_A <- function(n, # number of species
                     aij_dist = rexp, # distribution for matrix elements (exp with lambda = 1 by default),
                     adjust_symmetry = TRUE, # if TRUE, each ij/ji pair is summed and redistributed according to a symmetric Beta R.V. with parameter given by beta_param 
                     beta_param = 1, # parameter of the beta distribution; values > 1 make interactions more symmetric
                     norm_rows = FALSE # should rows be normalized to sum to one?
){
  # sample a random matrix with iid elements 
  # followed by optional adjustment of symmetry (using symmetric Beta) and/or row normalization
  
  A <- matrix(aij_dist(n^2), n, n)
  
  if(adjust_symmetry){
    
    beta_mat <- matrix(rbeta(n^2, beta_param, beta_param), n, n)
    beta_mat[lower.tri(beta_mat)] <- 1 - t(beta_mat)[lower.tri(beta_mat)]
    diag(beta_mat) <- 1/2
    
    A <- (A + t(A)) * beta_mat # redistribute (A_ij + A_ji) according to iid Betas
  }
  if(norm_rows) A <- diag(1 / rowSums(A)) %*% A # normalize rows to sum to 1
  
  return(A)
}


sample_r <- function(n, # number of species
                     r_dist = rexp, # distribution for individual growth rates (exp with lambda = 1 by default)
                     normalize = FALSE # should growth rates be normalized to sum to one?
){
  # sample a random vector of growth rates 
  # iid from input distribution (exp by default)
  
  r <- r_dist(n)
  if(normalize) r <- r / sum(r)
  
  return(r)
}


sample_A_r_pair <- function(x, # vector of species abundances
                            require_stable = TRUE, # keep sampling until a stable equilibrium is found?
                            method, # sampling method (one of "sample_A_solve_r", "sample_r_renormalize_A", "assemble_from_pool")
                            params = c(1,1), # method-specific parameters,
                            thresh = 10^-6, # persistence threshold (only for "assemble_from_pool" method)
                            ...
){
  # sample an interaction matrix (A) and growth rate vector (r) pair for an input SAD (x) such that A x = r
  # optional rejection sampling is used to find A, r pairs corresponding to stable equilibrium
  # sampling methods include:
  ##  "sample_A_solve_r": sample A iid and then solve for r = A x
  ##  "sample_r_renormalize_A": sample A' and r iid, then define A = D(r) D(1/rowSums(A)) A D(1/x)
  ##  "assemble_from_pool": sample A' and r' iid, then integrate GLV dynamics to find persistent set of species. Return A and r for the persistent set (n_final <= n_initial)
  ### Note that "assemble_from_pool" is much slower and usually generates very small assembled communities
  
  n <- length(x)
  stability_flag <- require_stable
  
  repeat{ # do loop for stability rejection sampling (if require_stable = FALSE, exit after one pass)
    
    if(method == "sample_A_solve_r"){
      
      A <- sample_A(n, ...)
      r <- A %*% x # solve r s.t. A x = r
    }
    
    if(method == "sample_r_renormalize_A"){
      
      r <- sample_r(n, r_dist = function(x) rpower(x, n, 1)) # use power distribution to sample uniformly over parameter space
      A <- diag(as.vector(r)) %*% sample_A(n, norm_rows = TRUE, ...) %*% diag(1 / x) # sample random (iid) A and then normalize s.t. A x = r
    }
    
    if(method == "assemble_from_pool"){
      
      n_greater_than_two <- FALSE # community size flag to check if assembled community has more than two species
      while(!n_greater_than_two){ # keep sampling to find community with more than two species
        
        x <- rep(1/n, n) # initial condition for numerical integration (all species equally abundant)
        flag <- FALSE # flag for convergence to equilibrium
        counter <- 0 # counter for restarts of the numerical integration
        
        while(!flag & counter < 10){ # keep integrating until the dynamics converge to equilibrium OR the numerical integration has continued 10 times (if so, re-try with new parameters)
          
          A <- sample_A(n, ...)
          r <- sample_r(n)
          
          parameters <- list(A = -as.vector(A), r = r, n = n)
          times <- seq(0, 100, by = 1)
          
          tryCatch(out <- ode(y = x, time = times, func = glv, parms = parameters, method = "ode45"),
                   error = function(cond) x <- rep(0, n),
                   warning = function(cond) x <- rep(0, n))
          out[out < thresh] <- 0 # threshold abundances 
          
          x <- out[nrow(out), -1]
          
          if(all(abs(out[nrow(out), -1] - out[nrow(out) - 1, -1]) / pmax(out[nrow(out), -1], thresh) < thresh)){ # check if change in abundances between final two timesteps is below threshold
            flag <- TRUE # if change is below threshold, convergence has occurred, so exit loop
          }
          counter <- counter + 1
        }
        
        if(counter == 10) x <- rep(0, n) # if counter has reached 10, numerical integration failed to converge to equilibrium
        
        present <- x > 0 # vector of persistent species
        
        if(sum(present) > 2){ # check if more than two species persisted
          
          n_greater_than_two <- TRUE
          
          # subset to persistent species
          A <- A[present, present]
          r <- r[present]
          x <- x[present]
        }
      }
    }
    
    # if require_stable = TRUE, check for stability
    # break if stable or resample and try again
    if(!require_stable | all(Re(eigen(-diag(x) %*% A)$val) < 0)) break
  }
  
  r <- as.vector(r) # ensure r is a vector (not matrix)
  
  return(list(A = A, r = r))
}

rpower <- function(n, p, max = 1){
  # sample from the power distribution on (0, max) with power p
  # to sample parameters uniformly (subject to r_max), p = n (number of species)
  
  x <- runif(n)
  out <- max * x^(1/p) # use inverse transform sampling
  
  return(out)
}

glv <- function(t, x, parameters, thresh = 10^-8) {
  with(as.list(c(x, parameters)), {
    
    # GLV dynamics for numerical integration with deSolve
    
    A <- matrix(A, nrow = n)
    
    x[x < thresh] <- 0
    dx <- x * (r + A %*% x)
    
    return(list(dx))
  })
}


##### Data analysis functions #####

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


count_total_triangles <- function(Q){
  
  # count the number of triangles formed by exclusions in the network of pairwise outcomes
  
  Q <- replace(Q, (is.na(Q) | Q == -1), 0) # replace NA and -1 by 0 (remaining 1s are directed edges)
  g <- graph_from_adjacency_matrix(Q) # convert matrix to directed graph
  t <- length(triangles(g)) / 3 # count directed triangles
  
  return(t)
}


count_intransitive_triangles <- function(Q){
  
  # count the number of directed triangles (rock-paper-scissors motifs) in the network of pairwise outcomes
  
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
  C <- sum(Q == 0)
  
  return(C / (n * (n - 1))) # return fraction of all pairwise outcomes that are coexistence (NOTE: 0s are counted twice for each pairwise coexistence)
}


test_for_loops <- function(Q){
  
  # check if there are *any* directed (intransitive) loops of pairwise exclusions 
  
  Q <- replace(Q, (is.na(Q) | Q == -1), 0) # replace NA and -1 by 0
  expQ <- expm(Q > 0) # compute matrix exponential -- loops will cause diagonal > 1
  
  return(!all(diag(expQ) == 1)) # if all diagonal == 1, there are no directed loops
}


calculate_inter_vs_intra_ratio <- function(A){
  
  # calculate the ratio of mean inter- vs. mean intraspecific interactions 
  # (i.e. ratio of mean off-diagonal to mean diagonal elements of the interaction matrix A)
  
  offdiag <- A[(1 - diag(nrow(A))) == 1] # get just the off diagonal elements of A
  diag <- A[diag(nrow(A)) == 1] # get just the diagonal elements of A
  mean_ratio <- mean(offdiag) / mean(diag) # ratio of means
  
  return(mean_ratio)
}


calculate_positive_effects_frac <- function(A, r){
  
  # calculate the number of indirect (net) effects (i.e. elements of -A^-1) that are positive
  
  A_inv <- -solve(A)
  diag(A_inv) <- 0
  pos <- sum(A_inv > 0) / (2 * choose(nrow(A), 2))
  
  return(pos)
}


build_random_directed_network <- function(n, C){
  
  # build a competitive outcome network with choose(n,2) * C edges (competitive exclusions) placed and oriented uniformly at random
  
  # get the linear (vectorized) index for each matrix position
  ids <- matrix(1:n^2, nrow = n)
  lower_ids <- ids[lower.tri(ids)]
  upper_ids <- t(ids)[lower.tri(ids)]
  
  idx <- sample(1:choose(n,2), choose(n,2) * C) # sample choose(n,2) * C positions in the adjacency matrix to have edges
  lower <- sample(c(T, F), choose(n, 2) * C, replace = T) # for each edge, choose the orientation (i.e. does the 1 appear in the top or bottom of the matrix) randomly
  
  winners <- c(lower_ids[idx[lower]], upper_ids[idx[!lower]]) # collection of indices corresponding to each pairwise exclusion (i.e. an index is included if corresponding g_ij = 1 meaning i beats j)
  
  zeros <- rep(0, n^2)
  zeros[winners] <- 1 # add 1s for winners of each competitive exclusion
  
  g <- matrix(zeros, nrow = n) # reshape into an n x n matrix
  
  return(g)
}

##### Visualization code #####

make_r_from_x <- function(x){
  
  # construct a sequence of growth rates given equilibrium abundances (Eq. 3 in manuscript)
  
  cumulative_x <- c(0, cumsum(x)[-length(x)])
  r <- c(cumulative_x + 0.5 * x)
  
  return(r)
}


make_pairwise_outcomes_plot <- function(r, # sequence of growth rates
                                        color, # color for filling in competitive outcomes
                                        m = NA, # a local extinction value (for competition-colonization trade-off)
                                                # if NA, the model is assumed to be growth-competition trade-off
                                                # if a value is supplied, the competition-colonization trade-off is used
                                        axis.title = FALSE
                                        ){
  # plot the outcome (competitive exclusion or coexistence) for every pair of species in the growth-competition or competition-colonization trade-off models
  # generates a square for each species pair, which is filled to indicate coexistence or empty to indicate exclusion
  
  n <- length(r)
  
  if(!is.na(m)){ # if m is supplied, convert growth rates to colonization rates and use competition-colonization trade-off model
    
    c <- r + m
    A <- diag(c) %*% lower.tri(diag(n)) + lower.tri(diag(n)) %*% diag(c) # build interaction matrix according to competition-colonization model
    diag(A) <- c
    
  }else{
    
    A <- lower.tri(diag(n)) + 0 # build interaction matrix according to growth-competition model
    diag(A) <- 0.5
  }
  
  outcomes <- build_pairwise_outcome_matrix(A, r) # get all pairwise competitive outcomes
  
  long_form <- as.data.frame.table(outcomes)
  long_form[1:2] <- lapply(long_form[1:2], as.numeric) # make species indices numeric values
  p <- long_form %>% 
    filter(Var1 > Var2) %>%
    add_row(Var1 = NA, Var2 = NA, Freq = -1, .before = 1) %>% # this is a hack to make sure factors are consistent 
    ggplot() + 
    aes(xmin = Var1 - 0.4, xmax = Var1 + 0.4, ymin = Var2 - 0.4, ymax = Var2 + 0.4, 
        alpha = as.factor(Freq)) + 
    geom_rrect(color = color, fill = color) + 
    scale_alpha_manual(values = c(0, 1)) + 
    scale_x_continuous(breaks = c(2, ceiling(n/2), n), expand = c(0.02,0.02)) + 
    scale_y_continuous(breaks = c(1, ceiling(n/2)-1, n - 1), expand = c(0.02, 0.02), position = "right") +
    theme_classic() +
    theme(legend.position = "none",
          axis.line = element_blank())
  
  if(axis.title){
    p <- p + xlab("Competitive rank") + ylab("Competitive rank")
  }else{
    p <- p + theme(axis.title = element_blank())
  }
  
  return(p)
}


## functions below are for plotting rounded rectangles -- code adapted from https://rdrr.io/github/hrbrmstr/ggchicklet/man/geom_rrect.html

geom_rrect <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       radius = grid::unit(1, "pt"),
                       ...,
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRrect,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      radius = radius,
      na.rm = na.rm,
      ...
    )
  )
}

GeomRrect <- ggplot2::ggproto(
  "GeomRrect", ggplot2::Geom,
  default_aes = ggplot2::aes(
    colour = NA, fill = "grey35", size = 0.5, linetype = 1, alpha = NA
  ),
  
  required_aes = c("xmin", "xmax", "ymin", "ymax"),
  
  draw_panel = function(self, data, panel_params, coord,
                        radius = grid::unit(6, "pt")) {
    
    coords <- coord$transform(data, panel_params)
    
    lapply(1:length(coords$xmin), function(i) {
      
      grid::roundrectGrob(
        coords$xmin[i], coords$ymax[i],
        width = (coords$xmax[i] - coords$xmin[i]),
        height = (coords$ymax[i] - coords$ymin)[i],
        r = radius,
        default.units = "native",
        just = c("left", "top"),
        gp = grid::gpar(
          col = coords$colour[i],
          fill = alpha(coords$fill[i], coords$alpha[i]),
          lwd = coords$size[i] * .pt,
          lty = coords$linetype[i],
          lineend = "butt"
        )
      )
      
    }) -> gl
    
    grobs <- do.call(grid::gList, gl)
    
    ggplot2:::ggname("geom_rrect", grid::grobTree(children = grobs))
    
  },
  draw_key = ggplot2::draw_key_polygon
)
