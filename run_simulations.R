# Simulation functions to accompany "Multispecies coexistence emerges from pairwise competitive exclusions in communities with competitive hierarchy" by Z.R. Miller and D. Max
# Questions about code? Contact: zachary.miller@yale.edu

source("./functions.R")

### Run numerical simulations of random GLV communities (corresponding to Figs. 4 and 5 and S4-6) ###

## Outputs are: 
## results = a tibble of summary statistics used in the analysis
## A_list = a list of all interaction matrices
## r_list = a list of all growth rate vectors

## Set-up
method <- "sample_A_solve_r" # method for sampling A and r (sample_A_solve_r is used in the manuscript)
n_reps <- 1000 # number of parameter combinations (10000 is used for manuscript figures; this will take 10+ hours. Run with n_reps = 100 or fewer for fast results)
community_sizes <- 4:10 # note: sampling for n > 10 becomes VERY slow

# initialize a tibble for simulation results
results <- tibble(rep = numeric(), # replicate id for each n
                  n = numeric(), # number of species
                  realized_n = numeric(), # number of persistent species (= n unless method = "assemble_from_pool")
                  exclusion_frac = numeric(), # fraction of pairwise outcomes that are competitive exclusion
                  coexistence_frac = numeric(), # fraction of pairwise outcomes that are coexistence (note 1 - exclusion_frac - coexistence_frac = bistability_frac)
                  loops = logical(), # are there any (intransitive) loops in the pairwise outcome network?
                  directed_triangles = numeric(), # number of directed (intransitive) triangles
                  total_triangles = numeric(), # number of triangles of any orientation (e.g. set of 3 species where all pairwise outcomes are competitive exclusion)
                  mean_ratio = numeric(), # ratio of mean inter- vs. mean intraspecific interaction strengths
                  positive_effects_frac = numeric()) # fraction of indirect effects (elements of A^-1) that are positive (excluding diagonal elements)
                  
A_list <- vector(mode = "list", length = length(community_sizes) * n_reps)
r_list <- vector(mode = "list", length = length(community_sizes) * n_reps)
counter <- 1

for(n in community_sizes){ # loop over communities of different sizes
  
  print(n)
  
  for(i in 1:n_reps){
    
    x <- sample_SAD(n, lambda = 1)
    Ar_pair <- sample_A_r_pair(x,
                               adjust_symmetry = FALSE, # symmetry not adjusted for this manuscript
                               method = method,
                               require_stable = TRUE
    )

    A <- Ar_pair$A
    r <- Ar_pair$r
    
    pairwise_outcomes <- build_pairwise_outcome_matrix(A, r)

    results <- results %>% add_row(rep = i,
                                   n = n,
                                   realized_n = length(r),
                                   exclusion_frac = count_fraction_excluded(pairwise_outcomes),
                                   coexistence_frac = count_fraction_coexist(pairwise_outcomes),
                                   loops = test_for_loops(pairwise_outcomes),
                                   directed_triangles = count_intransitive_triangles(pairwise_outcomes),
                                   total_triangles = count_total_triangles(pairwise_outcomes),
                                   mean_ratio = calculate_inter_vs_intra_ratio(A),
                                   positive_effects_frac = calculate_positive_effects_frac(A)
    )
    
    # 
    A_list[[counter]] <- A
    r_list[[counter]] <- r
    counter <- counter + 1
    
  }
}
