n <- 6 # number of species
method <- "sample_A_solve_r" # method for sampling A and r
n_reps <- 500 # number of parameter combinations

results <- tibble(rep = numeric(),
                  n = numeric(),
                  exclusion_frac = numeric(),
                  coexistence_frac = numeric(),
                  hierarchy_score = numeric(),
                  loops = logical(),
                  positive_effects_frac = numeric(),
                  directed_triangles = numeric())


for(i in 1:n_reps){
  
  x <- sample_SAD(n, lambda = 1)
  Ar_pair <- sample_A_r_pair(x, method = method)
  
  A <- Ar_pair$A
  r <- Ar_pair$r
  pairwise_outcomes <- build_pairwise_outcome_matrix(A, r)
  
  results <- results %>% add_row(rep = i,
                                 n = n,
                                 exclusion_frac = count_fraction_excluded(pairwise_outcomes),
                                 coexistence_frac = count_fraction_coexist(pairwise_outcomes),
                                 hierarchy_score = compute_hierarchy_score(pairwise_outcomes),
                                 loops = test_for_loops(pairwise_outcomes),
                                 positive_effects_frac = count_fraction_positive_net_effects(A),
                                 directed_triangles = count_intransitive_triangles(pairwise_outcomes)
                                 )
}

results %>% ggplot() + 
  aes(x = hierarchy_score) + 
  geom_histogram() + 
  theme_bw()

results %>% ggplot() + 
  aes(x = exclusion_frac, y = positive_effects_frac) + 
  geom_jitter() + 
  theme_bw()
