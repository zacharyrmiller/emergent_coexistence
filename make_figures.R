# Simulation functions to accompany "Multispecies coexistence emerges from pairwise competitive exclusions in communities with competitive hierarchy" by Z.R. Miller and D. Max
# Questions about code? Contact: zachary.miller@yale.edu

library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(pammtools)
library(deSolve)

source("./functions.R")

### Make all manuscript figures ###

##### for Fig 1 (conceptual figure; full figure assembled in Illustrator) #####

n <- 4
r <- c(0.5, 1.2, 1.7, 2.15) # r = A %*% c(1, 0.4, 0.6, 0.3) so eq. values approximately match figure

A <- lower.tri(diag(n)) + 0
diag(A) <- 0.5

parameters <- list(A = -as.vector(A), r = r, n = n)
times <- seq(0, 60, by = 0.1) 

x <- c(0.4, 0.3, 0.2, 0.1) # initial conditions for panel d 
# x <- c(0, 0.3, 1.2, 0) # initial conditions for panel f (left)
# x <- c(0, 0.3, 0, 1.2) # initial conditions for panel f (right)

out <- ode(y = x, time = times, func = glv, parms = parameters, method = "ode45")

p <- out %>% as.data.frame() %>%
  pivot_longer(-time) %>%
  ggplot() + aes(x = time, y = value,
                 group = name, color = as.factor(name)) + 
  geom_line(size = 1) + 
  ylim(c(0.01, 1.2)) + 
  xlab("Time") + ylab("Density") +
  scale_colour_manual(values = c("#1c4587ff","#3d85c6ff","#a4c2f4ff","#cfe2f3ff")) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())

# ggsave("./full_community.png", p, device = "png", 
#       width = 2.2, height = 1.2, 
#       units = "in", dpi = 400)


##### Fig 2 #####

n <- 9

## Generate three illustrative growth rate sequences (as discussed in text)

x <- 0.5^(1:n)
r_concave <- make_r_from_x(x)

r1 <- 1 / 2^(n-1)
r_convex <- r1 * 2^(0:(n-1))

const <- 1 / (n - 0.5)
r_linear <- make_r_from_x(rep(const, n))

## Make plot panels

# combine growth rates for plotting 
r_df <- tibble(rank = rep(1:n, 3), 
               r = c(r_concave, r_convex, r_linear),
               type = rep(1:3, each = n))

# plot different growth rate sequences vs. competitive rank
p_a <- r_df %>% ggplot() + aes(x = rank, y = r, 
                               color = as.factor(type), shape = as.factor(type)) + 
  geom_point(size = 2) +
  scale_shape_manual(values = c(19, 17, 4)) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_x_continuous(breaks = c(1, ceiling(n/2), n), name = "Competitive rank") + 
  scale_y_continuous(breaks = c(0, 1), limits = c(0,1), name = "Growth rate (r)") + 
  theme_classic() + 
  theme(legend.position = "none")

# make pairwise outcome plots for each set of growth rates
p_1 <- make_pairwise_outcomes_plot(r_concave, color = brewer.pal(3, "Dark2")[1])
p_2 <- make_pairwise_outcomes_plot(r_convex, color = brewer.pal(3, "Dark2")[2])
p_3 <- make_pairwise_outcomes_plot(r_linear, color = brewer.pal(3, "Dark2")[3])

# plot minimal coexistence fraction vs. number of species
p_c <- tibble(x = seq(3, 100, by = 3),
              y = (x - 1) / choose(x, 2)) %>%
  ggplot() + aes(x = x, y = y) + 
  geom_point(color = brewer.pal(3, "Dark2")[1]) + 
  annotate(geom = "point", x = n, y = (n - 1) / choose(n,2), 
           color = brewer.pal(3, "Dark2")[1], size = 3, pch = 5) +
  scale_y_log10() + 
  xlab("Number of species") + ylab("Fraction of pairs coexisting") +
  theme_classic()

# combine panels
p <- ggarrange(p_a,
               ggarrange(p_3, p_1, p_2, ncol = 2, nrow = 2) + 
                 theme(plot.margin = unit(c(0, 4, 13, 6), "pt")),
               p_c,
               ncol = 3, nrow = 1, labels = "auto"
)

# add annotations (position values are chosen by hand for a nice appearance)
p <- p + annotate("rrect", xmin = 0.515, xmax = 0.53, ymin = 0.45, ymax = 0.495,
                  color = brewer.pal(3, "Dark2")[1], fill = brewer.pal(3, "Dark2")[1],
                  radius = grid::unit(2, "pt")) +
  annotate("text", x = 0.592, y = 0.475, label = "Coexistence") + 
  annotate("rrect", xmin = 0.515, xmax = 0.53, ymin = 0.39, ymax = 0.435,
           color = brewer.pal(3, "Dark2")[1], fill = "white",
           radius = grid::unit(2, "pt")) +
  annotate("text", x = 0.58, y = 0.413, label = "Exclusion") + 
  annotate("text", x = 0.555, y = 0.09, label = "Competitive rank") + 
  annotate("segment", x = 0.505, y = 0.3, xend = 0.505, yend = 0.15,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) + 
  annotate("segment", x = 0.46, y = 0.06, xend = 0.41, yend = 0.06,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))

# ggsave("./figure_2.png", p, device = "png", 
#       width = 8, height = 2.7, 
#       units = "in", dpi = 600)


##### Fig 3 #####

n <- 40

## Sample random growth rates from (sorted) uniform distribution (Fig. 3a-c) or triangular distribution (Fig. 3d-f)

sample_triangle <- function(n){
  # sample from triangular distribution using inverse transform sampling
  
  raw <- runif(n)
  transformed <- sqrt(raw)
  return(transformed)
}

# uniform pool
set.seed(88)
r_pool <- sort(runif(n))

# triangular pool
# set.seed(34)
# r_pool <- sort(sample_triangle(n))

# build interaction matrix according to growth-competition model
A <- lower.tri(diag(n)) + 0
diag(A) <- 0.5

parameters <- list(A = -as.vector(A), r = r_pool, n = n)
times <- seq(0, 10000, by = 0.1) 
x <- rep(0.1, n) # initial all species at equal abundances

out <- ode(y = x, time = times, func = glv, parms = parameters, method = "ode45")

# plot time series of community dynamics
p_b <- out %>% as.data.frame() %>%
  pivot_longer(-time) %>%
  ggplot() + aes(x = time, y = value,
                 group = name, color = as.numeric(name)) + 
  geom_line(size = 1) + 
  xlab("Time") + ylab("Density") +
  scale_y_log10(limits = c(10^-7, NA), labels = scales::label_log(), expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0.01), 
                     breaks = c(0, 2000, 4000, 6000, 8000)) + 
  scale_colour_gradientn(colours = c("#1c4587ff","#3d85c6ff","#a4c2f4ff","#cfe2f3ff")) +
  theme_classic() + 
  theme(legend.position = "none")

# plot growth rate vs. competitive rank
persistent <- out[nrow(out), -1] > 10^-6 # which species persisted in the final community
p_a <- tibble(rank = 1:n, r_pool = r_pool, persistent = persistent) %>%
  ggplot() + 
  aes(x = rank, y = r_pool, color = as.numeric(rank), 
      shape = persistent) + 
  geom_point(size = 2) +
  scale_color_gradientn(colours = c("#1c4587ff","#3d85c6ff","#a4c2f4ff","#cfe2f3ff")) + 
  scale_x_continuous(name = "Competitive rank") + 
  scale_y_continuous(name = "Growth rate (r)", limits = c(0, 1)) + 
  scale_shape_manual(values = c(4, 19)) +
  theme_classic() + 
  theme(legend.position = "none")

# plot pairwise outcomes among species in the assembled community
r_assembled <- r_pool[persistent] # the set of persisting species
p_c <- make_pairwise_outcomes_plot(r_assembled, color = "#3d85c6ff", axis.title = TRUE)

p_a <- p_a + theme(axis.title.x = element_blank())
p_b <- p_b + theme(axis.title.x = element_blank())
p_c <- p_c + theme(axis.title.x = element_blank())

# combine panels
p <- ggarrange(p_a, p_b, p_c,
               nrow = 1,
               labels = c("a", "b", "c"))

# add annotations (position values are chosen by hand for a nice appearance)
p <- p + 
  annotate("rrect", xmin = 0.7, xmax = 0.7135, ymin = 0.8, ymax = 0.84,
           color = "#3d85c6ff", fill = "#3d85c6ff",
           radius = grid::unit(2, "pt")) +
  annotate("text", x = 0.782, y = 0.82, label = "Coexistence") + 
  annotate("rrect", xmin = 0.7, xmax = 0.7135, ymin = 0.73, ymax = 0.77,
           color = "#3d85c6ff", fill = "white",
           radius = grid::unit(2, "pt")) +
  annotate("text", x = 0.77, y = 0.75, label = "Exclusion")

p_top <- p
# p_bottom <- p

# ggsave("./figure_3_SI.png", ggarrange(p_top, p_bottom, nrow = 2), 
#       device = "png", width = 7, height = 4.25, 
#       units = "in", dpi = 500)


##### Fig 4 #####

# read in results from numerical simulations of random GLV communities
results <- read_csv(".results_sample_A_solve_r_n4-10_10000_reps.csv")

# plot distribution of coexistence fractions (violin plot) for each n
p <- results %>% ggplot() + 
  aes(x = as.factor(n), y = coexistence_frac, color = as.factor(n), group = n) + 
  geom_violin(draw_quantiles = c(0.25, 0.75), adjust = 2, linetype = "dashed") + # first plot quartiles
  geom_violin(draw_quantiles = c(0.5), adjust = 2, fill = "transparent") + # overlay median
  scale_color_manual(values = brewer.pal(7, "Dark2"))+
  xlab("Number of species") + ylab("Coexistence fraction") + 
  theme_classic() + 
  theme(legend.position = "none")

# ggsave("./figure_4.png", p, 
#       device = "png", width = 5, height = 3, 
#       units = "in", dpi = 500)


##### Fig 5 #####

show_n <- c(5, 8, 10) # which values of n should be plotted

## panel a: Mean number of directed triangles observed vs. null distribution (after Chang et al. Science 2023; Fig. 3)

# generate null distributions of number of directed (intransitive) triangles, conditioned on total triangles
null_reps <- 10000
tt <- results %>% group_by(n) %>% summarize(tt = sum(total_triangles)) %>% pull(tt) # total triangles by n
null_distribution = data.frame(n = rep(4:10, each = null_reps), 
                               null_directed_triangles = unlist(sapply(tt, function(x) rbinom(null_reps, x, 1/4), # null distribution is Binom(tt, 1/4)
                                                                       simplify = FALSE)))

# plot histogram of null distribution vs observed
pa <- null_distribution %>% 
  filter(n %in% show_n) %>%
  ggplot() + 
  aes(x = null_directed_triangles / 10000) +
  geom_histogram(aes(y =..density..), binwidth = 0.002, fill = "dodgerblue") +
  geom_vline(data = results %>% 
               filter(n %in% show_n) %>% 
               group_by(n) %>% summarize(dt = sum(directed_triangles)),
             aes(xintercept = dt / 10000), linetype = "dashed", linewidth = 1, color = "red") +
  facet_wrap(.~n, scales = "free", nrow = 1,
             labeller = labeller(n = function(x) paste0("n = ", x))) + 
  expand_limits(x = 0) + 
  xlab("Intransitive triplets") + ylab("Density") + 
  theme_classic() + 
  theme(legend.position = "none", panel.border = element_rect(fill = "NA", linewidth = 1))

## panel b: mean number of directed triangles vs. exclusion fraction for simulations and null model 

# generate null distributions of number of directed triangles for each observed exclusion fraction and n combination
null_reps <- 10000
null_t_ratio <- results %>% select(c(n, exclusion_frac)) %>%
  filter(exclusion_frac > 0, n %in% show_n) %>% distinct() # distinct (non-zero) exclusion fractions observed in simulations

for(i in 1:nrow(null_t_ratio)){ # for each distinct exclusion fraction (and n), generate many random networks 
  
  null_dist <- rep(NA, null_reps)
  for(j in 1:null_reps){
    null_matrix <- build_random_directed_network(null_t_ratio$n[i], null_t_ratio$exclusion_frac[i])
    null_dist[j] <- count_intransitive_triangles(null_matrix) 
  }
  
  # save the mean and quartiles for the null distribution
  null_t_ratio$mean[i] <- mean(null_dist)
  null_t_ratio$lower[i] <- quantile(null_dist, 0.25)
  null_t_ratio$upper[i] <- quantile(null_dist, 0.75)
}

# plot directed triangles vs. exclusion fraction
pb <- results %>% group_by(n, exclusion_frac) %>%
  filter(n %in% show_n) %>%
  summarise(lower = quantile(directed_triangles, 0.25),
            upper = quantile(directed_triangles, 0.75),
            mean = mean(directed_triangles), 
            k = n()) %>%
  ggplot() + aes(x = exclusion_frac) + 
  geom_step(data = null_t_ratio %>% # null model mean
              filter(n %in% show_n), aes(y = mean), direction = "mid", color = "dodgerblue", size = 1.2) + 
  geom_stepribbon(data = null_t_ratio, aes(ymin = lower, ymax = upper), # null model interquartile range
                 fill = "dodgerblue", alpha = 0.4, direction = "mid") +
  geom_point(aes(y = mean, size = k > 3), color = "red") + # simulation mean
  geom_errorbar(aes(ymin = ifelse(k > 3, lower, mean), # simulation interquartile range (only if > 3 observations)
                    ymax = ifelse(k > 3, upper, mean)), 
                width = 0, color = "red") +
  scale_size_manual(values = c(0.75, 1.5)) + 
  scale_x_continuous(breaks = seq(0, 1, 0.2)) + 
  scale_y_continuous(breaks = seq(0, 20, 4)) +
  facet_wrap(.~n, scales = "free", nrow = 1) + 
  xlab("Exclusion fraction") + ylab("Intransitive triplets") + 
  theme_classic() + 
  theme(legend.position = "none", 
        panel.border = element_rect(fill = "NA", linewidth = 1),
        strip.background = element_blank(),
        strip.text = element_blank())

### panel c : probability of hierarchy (no intransitive loops of any length) vs. exclusion fraction for simulations and null model

# count probability of loop-less networks in null models for each observed exclusion fraction and n combination
null_reps <- 10000
null_loop_ratio <- results %>% select(c(n, exclusion_frac)) %>% 
  filter(exclusion_frac > 0) %>% distinct()  # distinct (non-zero) exclusion fractions observed in simulations

for(i in 1:nrow(null_loop_ratio)){ # for each distinct exclusion fraction (and n), generate many random networks 
  
  null_dist <- rep(0, null_reps)
  for(j in 1:null_reps){
    null_matrix <- build_random_directed_network(null_loop_ratio$n[i], null_loop_ratio$exclusion_frac[i])
    null_dist[j] <- test_for_loops(null_matrix)
  }
  
  # save probability of intransitive loops across null model ensemble
  null_loop_ratio$mean[i] <- mean(null_dist)
}

# plot probability of hierarchy vs. exclusion fraction 
pc <- results %>% group_by(n, exclusion_frac) %>% 
  filter(n %in% show_n) %>% 
  mutate(k = n(),
         mean_loop_ratio = mean(loops), 
         mean_loop_ratio = ifelse(k > 3, mean_loop_ratio, NA)
  ) %>%
  ggplot() + aes(x = exclusion_frac) + 
  geom_step(aes(y = 1 - mean_loop_ratio), direction = "mid", color = "red", size = 1.2) + 
  geom_step(data = null_loop_ratio %>% filter(n %in% show_n), aes(y = 1 - mean), 
            direction = "mid", color = "dodgerblue", size = 1.2) + 
  facet_wrap(.~n, scales = "free", nrow = 1) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  xlab("Exclusion fraction") + ylab("Probability of hierarchy") + 
  #geom_text(data = data.frame(n = 10, exclusion_frac = 0.45), y = 0.7, label = "GLV", color = "red", size = 5) + 
  #geom_text(data = data.frame(n = 10, exclusion_frac = 0.2), y = 0.2, label = "Null", color = "dodgerblue", size = 5) + 
  theme_classic() + 
  theme(legend.position = "none", 
        panel.border = element_rect(fill = "NA", linewidth = 1),
        strip.background = element_blank(),
        strip.text = element_blank())

# combine panels
p <- ggarrange(plotlist = list(pa, pb, pc), nrow = 3, labels = "auto")

ggsave("../figure_5.png", p,
      device = "png", width = 7.5, height = 7.5,
      units = "in", dpi = 500)



##### SI Figures ##### 

##### Fig S1 ##### 

n_reps <- 10^7

results <- tibble(n = numeric(), type = character(), probability = numeric())
for(n in 4:15){
  
  count_all_pairs <- 0
  count_maximal_exclusion <- 0
  for(i in 1:n_reps){
    
    x <- rexp(n) # sample equilibrium abundances iid
    
    all_pairs <- all(x[-(1:2)] > 2 * cumsum(x)[-((n-1):n)]) # check condition for all pairs to coexist
    maximal_exlusion <- 2 * x[1] > 2 * sum(x[-(1:2)]) - x[n] # check condition for maximal exclusion
    
    # add to running count for both
    count_all_pairs <- count_all_pairs + all_pairs
    count_maximal_exclusion <- count_maximal_exclusion + maximal_exlusion
  }
  
  results <- results %>% 
    add_row(n = n, type = "All pairs coexist", probability = count_all_pairs / n_reps) %>%
    add_row(n = n, type = "Maximal exclusion", probability = count_maximal_exclusion / n_reps)
}

# plot probability of each outcome vs. n
p <- results %>% 
  mutate(probability = ifelse(probability == 0, 10^-10, probability)) %>% 
  ggplot() + 
  aes(x = n, y = probability, group = type, color = type) + 
  geom_point() + geom_line() + 
  geom_line(data = data.frame(n = 4:15, # plot theoretical upper bound (Eq. S15)
                              probability = 1 / (factorial(floor((4:15)/2)) * factorial(ceiling((4:15)/2))), 
                              type = "All pairs coexist"), 
            linetype = "dashed") + 
  scale_y_log10() +
  coord_cartesian(ylim = c(10^-6, 1)) + 
  xlab("Number of species (n)") + 
  ylab("Probability") + 
  theme_classic() + 
  theme(legend.title = element_blank(), legend.position = "inside",
        legend.position.inside = c(0.75, 0.8))

# ggsave("./figure_S1.png", p, 
#       device = "png", width = 3.5, height = 3.5, 
#       units = "in", dpi = 500)


##### Fig S2 #####

p_vec <- c(-0.9, -0.5, -0.25, 0, 1, 2) # exponents to plot
n_vec <- 2^(2:11) # values of n to plot

results <- tibble(p = numeric(), n = numeric(), CF = numeric(), approx_CF = numeric())
for(n in n_vec){
  
  # build interaction matrix according to growth-competition model
  A <- lower.tri(diag(n)) + 0
  diag(A) <- 0.5
  
  for(p in p_vec){
    
    x <- ((1:n) / n)^p # abundances follow power law
    r <- 2 * c(0, cumsum(x)[-n]) + x # generate corresponding growth rates
    
    pairwise_outcomes <- build_pairwise_outcome_matrix(A, r)
    CF <- count_fraction_coexist(pairwise_outcomes) # actual coexistence fraction
    approx_CF <- 0.5^(1 / (p + 1)) # approximate coexistence fraction (Eq. S31)
    
    results <- results %>% 
      add_row(p = p, n = n, CF = CF, approx_CF = approx_CF)
  }
}

# plot actual and approximate coexistence fraction
p <- results %>% ggplot() + 
  aes(x = n, color = as.factor(p), group = p) + 
  geom_line(aes(y = approx_CF), linetype = "dashed") + 
  geom_line(aes(y = 2 / n + (1 - 2 / n) * approx_CF)) + # "adjusted" approximation (Eq. S28 applied to CF)
  geom_point(aes(y = CF)) + 
  scale_x_log10(name = "Number of species (n)") + ylab("Coexistence fraction") + 
  scale_color_viridis_d(name = "p = ") + 
  theme_classic()

# ggsave("./figure_S2.png", p, 
#       device = "png", width = 4, height = 3, 
#       units = "in", dpi = 500)


##### Fig S3 #####

n <- 40

## Sample random growth rates from (sorted) uniform distribution (Fig. S3a-c) or triangular distribution (Fig. S3d-f)

# uniform pool
set.seed(88)
r_pool <- sort(runif(n))

# triangular pool
# set.seed(34)
# r_pool <- sort(sample_triangle(n))

# convert growth rates to competition and local extinction rates, and build interaction matrix according to competition-colonization model
m <- 1
c_pool <- r_pool + m
A <- diag(c_pool) %*% lower.tri(diag(n)) + lower.tri(diag(n)) %*% diag(c_pool)
diag(A) <- c_pool

parameters <- list(A = -as.vector(A), r = r_pool, n = n)
times <- seq(0, 2*10000, by = 0.1) 
x <- rep(0.1, n) # initial all species at equal abundances

out <- ode(y = x, time = times, func = glv, parms = parameters, method = "ode45")

# plot time series of community dynamics 
p_b <- out %>% as.data.frame() %>%
  pivot_longer(-time) %>%
  ggplot() + aes(x = time, y = value,
                 group = name, color = as.numeric(name)) + 
  geom_line(size = 1) + 
  xlab("Time") + ylab("Density") +
  scale_y_log10(limits = c(10^-7, NA), labels = scales::label_log(), expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0.01), 
                     breaks = c(0, 4000, 8000, 12000, 16000)) + 
  scale_colour_gradientn(colours = c("#1c4587ff","#3d85c6ff","#a4c2f4ff","#cfe2f3ff")) +
  theme_classic() + 
  theme(legend.position = "none")

# plot growth rate vs. competitive rank
persistent <- out[nrow(out), -1] > 10^-6
p_a <- tibble(rank = 1:n, r_pool = r_pool, persistent = persistent) %>%
  ggplot() + 
  aes(x = rank, y = r_pool, color = as.numeric(rank), 
      shape = persistent) + 
  geom_point(size = 2) +
  scale_color_gradientn(colours = c("#1c4587ff","#3d85c6ff","#a4c2f4ff","#cfe2f3ff")) + 
  scale_x_continuous(name = "Competitive rank") + 
  scale_y_continuous(name = "Growth rate (r)", limits = c(0, 1)) + 
  scale_shape_manual(values = c(4, 19)) +
  theme_classic() + 
  theme(legend.position = "none")

# plot pairwise outcomes among species in the assembled community
r_assembled <- r_pool[persistent]
p_c <- make_pairwise_outcomes_plot(r_assembled, m = 1, color = "#3d85c6ff", axis.title = TRUE)

p_a <- p_a + theme(axis.title.x = element_blank())
p_b <- p_b + theme(axis.title.x = element_blank())
p_c <- p_c + theme(axis.title.x = element_blank())

# combine panels
p <- ggarrange(p_a, p_b, p_c,
               nrow = 1,
               labels = c("a", "b", "c"))

# add annotations (position values are chosen by hand for a nice appearance)
p <- p + 
  annotate("rrect", xmin = 0.7, xmax = 0.7135, ymin = 0.8, ymax = 0.84,
           color = "#3d85c6ff", fill = "#3d85c6ff",
           radius = grid::unit(2, "pt")) +
  annotate("text", x = 0.782, y = 0.82, label = "Coexistence") + 
  annotate("rrect", xmin = 0.7, xmax = 0.7135, ymin = 0.73, ymax = 0.77,
           color = "#3d85c6ff", fill = "white",
           radius = grid::unit(2, "pt")) +
  annotate("text", x = 0.77, y = 0.75, label = "Exclusion")

p_top <- p
# p_bottom <- p

# ggsave("./figure_S3.png", ggarrange(p_top, p_bottom, nrow = 2), 
#       device = "png", width = 7, height = 4.25, 
#       units = "in", dpi = 500)


##### Fig S4 #####

# read in results from numerical simulations of random GLV communities
results <- read_csv(".results_sample_A_solve_r_n4-10_10000_reps.csv")

p <- results %>%
  ggplot() + 
  aes(x = log(mean_ratio), y = coexistence_frac) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(.~n, scales = "free") + 
  geom_smooth(color = "red", fill = "red", alpha = 0.2) + 
  ylim(c(0,1)) + 
  xlab("Log ratio of mean inter- vs. mean intraspecific competition") + ylab("Coexistence fraction") + 
  theme_classic()

# ggsave("../figure_5_SI.png", plot = p, 
#       device = "png", width = 6, heigh = 6, units = "in", dpi = 500)


##### Fig S5 #####

# generate null distribution of interaction matrices

n_reps <- 10000 
community_sizes <- 4:10

positive_effects_null <- c()
null_A_list <- vector(mode = "list", length = length(community_sizes) * n_reps)
counter <- 1

for(n in community_sizes){ # loop over communities of different sizes
  
  print(n)
  
  for(i in 1:n_reps){
    
    A <- matrix(rexp(n^2), nrow = n)
    
    positive_effects_null <- append(positive_effects_null, calculate_positive_effects_frac(A))
    null_A_list[[counter]] <- A
    counter <- counter + 1
    
  }
}


p <- results %>% mutate(positive_effects_null = positive_effects_null) %>% 
  ggplot() + 
  geom_histogram(aes(x = positive_effects_frac), 
                 color = "red", fill = "red", alpha = 0.25) + 
  geom_histogram(aes(x = positive_effects_null), 
                 color = "dodgerblue", fill = "dodgerblue", alpha = 0.25) +
  geom_vline(xintercept = 0.5, linetype = "dashed") + 
  geom_vline(data = df %>% group_by(n) %>% 
               summarize(mean = mean(positive_effects_frac)),
             aes(xintercept = mean),
             color = "red") + 
  geom_vline(data = df %>% group_by(n) %>% 
               summarize(mean = mean(positive_effects_null)),
             aes(xintercept = mean),
             color = "dodgerblue") + 
  #annotate(geom = "text", x = 0.59, y = 11.5, label = "Coexisting", color = "red") + 
  #annotate(geom = "text", x = 0.42, y = 7, label = "Null", color = "dodgerblue") + 
  xlab("Fraction of positive net effects") + ylab("Count") +
  facet_wrap(.~n, scales = "free") + 
  theme_classic()

# ggsave("../figure_6_SI.png", plot = p, 
#       device = "png", width = 7, heigh = 5.5, units = "in", dpi = 500)

