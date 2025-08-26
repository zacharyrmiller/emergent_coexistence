Code to reproduce simulations and generate all figures from "Multispecies coexistence emerges from pairwise exclusions in communities with competitive hierarchy" (Miller & Max 2025) published in Ecology Letters.

Functions.R contains all functions used to sample model parameters, analyze simulation data, and visualize results
run_simulations.R runs a long loop to sample many GLV communities for random community simulations (corresponding to Figs. 4 and 5 in the text)
make_figures.R contains code to generate all figures in the text (including SI figures)

Some figures (4, 5, S4-S6) require results from numerical simulations. The simulation set used in the manuscript is provided in the file results_sample_A_solve_r_n4-10_10000_reps.csv (so that it is not necessary to re-run the run_simulations script in order to reproduce figures in the paper).
