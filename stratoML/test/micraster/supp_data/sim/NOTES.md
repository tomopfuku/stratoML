# Commands for running and analyzing simulations

NOTE: all of these commands assume you are in the directory `stratoML/stratoML/test/sim` from the stratoML github repository. The files are provided here for posterity, but dependencies from stratoML are needed to run them.

## Simulate 30 trait matrices
python simulate_fastas.py test.tre strat.tab 30 3 sim_matrices_30_3/

## Simulate 60 trait matrices
python simulate_fastas.py test.tre strat.tab 60 3 sim_matrices_60_3/

## Estimate rates from simulated matrices and create scatterplots

python estimate_mfc_rates.py test.tre sim_matrices_30_3 strat.tab hr97

python estimate_mfc_rates.py test.tre sim_matrices_60_3 strat.tab hr97

## Calculate likelihood of all trees for each simulated matrix and plot accuracy
python calc_tree_like.py all_test.trees sim_matrices_30_3 strat.tab hr97

python calc_tree_like.py all_test.trees sim_matrices_60_3 strat.tab hr97

## Calculate number of character states across all simulated matrices

python calc_num_states.py sim_matrices_30_3/

python calc_num_states.py sim_matrices_60_3/
