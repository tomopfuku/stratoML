## MICRASTER

MCMC_MAX_THREADS=8 python ../../main_single_tree_MCMC_parallel.py best_mfc2.tre micrasteridae-character_matrix.fa micrasteridae-stratigraphy_abs.stag.csv  hr97 mfc2 50000

python ../../plot_clado_loss_rates.py micrasteridae-character_matrix.fa_mcmc_samples.csv 0.1062790198482735


## BARY

MCMC_MAX_THREADS=8 python ../../main_single_tree_MCMC_parallel.py bary.best_mfc.tre botryocrinidae_gahn_kammer.form.fa bary_ranges.csv  hr97 mfc2 10000

python ../../plot_clado_loss_rates.py botryocrinidae_gahn_kammer.form.fa_mcmc_samples.csv 0.06006928815485987


## DEND



python ../../plot_clado_loss_rates.py dend.fa_mcmc_samples.csv 0.4393345299822384
