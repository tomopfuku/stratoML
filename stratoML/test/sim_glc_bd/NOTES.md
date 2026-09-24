## MICRASTER

### MCMC to reconstruct rates

python ../../main_single_tree_MCMC_parallel.py \
  --trees best_mfc2.tre \
  --traits micrasteridae-character_matrix.fa \
  --strat-data micrasteridae-stratigraphy_abs.stag.csv \
  --strat-model hr97 \
  --morph-model mfc2 \
  --generations 50000 \
  --max-threads 8

### Plot rates

python ../../plot_clado_loss_rates.py micrasteridae-character_matrix.fa_mcmc_samples.csv 0.1062790198482735

### Model comparison

python ../../main_single_tree_like3.py --tree best_mfc2.tre --traits micrasteridae-character_matrix.fa --strat-data micrasteridae-stratigraphy_abs.stag.csv --strat-model hr97 --morph-model mfc2

BDS rates: [0.10627898767685667, 0.11350275485899007, 1.388630828099045]
no jump params: [0.04241621 0.09571556 0.06693704]
GAIN: 0.042416214640003336
LOSS: 0.09571555581392088
LAMBDA SUB: 0.07114000549888232
no jump AIC: 6375.100937658214
no cladogenesis params: [0.02672896 0.15692895]
no cladogenesis AIC: 6382.510072111914
cladogenesis AIC weight: 0.9759802801216785
no cladogenesis AIC weight: 0.024019719878321487



## BARY


### MCMC to reconstruct rates

python ../../main_single_tree_MCMC_parallel.py \
  --trees bary.best_mfc.tre \
  --traits botryocrinidae_gahn_kammer.form.fa \
  --strat-data bary_ranges.csv \
  --strat-model hr97 \
  --morph-model mfc2 \
  --generations 50000 \
  --max-threads 8

### Plot rates

python ../../plot_clado_loss_rates.py botryocrinidae_gahn_kammer.form.fa_mcmc_samples.csv 0.06006928815485987


### Model comparison

python ../../main_single_tree_like3.py --tree bary.best_mfc.tre --traits botryocrinidae_gahn_kammer.form.fa --strat-data bary_ranges.csv --strat-model hr97 --morph-model mfc2

BDS rates: [0.06006928778990395, 0.06829334778031942, 1.4615270767325743]
no jump params: [0.01099953 0.01981585 0.07898446]
GAIN: 0.010999532842649729
LOSS: 0.019815853125057518
LAMBDA SUB: 0.047445403160177806
no jump AIC: 2456.4604225006456
no cladogenesis params: [0.01824528 0.08071898]
no cladogenesis AIC: 2481.7127830030663
cladogenesis AIC weight: 0.9999967151301822
no cladogenesis AIC weight: 3.284869817839431e-06



## DEND

python ../../main_single_tree_MCMC_parallel.py \
  --trees dinky.tre \
  --traits dend.fa \
  --strat-data dend_fadlad.csv \
  --strat-model hr97 \
  --morph-model mfc2 \
  --generations 50000 \
  --max-threads 8


### Plot histogram of clado vs anagenetic loss rates

python ../../plot_clado_loss_rates.py dend.fa_mcmc_samples.csv 0.4393345299822384

### Model comparison:

python ../../main_single_tree_like3.py --tree dinky.tre --traits dend.fa --strat-data dend_fadlad.csv --strat-model hr97 --morph-model mfc2

BDS rates: [0.4393346315470331, 0.4626443069085857, 2.421481649977755]
no jump params: [0.3477552  0.55296046 0.04739941]
GAIN: 0.34775519762542295
LOSS: 0.552960457348989
LAMBDA SUB: 0.20824204173120472
no jump AIC: 3950.6825618819676
no cladogenesis params: [0.34826601 0.72244699]
no cladogenesis AIC: 3953.2569619337496
cladogenesis AIC weight: 0.7836728874633027
no cladogenesis AIC weight: 0.2163271125366973



### SIMULATIONS

python run_sim_likelihood_batch.py --tree best_mfc2.tre --ranges micrasteridae-stratigraphy_abs.stag.csv

python run_no_clado_false_positive_batch.py --tree best_mfc2.tre --ranges micrasteridae-stratigraphy_abs.stag.csv --num-sims 100




python ../../main_glc_asr.py best_mfc2.tre micrasteridae-character_matrix.fa micrasteridae-stratigraphy_abs.stag.csv hr97 glc
