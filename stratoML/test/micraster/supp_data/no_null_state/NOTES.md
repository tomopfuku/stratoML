## Files

micrasteridae_supp.xlsx -- this file contains all of the original morphological character data and stratigraphic range data. The data is organized into several sheets as follows:

    character_matrix -- this is a representation of the character matrix in spreadsheet format. 

    trait_descriptions -- descriptions of each morphological character and its character states

    species_temporal_ranges -- this approximates the species temporal ranges in discrete bins.

    timetable_for_bins -- start and end times of each bin in absolute time. 

micrasteridae-character_matrix.16.fa -- character matrix in FASTA format as used for the analyses.

micrasteridae-stratigraphy_abs.csv -- micraster stratigraphic ranges represented in absolute time

micrasteridae-stratigraphy_abs.stag.csv -- micraster stratigraphic ranges represented in absolute time with simultaneous FADs trivially adjusted to avoid numerical issues

top_trees_mfc2 -- all highest scoring Micraster topologies. This pool of trees was used as an approximation of the likelihood surface to calculate branch support.

best_mfc2.tree_table -- highest scoring topology, arranged into a tabular format with stratigraphic ranges and branch support values

best_mfc2.tre -- highest scoring topology in newick format

best_mfc2.2nd.tree_table -- highest scoring topology, arranged into a tabular format with stratigraphic ranges and branch support values

best_mfc2.2nd.tre -- second highest scoring topology in newick format

best_mfc2.3rd.tree_table -- highest scoring topology, arranged into a tabular format with stratigraphic ranges and branch support values

best_mfc2.3.tre -- third highest scoring topology in newick format



## Commands

NOTE: The relative path for all of these commands assums that scripts will be run from the sub-folder in the `tests/micraster` directory in the stratoML repository (https://github.com/tomopfuku/stratoML). They can all be exchanged with the full path to the folder containing all modules in the repository.

### run tree search starting with stratophenetic tree
python3 ../../main_newick.py micrasteridae-stratophenetic.rr.tre micrasteridae-character_matrix.fa micrasteridae-stratigraphy_abs.stag.csv


### run tree search starting with parsimony tree
python3 ../../main_newick.py micrasteridae-parsimony.rr.txt micrasteridae-character_matrix.fa micrasteridae-stratigraphy_abs.stag.csv



### calculate support for best tree and plot as table

python plot_tree_table.py best_mfc2.tre micrasteridae-stratigraphy_abs.csv top_trees_mfc2 > best_mfc2.tree_table 



### calculate support for tree 2

python plot_tree_table.py best_mfc2.2nd.tre micrasteridae-stratigraphy_abs.csv top_trees_mfc2 > best_mfc2.2nd.tree_table 


### calculate support for tree 2


python plot_tree_table.py best_mfc2.3rd.tre micrasteridae-stratigraphy_abs.csv top_trees_mfc2 > best_mfc2.2nd.tree_table 


### calculate number of states displayed by each micraster taxon
python get_n_states.py micrasteridae-character_matrix.fa best_hr.tree_table > best_hr.NSTATES
