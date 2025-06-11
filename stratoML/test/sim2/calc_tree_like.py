import sys
sys.path.append("../../")
import node
import tree_reader,read_fasta,tree_utils
import numpy as np
import qmat
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("usage: "+ sys.argv[0]+ " <newick> <simulated trait fastas> <stratigraphic data> <stratigraphic model>")
        sys.exit()

    nwks = open(sys.argv[1],"r").readlines()
    num_trees = len(nwks)

    pathlist = Path(sys.argv[2]).rglob('*.fa')
    corr_tree = []
    rates = []
    corr_tree_index = 0

    pref_trees = []
    for it, fl in enumerate(pathlist):
        #if it == 11:
        #    break
        flstr = str(fl).strip().split("/")[1]
        pref = flstr.strip().replace(".fa","")
        spls = pref.strip().split("_")
        rate = float(spls[1])
        rates.append(rate)
        traits,ss = read_fasta.read_fasta(fl)
        retraits  = read_fasta.recode_poly_traits(traits,ss)
        bestAIC = 10000000000.
        bestIND = None
        for i, line in enumerate(nwks):
            nwk = line.strip().split()[-1]
            tree = tree_reader.read_tree_string(nwk)
            tree_utils.map_strat_to_tree(tree,sys.argv[3])    
            tree_utils.map_tree_disc_traits(tree,retraits,ss)
            tree_utils.fix_obs_lv(tree) 
            qmats = qmat.Qmat(0.05,0.01)
            for n in tree.iternodes():
                n.update_pmat(qmats,max(ss))


            aic,traitll,bdsll = tree_utils.calc_tree_ll2(tree,qmats,ss,"hr97")
            #aic = tree_utils.single_tree_aic2(tree,ss,sys.argv[4])
            if aic < bestAIC:
                bestAIC = aic
                bestIND = i
        corr_tree.append( bestIND == corr_tree_index )
        print(fl,bestIND == corr_tree_index, bestIND, rate)
        pref_trees.append(bestIND)

    rates_trees = pd.DataFrame({"corr_tree":corr_tree,"rate":rates})
    print(rates_trees) 
    sns.boxplot(data=rates_trees, x = "corr_tree", y ="rate",color="grey").set(xlabel="Correct Tree")
    plt.show()

    counts = {}
    for i in range(num_trees):
        n = pref_trees.count(i)
        prop = n / float(len(pref_trees))
        counts[i] = prop

    print(counts)

    plt.bar(range(len(counts)), list(counts.values()), align='center', color = "grey")
    plt.xticks(range(len(counts)), list(counts.keys()))
    plt.show()
