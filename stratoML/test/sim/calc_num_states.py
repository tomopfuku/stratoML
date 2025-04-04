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
    if len(sys.argv) != 2:
        print("usage: "+ sys.argv[0]+ " <simulated trait fastas>")
        sys.exit()

    pathlist = Path(sys.argv[1]).rglob('*.fa')
    corr_tree = []
    rates = []
    corr_tree_index = 0

    pref_trees = []
    n_char = {1:0,2:0,3:0}
    for it, fl in enumerate(pathlist):
        opfl = open(fl,"r")
        for line in opfl:
            if line[0] == ">":
                continue
            traits = line.strip().split()
            for tr in traits:
                s = tr.strip().split("|")
                n = len(s)
                n_char[n] += 1

    prop = [i / sum(n_char.values()) for i in n_char.values()]
    
    for i in n_char:
        n_char[i] = prop[i-1]
    print(n_char)

    plt.bar(list(n_char.keys()),list(n_char.values()),color="grey")
    plt.xticks(list(n_char.keys()), list(n_char.keys()))
    plt.show()

    fig, ax = plt.subplots()
    colors = ['#ff9999','#99ff99','#ffcc99']
    ax.pie(list(n_char.values()), labels=n_char.keys(), autopct='%1.1f%%',colors=colors)
    plt.show()
    """
    width = 0.5
    fig, ax = plt.subplots()
    bottom = 0.

    for n_st, prop in n_char.items():
        p = ax.bar(0,prop,bottom=bottom)
        bottom+=prop
    plt.show()
    """
