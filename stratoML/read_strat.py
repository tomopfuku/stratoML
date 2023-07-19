import numpy as np
import sys
    
def map_strat_to_tree(tree, flnm):
    fl = open(flnm, "r")
    fl.readline()
    ranges = {}
        
    for line in fl:
        spls = line.strip().split()
        spnm = spls[0]
        fad = float(spls[1])
        lad = float(spls[2])
        ranges[spnm] = np.array([fad,lad])

    for node in tree.iternodes():
        if node.istip:
            try:
                strat = ranges[n.label]
                n.strat = strat
            except:
                print(n.label," is in the tree but was not found in the stratigraphic range data")
                sys.exit()

 