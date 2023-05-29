import sys
import numpy as np
#import read_fasta, tree_reader

def map_tree_disc_traits(tree,traits,ss):
    for n in tree.iternodes(1):
        if n.istip:
            try:
                curtr = traits[n.label]
                n.add_disc_traits(curtr,ss)
                #for i in n.disc_traits:
                #    print(list(i))
            except:
                print(n.label, "is present in the tree, but has no traits in the fasta")
                sys.exit()
            #print(curtr)
            #print(list(n.disc_starts))
            #print(list(n.disc_states))
        else:
            n.disc_traits = np.zeros((len(curtr),128))

    