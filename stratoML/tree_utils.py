import sys
#import read_fasta, tree_reader

def map_tree_disc_traits(tree,traits):
    for n in tree.iternodes():
        if n.istip:
            try:
                curtr = traits[n.label]
                n.add_disc_traits(curtr)
            except:
                print(n.label, "is present in the tree, but has no traits in the fasta")
                sys.exit()
            #print(curtr)
            #print(list(n.disc_starts))
            #print(list(n.disc_states))

    