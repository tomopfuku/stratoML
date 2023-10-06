import sys
import node
import tree_reader,read_fasta,tree_utils
import numpy as np
import qmat

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: "+ sys.argv[0]+ " <newick> <trait fasta file> <stratigraphic data>")
        sys.exit()

    traits,ss = read_fasta.read_fasta(sys.argv[2])
    retraits  = read_fasta.recode_poly_traits(traits,ss)
    nwk = open(sys.argv[1],"r").readline().strip()
    tree = tree_reader.read_tree_string(nwk)
    tree_utils.map_strat_to_tree(tree,sys.argv[3])    
    tree_utils.map_tree_disc_traits(tree,retraits,ss)
    
    qmats = qmat.Qmat(0.01,0.05)
    for n in tree.iternodes():
        n.update_pmat(qmats,max(ss))

    print(tree.get_newick_repr())
    tree_utils.tree_search2(tree,ss,"bds")

    print(tree.get_newick_repr())
