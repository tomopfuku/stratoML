import sys
import node
import tree_reader,read_fasta,tree_utils
import numpy as np
import qmat

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("usage: "+ sys.argv[0]+ " <newick> <trait fasta file> <stratigraphic data> <stratigraphic model>")
        sys.exit()

    traits,ss = read_fasta.read_fasta(sys.argv[2])
    retraits  = read_fasta.recode_poly_traits(traits,ss)
    for line in open(sys.argv[1],"r"):
        nwk = line.strip().split()[-1]
        tree = tree_reader.read_tree_string(nwk)
        tree_utils.map_strat_to_tree(tree,sys.argv[3])    
        tree_utils.map_tree_disc_traits(tree,retraits,ss)
        
        qmats = qmat.Qmat(0.01,0.05)
        for n in tree.iternodes():
            n.update_pmat(qmats,max(ss))

        
        aic = tree_utils.single_tree_aic(tree,ss,sys.argv[4])
        print(aic,tree.get_newick_repr()+";")

