import sys
import node
import tree_reader,read_fasta,tree_utils
import numpy as np
import mfc

"""test = node.Node()
print(test)
test.add_disc_traits([[0,1],[0],[0,1,2]])
print(len(test.disc_traits))
print(list(test.disc_traits[0]))
print(list(test.disc_traits[1]))"""

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: "+ sys.argv[0]+ " <tree table file> <trait fasta file>")
        sys.exit()

    traits,ss = read_fasta.read_fasta(sys.argv[2])
    #print(traits)
    print(ss)
    tree = tree_reader.read_tree_table(sys.argv[1])
    print(tree)
    tree_utils.map_tree_disc_traits(tree,traits)
    qmat,smap = mfc.create_q_matrix(3, 0.1, 0.3)
    pmat = mfc.calc_p_matrix(qmat,.5)
    for i in pmat:
        print(list(i))
