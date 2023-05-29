import sys
import node
import tree_reader,read_fasta,tree_utils
import numpy as np
import mfc
import qmat
import time
from scipy.optimize import minimize
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
    retraits = read_fasta.recode_poly_traits(traits,ss)
    tree = tree_reader.read_tree_table(sys.argv[1])
    tree_utils.map_tree_disc_traits(tree,retraits,ss)
    qmats = qmat.Qmat(0.01,0.05)
    for n in tree.iternodes():
        #if n.istip:
        #    for i in n.disc_traits:
        #        print(list(i))
        n.update_pmat(qmats,max(ss))

    #for n in tree.iternodes():
    #    print([sum(list(i)) for i in list(n.pmats[0])])
    #sys.exit()

    treell = mfc.evaluate_m_l(np.array([0.006,0.04]),tree,qmats,ss) 
    print(treell)



    treell = mfc.evaluate_m_l(np.array([0.01,0.05]),tree,qmats,ss) 
    print(treell)

    treell = mfc.evaluate_m_l(np.array([6.996e-03,4.061e-02]),tree,qmats,ss) 
    print(treell)

    
    treell = mfc.evaluate_m_l(np.array([0.009,0.02]),tree,qmats,ss) 
    print(treell)


    times = []
    for _ in range(10):
        start = time.time()
        res = minimize(mfc.evaluate_m_l,x0=np.array([0.1,0.0002]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.2),(0.00001,0.2)))
        stop = time.time()
        times.append(stop-start)
    print("median time per iteration: ", np.median(times), " seconds")


    #res = minimize(mfc.evaluate_m_l,x0=[0.001,0.004],args=(tree,qmats,ss),method="COBYLA")
    print(res)

    print(mfc.mfc_treell(tree,ss))
    
    #for n in tree.iternodes():
    #    print([list(i) for i in list(n.pmats[0])])
    """  
    times = []
    for _ in range(1000):
        start = time.time()
        treell = mfc.mfc_treell(tree,ss)
        stop = time.time()
        times.append(stop-start)
    print("median time per iteration: ", np.median(times))
    """

