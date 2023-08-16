import sys
import node
import tree_reader,read_fasta,tree_utils
import numpy as np
import mfc
import qmat
import stratlike
import time
from scipy.optimize import minimize
from scipy.optimize import dual_annealing
from scipy.optimize import basinhopping 
import bd
"""test = node.Node()
print(test)
test.add_disc_traits([[0,1],[0],[0,1,2]])
print(len(test.disc_traits))
print(list(test.disc_traits[0]))
print(list(test.disc_traits[1]))"""


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: "+ sys.argv[0]+ " <newick> <trait fasta file> <stratigraphic data>")
        sys.exit()

    traits,ss = read_fasta.read_fasta(sys.argv[2])
    #print(traits)
    retraits = read_fasta.recode_poly_traits(traits,ss)
    nwk = open(sys.argv[1],"r").readline().strip()
    tree = tree_reader.read_tree_string(nwk)
    tree_utils.map_strat_to_tree(tree,sys.argv[3])    
    tree_utils.map_tree_disc_traits(tree,retraits,ss)
    qmats = qmat.Qmat(0.01,0.05)
    for n in tree.iternodes():
        n.update_pmat(qmats,max(ss))


    print(tree.get_newick_repr())
    tree_utils.tree_search2(tree,ss)

    print(tree.get_newick_repr())
    #changed = tree_utils.find_best_spr(tree,ss)
    #print(changed)
    #print(tree.get_newick_repr())
    #tree_utils.search_ancestors(tree,ss)
    #tree_utils.search_bifurcating(tree,ss)
    
    """
    ## optimize bds model params:
    pqr_start = np.array([0.5,0.5,1.0])
    res = basinhopping(stratlike.bds_neg_ll,x0=pqr_start,niter=20,minimizer_kwargs={"method":"L-BFGS-B","args":(tree),"bounds":((0.00001,10),(0.00001,10),(0.00001,20))})
    p = res.x[0]
    q = res.x[1]
    r = res.x[2]
    ll = stratlike.bds_loglike(p,q,r,tree)

    for _ in range(1):
        stratlike.bds_dates(p,q,r,tree)
        bdsll = stratlike.bds_loglike(p,q,r,tree)
        for n in tree.iternodes():
            n.update_pmat(qmats,max(ss))
        res_tr = minimize(mfc.evaluate_m_l,x0=np.array([0.1,0.0002]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.2),(0.00001,0.2)))

    tree_ll = bdsll+-res_tr.fun 
    print("mut. rate: ", res_tr.x[0])
    print("fix. rate: ", res_tr.x[1])
    #print("pres. lambda: ", res_st.x[0])
    print("tree1 loglike: ",tree_ll, "morph. ll: ", -res_tr.fun, "bds ll: ", bdsll)

    print(tree.get_newick_repr())

    pluck_node, prev_par, sib = tree_utils.random_spr(tree)

    stratlike.calibrate_brlens_strat(tree)
    qmats = qmat.Qmat(0.01,0.05)
    for n in tree.iternodes():
        n.update_pmat(qmats,max(ss))

    res = basinhopping(stratlike.bds_neg_ll,x0=pqr_start,niter=20,minimizer_kwargs={"method":"L-BFGS-B","args":(tree),"bounds":((0.00001,10),(0.00001,10),(0.00001,20))})
    p = res.x[0]
    q = res.x[1]
    r = res.x[2]
    ll = stratlike.bds_loglike(p,q,r,tree)
    print(ll)
    for _ in range(1):
        start = time.time()
        stratlike.bds_dates(p,q,r,tree)
        bdsll = stratlike.bds_loglike(p,q,r,tree)
        for n in tree.iternodes():
            n.update_pmat(qmats,max(ss))
        res_tr = minimize(mfc.evaluate_m_l,x0=np.array([0.1,0.0002]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.2),(0.00001,0.2)))
        stop = time.time()
        times.append(stop-start)

    tree_ll = bdsll+-res_tr.fun 
    print("mut. rate: ", res_tr.x[0])
    print("fix. rate: ", res_tr.x[1])
    #print("pres. lambda: ", res_st.x[0])
    print("tree2 loglike: ",tree_ll, "morph. ll: ", -res_tr.fun, "bds ll: ", bdsll)
    print(tree.get_newick_repr())

    tree_utils.reattach_to_orig_parent(pluck_node,prev_par,sib)
    print(tree.get_newick_repr())

    stratlike.calibrate_brlens_strat(tree)
    qmats = qmat.Qmat(0.01,0.05)
    for n in tree.iternodes():
        n.update_pmat(qmats,max(ss))

    res = basinhopping(stratlike.bds_neg_ll,x0=pqr_start,niter=20,minimizer_kwargs={"method":"L-BFGS-B","args":(tree),"bounds":((0.00001,10),(0.00001,10),(0.00001,20))})
    p = res.x[0]
    q = res.x[1]
    r = res.x[2]
    for _ in range(1):
        stratlike.bds_dates(p,q,r,tree)
        bdsll = stratlike.bds_loglike(p,q,r,tree)
        for n in tree.iternodes():
            n.update_pmat(qmats,max(ss))
        res_tr = minimize(mfc.evaluate_m_l,x0=np.array([0.1,0.0002]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.2),(0.00001,0.2)))

    tree_ll = bdsll+-res_tr.fun 
    print("mut. rate: ", res_tr.x[0])
    print("fix. rate: ", res_tr.x[1])
    #print("pres. lambda: ", res_st.x[0])
    print("tree1 loglike: ",tree_ll, "morph. ll: ", -res_tr.fun, "bds ll: ", bdsll)
    print(tree.get_newick_repr())
    """


    #print(tree.get_newick_repr())
    #res = stratlike.bds_node_starts_neg_ll(np.array(sp_times),p,q,r,tree)
    #res = basinhopping(stratlike.bds_node_starts_neg_ll,x0=np.array(sp_times),niter=100,minimizer_kwargs={"method":"Nelder-Mead","args":(p,q,r,tree)})
    #print(res)
    """ for n in tree.iternodes(0):
        if n.istip:
            print("\n\nNEW",n.strat[0]-n.strat[1],p,q,r)
            stratlike.single_range_mle(p,q,r,n)
        else:
            print("\n\nHYPOTHETICAL ANC")
            stratlike.hyp_anc_mle(p,q,0.3)
        #print(n.label,n.lower,n.strat[0])"""

    """  
    times = []
    for _ in range(1000):
        start = time.time()
        treell = mfc.mfc_treell(tree,ss)
        stop = time.time()
        times.append(stop-start)
    print("median time per iteration: ", np.median(times))
    """

