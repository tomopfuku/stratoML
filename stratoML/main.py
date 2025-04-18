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
        n.update_pmat(qmats,max(ss))
    times = []

    """
    stratll = stratlike.poisson_loglike(1.0,tree)
    print(stratll)

    res = minimize(stratlike.poisson_neg_ll,x0=np.array([1.0]),args=(tree),method="Powell")

    for _ in range(1):
        start = time.time()
        res_tr = minimize(mfc.evaluate_m_l,x0=np.array([0.1,0.0002]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.2),(0.00001,0.2)))
        res_st = minimize(stratlike.poisson_neg_ll,x0=np.array([0.1]),args=(tree),method="L-BFGS-B",bounds=[(0.00001,10)])
        stop = time.time()
        times.append(stop-start)

    print("median time per iteration: ", np.median(times), " seconds")

    tree_ll = -res_tr.fun + -res_st.fun
    print("mut. rate: ", res_tr.x[0])
    print("fix. rate: ", res_tr.x[1])
    print("pres. lambda: ", res_st.x[0])
    print("tree loglike: ",tree_ll, "morph. ll: ", -res_tr.fun, "strat. ll: ", -res_st.fun)
    #print(mfc.mfc_treell(tree,ss))
    """
    ## optimize bds model params:
    print(stratlike.bds_neg_ll(np.array([0.8,0.4,1.3]),tree))    

    #minimize(stratlike.bds_neg_ll,x0=np.array([0.5,0.5,1.0]),args=(tree),method="L-BFGS-B",bounds=((0.00001,10),(0.00001,10),(0.00001,10)))
    pqr_start = np.array([0.5,0.5,1.0])
    #res = minimize(stratlike.bds_neg_ll,x0=pqr_start,args=(tree),method="L-BFGS-B",bounds=((0.00001,10),(0.00001,10),(0.00001,10)))
    start = time.time()
    res = basinhopping(stratlike.bds_neg_ll,x0=pqr_start,niter=20,minimizer_kwargs={"method":"L-BFGS-B","args":(tree),"bounds":((0.00001,10),(0.00001,10),(0.00001,20))})
    stop = time.time()
    #res = dual_annealing(stratlike.bds_neg_ll,args=tree,visit=1,accept=1,bounds=((0.00001,10),(0.00001,10),(0.00001,10)),no_local_search=True)
    #print(res)
    #print(stop-start)

    ## optimize speciation times under BDS model:
    """
    sp_times = []
    for i,n in enumerate(tree.iternodes(0)):
        if n == tree and n.istip == False:
            continue
        n.index = i
        sp_times.append(n.lower+0.01)
    """
    p = res.x[0]
    q = res.x[1]
    r = res.x[2]
    ll = stratlike.bds_loglike(p,q,r,tree)
    print(ll)
    #stratlike.bds_dates(p,q,r,tree)
    #ll = stratlike.bds_loglike(p,q,r,tree)
    #print("trait + bds ll:",ll+-res_tr.fun)

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
    print("tree1 loglike: ",tree_ll, "morph. ll: ", -res_tr.fun, "bds ll: ", bdsll)
    print(tree.get_newick_repr())

    tree_utils.tree_search(tree,ss)
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

