import argparse
import sys
import node
import tree_reader,read_fasta,tree_utils,stratlike,mfc
import numpy as np
import qmat
from scipy.optimize import minimize
import time

def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", required=True, help="Input Newick tree file.")
    parser.add_argument("--traits", required=True, help="Trait FASTA file.")
    parser.add_argument("--strat-data", required=True, help="Stratigraphic data file.")
    parser.add_argument("--strat-model", required=True, help="Stratigraphic model.")
    parser.add_argument("--morph-model", required=True, help="Morphologic model.")
    return parser.parse_args(argv)

if __name__ == "__main__":
    args = parse_args()
    
    traits, ss = read_fasta.read_fasta(args.traits)
    retraits  = read_fasta.recode_poly_traits(traits,ss)

   
    for line in open(args.trees,"r"):
        nwk = line.strip().split()[-1]
        tree = tree_reader.read_tree_string(nwk)

        tree_utils.map_strat_to_tree(tree,args.strat_data)    
        #stratlike.calibrate_brlens_strat(tree,0.3)
        tree_utils.map_tree_disc_traits(tree,retraits,ss)
        tree_utils.fix_obs_lv(tree) 
        #tree_utils.sort_children_by_age(tree)
        #tree_utils.init_budd_marginals(tree,len(ss))
        qmats = qmat.Qmat(0.01,0.05)

        for n in tree.iternodes():
            n.update_pmat(qmats,max(ss),"mid")


        #pmat = qmats.calc_single_p_mat(1.0, 2)
        #for i in pmat:
        #    print(list(i))

        #print(tree)
        #treell = -mfc.evaluate_m_l2(np.array([0.01,0.05]),tree,qmats,ss)
        #print("TREELL 1", treell) 
        #treell = -mfc.evaluate_m_l2(np.array([0.01,0.001]),tree,qmats,ss)
        #print("TREELL 2", treell)
        t1 = time.time()
        aic,traitll,bdsll = tree_utils.calc_tree_ll2(tree,qmats,ss,"hr97")
        t2 = time.time()
        #print(aic,traitll,bdsll)
        print(aic,nwk)
        #print("TIME OPTIMIZING",t2-t1)
        #tree_utils.tree_search4(tree,ss,qmats,"hr97",False)
        tree_utils.tree_search3(tree,ss,qmats,"hr97",False)

        """
        res_st = minimize(stratlike.poisson_neg_ll,x0=np.array([1.0]),args=(tree),method="Nelder-Mead")
        bdsll = -res_st.fun
        print("stratlike:",bdsll)
        nparam = 1.0 + 2.0
        res_tr = minimize(mfc.evaluate_m_l2,x0=np.array([0.01,0.01]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.5),(0.00001,0.5)))
        print("mfc rates",res_tr.x)
        traitll = -res_tr.fun
        print(traitll, len([n for n in tree.iternodes()]))
        tree_ll = traitll + bdsll

        nparam += float(len([n for n in tree.iternodes()]) - 1)
        aic = (2. * nparam) - (2. * tree_ll) 
        print("AIC",aic)
        #aic = tree_utils.single_tree_aic(tree,ss,args.strat_model,args.morph_model)
        #print(aic,tree.get_newick_repr()+";")
        """
