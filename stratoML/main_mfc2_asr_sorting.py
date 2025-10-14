import sys
import node
import tree_reader,read_fasta,tree_utils,stratlike,mfc
from main_mfc2_asr import *
import numpy as np
import qmat
import smaps
from scipy.optimize import minimize
import time

def par_anc_state_polymorph(tree, ss):
    for n in tree.iternodes():
        
        par_lv = n.timeslice_lv[n.parent_lv_index]
        par_lv_obs = n.parent.timeslice_lv[n.parent.midpoint_lv_index]
        desc_lv = n.timeslice_lv[n.midpoint_lv_index]
        par_traits = []
        par_traits_obs = []
        desc_traits = []
        for chari in range(1,len(par_lv)):
            char = par_lv[chari]
            max_index = mfc.index_of_max(char)
            par_traits.append(max_index)

        par_traits = recode_trait_vec(par_traits, ss)
        n_polymorph_par = get_total_char(par_traits)



if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("usage: "+ sys.argv[0]+ " <newick> <trait fasta file> <stratigraphic data> <stratigraphic model> <morphologic model>")
        sys.exit()
    
    traits,ss = read_fasta.read_fasta(sys.argv[2])
    retraits  = read_fasta.recode_poly_traits(traits,ss)

    print("taxon nstates_par_obs nstates_par_obs nstates")   
    for line in open(sys.argv[1],"r"):
        nwk = line.strip().split()[-1]
        tree = tree_reader.read_tree_string(nwk)
        tree_utils.map_strat_to_tree(tree,sys.argv[3])    
        tree_utils.map_tree_disc_traits(tree,retraits,ss)
        tree_utils.fix_obs_lv(tree) 

        qmats = qmat.Qmat(0.02,0.01)
        res_tr = minimize(mfc.evaluate_m_l2,x0=np.array([0.01,0.01]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.5),(0.00001,0.5)))
        print("mfc rates",res_tr.x)
        qmats = qmat.Qmat(res_tr.x[0],res_tr.x[1])
        mfc.compute_mfc2_ASRs(tree,qmats,ss) 
        #print(tree.timeslice_lv[0][1])
        par_anc_state_polymorph(tree, ss[1:])
