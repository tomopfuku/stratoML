import sys
import node
import tree_reader,read_fasta,tree_utils,stratlike,mfc
import numpy as np
import qmat
import smaps
from scipy.optimize import minimize
import time

def get_trait_dict(tree,ss):
    traits = {}
    for n in tree.iternodes():
        traits[n.label] = []
        cur_tr = n.disc_traits
        re_tr = []
        for i in cur_tr:
            cur_state = [j for j in range(len(list(i))) if i[j] == 1.0][0]
            re_tr.append(cur_state)
        #print(n.label,re_tr)
        re_re_tr = recode_trait_vec(re_tr,ss)
        traits[n.label] = re_re_tr
    return traits

def recode_trait_vec(trait_vec, ss):
    newtr = []
    for ind, i in enumerate(trait_vec):
        if ss[ind] == 1:
            newtr.append("0")
            continue
        
        smap = smaps.get_smap(ss[ind])
        pres_row = smap[i]
        new_coding = []
        for j,k in enumerate(pres_row):
            if k != 0:
                new_coding.append(str(j))
        new_coding = "|".join(new_coding)
        newtr.append(new_coding)
    return newtr

def pick_ML_anc_state(tree, ss):
    for n in tree.iternodes():
        for lv in n.timeslice_lv:
            for chari in range(1,len(lv)):
                char = lv[chari]
                max_index = mfc.index_of_max(char)
                for i in range(len(lv[chari])):
                    lv[chari][i] = 0.0
                lv[chari][max_index] = 1.0



def par_anc_state_polymorph(tree, ss):
    for n in tree.iternodes():
        if n == tree:
            continue

        par_lv = n.parent.timeslice_lv[n.parent_lv_index]
        par_lv_obs = n.parent.timeslice_lv[n.parent.midpoint_lv_index]
        desc_lv = n.timeslice_lv[n.midpoint_lv_index]
        par_traits = []
        par_traits_obs = []
        desc_traits = []
        for chari in range(1,len(par_lv)):
            char = par_lv[chari]
            max_index = mfc.index_of_max(char)
            par_traits.append(max_index)

            char = par_lv_obs[chari]
            max_index = mfc.index_of_max(char)
            par_traits_obs.append(max_index)

            char = desc_lv[chari]
            max_index = mfc.index_of_max(char)
            desc_traits.append(max_index)
        par_traits = recode_trait_vec(par_traits, ss)
        par_traits_obs = recode_trait_vec(par_traits_obs, ss)
        desc_traits = recode_trait_vec(desc_traits, ss)
        n_polymorph_par = get_total_char(par_traits)
        n_polymorph_par_obs = get_total_char(par_traits_obs)
        n_polymorph_desc = get_total_char(desc_traits)
        if n.label != "":
            print(n.label,n_polymorph_par,n_polymorph_par_obs,n_polymorph_desc)

def get_total_char(trait_vec):
    count = 0
    for i in trait_vec:
        spls = i.strip().split("|")
        count += len(spls)
    return count

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
        #stratlike.calibrate_brlens_strat(tree,0.3)
        tree_utils.map_tree_disc_traits(tree,retraits,ss)
        tree_utils.fix_obs_lv(tree) 
        #tree_utils.sort_children_by_age(tree)
        #tree_utils.init_budd_marginals(tree,len(ss))

        qmats = qmat.Qmat(0.02,0.01)
        res_tr = minimize(mfc.evaluate_m_l2,x0=np.array([0.01,0.01]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.5),(0.00001,0.5)))
        print("mfc rates",res_tr.x)
        qmats = qmat.Qmat(res_tr.x[0],res_tr.x[1])
        mfc.compute_mfc2_ASRs(tree,qmats,ss) 
        #print(tree.timeslice_lv[0][1])
        par_anc_state_polymorph(tree, ss[1:])
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
        #aic = tree_utils.single_tree_aic(tree,ss,sys.argv[4],sys.argv[5])
        #print(aic,tree.get_newick_repr()+";")
        """
