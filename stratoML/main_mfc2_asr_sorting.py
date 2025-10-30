import sys
import node
import tree_reader,read_fasta,tree_utils,stratlike,mfc
from main_mfc2_asr import *
import numpy as np
import qmat
import smaps
from scipy.optimize import minimize
from scipy.optimize import basinhopping
import time

def label_internal_nodes(tree):
    count = 0
    for n in tree.iternodes():
        if n.istip == False:
            lab = "n"+str(count)
            n.label = lab
            count += 1

def par_anc_state_polymorph(tree, ss):
    node_cis = stratlike.bds_CIs(p,q,r,tree)

    nl_cis = {}
    for n in node_cis:
        nl_cis[n.label] = node_cis[n]

    node_cis = nl_cis

    ncount = 0
    hemi_dict = {}
    codes = {}
    for n in tree.iternodes():
        hemi_prop = "NA" 
        if len(n.children) > 1:
            par_lv_obs = n.timeslice_lv[n.midpoint_lv_index]
            par_traits_obs = []
            for chari in range(1,len(par_lv_obs)):
                char = par_lv_obs[chari]
                max_index = mfc.index_of_max(char)
                par_traits_obs.append(max_index)
            par_traits_obs = recode_trait_vec(par_traits_obs, ss)

            desc_traits_all = []
            for ch in n.children:
                desc_traits = []
                ch_lv = ch.timeslice_lv[ch.midpoint_lv_index]
                for chari in range(1,len(par_lv_obs)):
                    char = ch_lv[chari]
                    max_index = mfc.index_of_max(char)
                    desc_traits.append(max_index)
                desc_traits = recode_trait_vec(desc_traits, ss)
                desc_traits_all.append(desc_traits)

            num_hemi = 0
            for i in range(len(par_traits_obs)):
                anc_tr = set(par_traits_obs[i].strip().split("|"))
                if len(anc_tr) < 2:
                    continue

                for j in range(len(desc_traits_all)):
                    ch_tr1 = set(desc_traits_all[j][i].strip().split("|"))
                    par_inter1 = ch_tr1.intersection(anc_tr)
                    if len(par_inter1) == 0:
                        continue

                    for k in range(len(desc_traits_all)):
                        if j >= k:
                            continue
                        ch_tr2 = set(desc_traits_all[k][i].strip().split("|"))
                        par_inter2 = ch_tr2.intersection(anc_tr)
                        if len(par_inter2) == 0:
                            continue

                        ch_intersect = ch_tr1.intersection(ch_tr2)
                        if len(ch_intersect) == 0:
                            num_hemi += 1
                            # TODO: do i need to break here to not double count hemiplasies when there are >2 desc?

            hemi_prop = round(num_hemi / float(len(par_traits_obs)) * 100)
        #n.label = "n"+str(ncount)
        codes[n.label] = str(ncount)
        ncount += 1
        hemi_dict[n.label] = str(hemi_prop)

    outfl = open(".".join(sys.argv[1].strip().split(".")[0:-1])+".tree_table.HEMI","w")
    outfl.write("Name,Code,Start,End,FAD,LAD,2.5HPD,97.5HPD,Parent,HemiPercent\n")
    for n in tree.inorder():
        curcode = codes[n.label]
        if n.parent != None:
            parcode = codes[n.parent.label]
        else:
            parcode = "NA"
        if n.istip:
            outfl.write(n.label+","+curcode+","+str(n.lower)+","+str(n.upper)+","+str(n.strat[0])+","+str(n.strat[1])+","+str(node_cis[n.label][0])+","+str(node_cis[n.label][1])+","+parcode+","+str(hemi_dict[n.label])+"\n")
        else:
            outfl.write(""+","+curcode+","+str(n.lower)+","+str(n.upper)+","+"NA"+","+"NA"+","+"NA"+","+"NA"+","+parcode+","+str(hemi_dict[n.label])+"\n")
    outfl.close()



if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("usage: "+ sys.argv[0]+ " <newick> <trait fasta file> <stratigraphic data> <stratigraphic model> <morphologic model>")
        sys.exit()
    
    traits,ss = read_fasta.read_fasta(sys.argv[2])
    retraits  = read_fasta.recode_poly_traits(traits,ss)

    for line in open(sys.argv[1],"r"):
        nwk = line.strip().split()[-1]
        tree = tree_reader.read_tree_string(nwk)
        tree_utils.map_strat_to_tree(tree,sys.argv[3])    
        tree_utils.map_tree_disc_traits(tree,retraits,ss)
        tree_utils.fix_obs_lv(tree) 
        stratlike.calibrate_brlens_strat(tree,0.4)
        pqr_start = np.array([0.5,0.5,1.0])
        res = basinhopping(stratlike.bds_neg_ll,x0=pqr_start,niter=20,minimizer_kwargs={"method":"L-BFGS-B","args":(tree),"bounds":((0.00001,10),(0.00001,10),(0.00001,20))})
        p = res.x[0]
        q = res.x[1]
        r = res.x[2]
        print("speciation:",p,"\nextinction:",q,"\npreservation",r,"\n")
        stratlike.bds_dates(p, q, r, tree)
        tree_utils.sort_children_by_age(tree)
        label_internal_nodes(tree)

        qmats = qmat.Qmat(0.02,0.01)
        res_tr = minimize(mfc.evaluate_m_l2,x0=np.array([0.01,0.01]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.5),(0.00001,0.5)))
        qmats = qmat.Qmat(res_tr.x[0],res_tr.x[1])
        mfc.compute_mfc2_ASRs(tree,qmats,ss) 
        #print(tree.timeslice_lv[0][1])
        par_anc_state_polymorph(tree, ss[1:])
