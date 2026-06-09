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

def par_anc_state_polymorph_DEP(tree, ss):
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

def par_anc_state_polymorph(tree, ss):
    node_cis = stratlike.bds_CIs(p,q,r,tree)

    nl_cis = {}
    for n in node_cis:
        nl_cis[n.label] = node_cis[n]

    node_cis = nl_cis

    ncount = 0
    hemi_dict = {}
    unsort_dict = {}
    uninf_dict  = {}
    codes = {}
    traits = {}
    for n in tree.iternodes():
        hemi_prop = "NA" 
        nhemi = 0
        nunsort = 0
        nuninf = 0
        par_lv_obs = n.timeslice_lv[n.midpoint_lv_index]
        par_traits_obs = []
        for chari in range(1,len(par_lv_obs)):
            char = par_lv_obs[chari]
            max_index = mfc.index_of_max(char)
            par_traits_obs.append(max_index)
        par_traits_obs = recode_trait_vec(par_traits_obs, ss)
        traits[n] = par_traits_obs

        if len(n.children) > 1:
            desc_traits_all = []
            for ch in n.children:
                desc_traits_all.append([])
                for desc in ch.leaves():
                    desc_traits = []
                    desc_lv = desc.timeslice_lv[desc.midpoint_lv_index]
                    for chari in range(1,len(par_lv_obs)):
                        char = desc_lv[chari]
                        max_index = mfc.index_of_max(char)
                        desc_traits.append(max_index)
                    desc_traits = recode_trait_vec(desc_traits, ss)
                    desc_traits_all[-1].append(desc_traits)

            for i in range(len(par_traits_obs)):
                anc_tr = set(par_traits_obs[i].strip().split("|"))
                all_subtrees = []
                for subtree in desc_traits_all:
                    all_subtrees.append([])
                    for desc_tr in subtree:
                        curtr = desc_tr[i]
                        splttr = set(curtr.strip().split("|"))
                        from_anc = sorted(list(splttr.intersection(anc_tr)))
                        if len(from_anc) > 0:
                            all_subtrees[-1].append("|".join(from_anc))
 
                if len(anc_tr) < 2:
                    uninf = find_uninf(anc_tr, all_subtrees)
                    if uninf:
                        nuninf += 1
                    continue
               
                discord = find_incomp_sort(all_subtrees)
                if discord == True:
                    nhemi += 1
                else:
                    unsorted = find_unsort(all_subtrees)
                    if unsorted:
                        nunsort += 1

        denom = float(len(par_traits_obs) - nuninf )
        if denom == 0:
            denom = 1
        hemi_prop = round((nhemi / denom) * 100)
        unsort_prop = round((nunsort / denom) * 100)
        print(n.label, hemi_prop, unsort_prop)
        #n.label = "n"+str(ncount)
        codes[n.label] = str(ncount)
        ncount += 1
        hemi_dict[n.label] = str(hemi_prop)
        unsort_dict[n.label] = str(unsort_prop)

    outfl = open(".".join(sys.argv[1].strip().split(".")[0:-1])+".tree_table.HEMI","w")
    outfl.write("Name,Code,Start,End,FAD,LAD,2.5HPD,97.5HPD,Parent,HemiPercent,UnsortPercent\n")
    for n in tree.inorder():
        curcode = codes[n.label]
        if n.parent != None:
            parcode = codes[n.parent.label]
            if n.istip and len(n.children) > 1:
                outfl.write(n.label+","+curcode+","+str(n.lower)+","+str(n.upper)+","+str(n.strat[0])+","+str(n.strat[1])+","+str(node_cis[n.label][0])+","+str(node_cis[n.label][1])+","+parcode+","+str(hemi_dict[n.label])+","+str(unsort_dict[n.label])+"\n")
            elif n.istip and len(n.children) < 2:
                outfl.write(n.label+","+curcode+","+str(n.lower)+","+str(n.upper)+","+str(n.strat[0])+","+str(n.strat[1])+","+str(node_cis[n.label][0])+","+str(node_cis[n.label][1])+","+parcode+","+"NA"+","+"NA"+"\n")
            else:
                outfl.write(""+","+curcode+","+str(n.lower)+","+str(n.upper)+","+"NA"+","+"NA"+","+"NA"+","+"NA"+","+parcode+","+str(hemi_dict[n.label])+","+str(unsort_dict[n.label])+"\n")

        else:
            parcode = "NA"
            outfl.write(","+curcode+","+str(n.lower)+","+str(n.upper)+","+str(n.strat[0])+","+str(n.strat[1])+","+str(node_cis[n.label][0])+","+str(node_cis[n.label][1])+","+parcode+","+"NA"+","+"NA"+"\n")
            continue
    outfl.close()

    """
    # uncomment to look at ancst recons
    for n in tree.inorder():
        curcode = codes[n.label]
        c = " ".join(traits[n])
        if n.istip:
            label = n.label
        else:
            label = "n"+str(curcode)
        print(">"+label)
        print(c)
    """

def find_uninf(par_tr, all_subtrees):
    discord = False
    states = []
    for j in range(len(all_subtrees)):
        for k in range(len(all_subtrees)):
            if j >= k:
                continue
            for st in all_subtrees[j]:
                states.append(st)
            for st in all_subtrees[k]:
                states.append(st)
    if len(list(set(states))) == 1:
        discord = True
    return discord


def find_unsort(all_subtrees):
    discord = False
    for j in range(len(all_subtrees)):
        for k in range(len(all_subtrees)):
            if j >= k:
                continue
            num_poly = 0
            st1 = []
            for st in all_subtrees[j]:
                spls = st.strip().split("|")
                if len(spls) > 1:
                    num_poly += 1
                st1 += spls
            st2 = []
            for st in all_subtrees[k]:
                spls = st.strip().split("|")
                if len(spls) > 1:
                    num_poly += 1
                st2 += spls

            st1 = set(st1)
            st2 = set(st2)
            shared = st1.intersection(st2)
            #print(shared)
            #shared = [char for char in shared if len(char.split("|"))>1]
            #allst = st1.union(st2)
            #allst = [char for char in allst if len(char.split("|"))>1]
            #print(allst)
            #unshared = set(allst) - set(shared)
            if len(shared) > 0 and num_poly > 0:
                discord = True
                break
    return discord


def find_incomp_sort(all_subtrees):
    discord = False
    for j in range(len(all_subtrees)):
        for k in range(len(all_subtrees)):
            if j >= k:
                continue
            st1 = set(all_subtrees[j])
            st2 = set(all_subtrees[k])
            shared = st1.intersection(st2)
            shared = [char for char in shared if len(char.split("|"))==1]
            #print(shared)
            allst = st1.union(st2)
            allst = [char for char in allst if len(char.split("|"))==1]
            unshared = set(allst) - set(shared)
            #print(unshared)
            if len(shared) > 0 and len(unshared) > 0:
                discord = True
                break
    return discord

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
