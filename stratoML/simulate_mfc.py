import sys
import tree_reader,tree_utils
import node
import qmat
import buddmat as bm
import numpy as np
from random import choices
from random import choice
import smaps


def print_memoryslice_2d(ms):
    for i in ms:
        print(list(i))

# this version implements a proper anagenetic model, with changes computed at each budding point
def sim_along_branch2(node, qmats, ss):
    for i in range(len(node.disc_traits)):
        if node.parent == None:
            par_tr = [0.0] * ( 2 ** ss )
            par_st = choice([i for i in range(1, 2 ** ss)])
            par_tr[par_st] = 1.0
        else:
            par_tr = node.parent.budd_marginals[node.index_from_parent][i]

        par_state = [j for j in range(len(par_tr)) if par_tr[j] == 1.0][0]
        all_scen = bm.get_buddmat(ss)
        cur_scen = [j for j in all_scen[par_state] if j != 0]
        sim_scen = choices(cur_scen,k=1)[0]
        midpoint = node.lower - ( (node.lower - node.upper) / 2. )
        timeslices = [node.lower] + [ch.lower for ch in node.children]
        midpoint_index = 0
        for j in timeslices:
            if j > midpoint:
                midpoint_index+=1

        timeslices.append(midpoint)#+=[node.lower,midpoint]
        timeslices.sort(reverse=True)
        timeseries = [sim_scen]

        for j in range(1,len(timeslices)):
            t = timeslices[j]
            dt = timeslices[j-1] - t 
            pmat = qmats.calc_single_p_mat(dt,ss)
            last_state = timeseries[j-1]
            trans_probs = pmat[last_state]
            trans_probs[0] = 0.0
            trans_probs = [k / sum(trans_probs) for k in trans_probs]
            states = [k for k in range(len(trans_probs))]
            sim_state = choices(states,k=1,weights=trans_probs)[0]
            timeseries.append(sim_state)
            cur_traits = [0.0] * ( 2 ** ss )
            cur_traits[sim_state] = 1.0

            if j < midpoint_index:
                node.budd_marginals[j-1][i]=np.array(cur_traits)
            elif j > midpoint_index:
                node.budd_marginals[j-2][i]=np.array(cur_traits)
            else:
                node.disc_traits[i] = np.array(cur_traits)
  
   
       
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

def recode_trait_vec(trait_vec,ss):
    smap = smaps.get_smap(ss)
    newtr = []
    for i in trait_vec:
        pres_row = smap[i]
        new_coding = []
        for j,k in enumerate(pres_row):
            if k != 0:
                new_coding.append(str(j))
        new_coding = "|".join(new_coding)
        newtr.append(new_coding)
    return newtr

def sim_traits_across_tree(tree,qmats,num_traits,ss):
    tree_utils.sort_children_by_age(tree) 
    for n in tree.iternodes(0):
        trait_probs = np.zeros((num_traits,2 ** ss),dtype=np.double)
        n.disc_traits = trait_probs
        #print_memoryslice_2d(n.disc_traits)
        #sys.exit()
        #print(n.label,len(n.children)) 
        if len(n.children) > 0:
            marginals = []
            for i in range(len(n.children)):
                marginals.append(trait_probs)

            n.budd_marginals = np.array(marginals)
        sim_along_branch2(n, qmats, ss)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("usage:" + sys.argv[0] + " <tree> <stratigraphic data> <num_traits> <max_num_states>")
        sys.exit()

    nwk = open(sys.argv[1],"r").readline()
    tree = tree_reader.read_tree_string(nwk)
    tree_utils.map_strat_to_tree(tree,sys.argv[2]) 
    num_traits = int(sys.argv[3])
    ss = int(sys.argv[4])

    qmats = qmat.Qmat(0.08,0.08)
    #for row in qmats.get_qmat(2):
    #    print(list(row))

    #for i in qmats.calc_p_mats(2.2,2)[0]:
    #    print(list(i))

    stratlike.calibrate_brlens_strat(tree,0.3)
    tree_utils.sort_children_by_age(tree) 
    sys.exit()
    sim_traits_across_tree(tree,qmats,num_traits,ss)

    trait_dict = get_trait_dict(tree,ss)

    for i in trait_dict:
        print(">"+i)
        print(" ".join(trait_dict[i]))

