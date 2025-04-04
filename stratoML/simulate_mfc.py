import sys
import tree_reader,tree_utils
import node
import qmat
import buddmat as bm
import numpy as np
from random import choices
import smaps

def sim_along_branch(node, qmats, ss):
    node.update_pmat(qmats,ss)
    for i in range(len(node.disc_traits)):
        conv_pmats = []
        for j in node.pmats[ss-2]:
            conv_pmats.append(list(j))

        if node.parent == None:
            par_tr = [0.0] * ( 2 ** ss )
            par_tr[1] = 1.0
        else:
            par_tr = node.parent.disc_traits[i]
        
        par_state = [j for j in range(len(par_tr)) if par_tr[j] == 1.0][0]
        all_scen = bm.get_buddmat(ss)
        cur_scen = [j for j in all_scen[par_state] if j != 0]
        sim_scen = choices(cur_scen,k=1)[0]
        trans_probs = conv_pmats[sim_scen]
        trans_probs[0] = 0.0
        trans_probs = [j / sum(trans_probs) for j in trans_probs]
        states = [j for j in range(len(trans_probs))]
        sim_state = choices(states,k=1,weights=trans_probs)[0]
        node.disc_traits[i][sim_state] = 1.0
        
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
    for n in tree.iternodes(0):
        trait_probs = np.zeros((num_traits,2 ** ss),dtype=np.double)
        n.disc_traits = trait_probs
        sim_along_branch(n, qmats, ss)


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
    for row in qmats.get_qmat(2):
        print(list(row))

    for i in qmats.calc_p_mats(2.2,2)[0]:
        print(list(i))

    sim_traits_across_tree(tree,qmats,num_traits,ss)

    trait_dict = get_trait_dict(tree,ss)

    for i in trait_dict:
        print(">"+i)
        print(" ".join(trait_dict[i]))

