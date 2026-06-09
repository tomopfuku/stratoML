import matplotlib.pyplot as plt
import numpy as np
from simulate_mfc import *

## this simulates along a "split" branch assuming a budding model to match the treatment under GLC_BD model
def sim_split_branch2(node, qmats, ss, trait_ind = None, random_start = True):
    if node.istip:
        print("ERROR: trying to simulate bifurcation along sampled ancestor")
        sys.exit()

    if trait_ind == None:
        trait_ind = range(len(node.disc_traits))

    for i in trait_ind:
        if node.parent == None:
            par_tr = [0.0] * ( 2 ** ss )
            if random_start == False:
                par_st = 1 #choice([i for i in range(1, 2 ** ss)])
            else:
                par_st = choice([i for i in range(1, 2 ** ss)])
            par_tr[par_st] = 1.0
            node.length = 0.0
            #node.disc_traits[i] = np.array(par_tr)
            #node.budd_marginals[0][i] = np.array(par_tr)
        else:
            par_tr = node.parent.budd_marginals[node.index_from_parent][i]

        par_state = [j for j in range(len(par_tr)) if par_tr[j] == 1.0][0]
        dt = node.length
        pmat = qmats.calc_single_p_mat(dt, ss)
        trans_probs = pmat[par_state]
        trans_probs[0] = 0.0
        trans_probs = [k / sum(trans_probs) for k in trans_probs]
        states = [k for k in range(len(trans_probs))]
        sim_state = choices(states,k=1,weights=trans_probs)[0]
        cur_traits = [0.0] * ( 2 ** ss )
        cur_traits[sim_state] = 1.0
        node.disc_traits[i] = np.array(cur_traits)
        
        all_scen = sm.get_spltmat(ss)
        cur_scen = [j for j in all_scen[sim_state] if j[0] != 0]
        #print(sim_state)
        #print(cur_scen)
        weights = np.zeros(len(cur_scen))
        l_sub = SUB_RATE 
        l_sym = 1 - l_sub
        weights[-1] = l_sym
        for w in range(len(weights)-1):
            weights[w] = l_sub / (len(weights) - 1)
        #print(weights)
        print("SPLIT")
        exit()
        sim_scen = choices(cur_scen,k=1,weights=weights)[0]
        desc0 = sim_scen[0]
        cur_traits = [0.0] * ( 2 ** ss )
        cur_traits[desc0] = 1.0
        node.budd_marginals[0][i]=np.array(cur_traits)

        desc1 = sim_scen[1]
        cur_traits = [0.0] * ( 2 ** ss )
        cur_traits[desc1] = 1.0
        node.budd_marginals[1][i]=np.array(cur_traits)


def sim_along_hyp_anc(node, qmats, ss, trait_ind = None, l_sub = SUB_RATE, l_jump = JUMP_RATE):
    if trait_ind == None:
        trait_ind = range(len(node.disc_traits))

    clado_diffs = 0
    ana_diffs = 0
    for i in trait_ind:
        if node.parent == None:
            par_tr = [0.0] * ( 2 ** ss )
            par_st = choice([i for i in range(1, 2 ** ss)])
            par_tr[par_st] = 1.0
        else:
            par_tr = node.parent.budd_marginals[node.index_from_parent][i]
        
        par_state = [j for j in range(len(par_tr)) if par_tr[j] == 1.0][0]
        smap = smaps.get_smap(ss)
        npar = np.add.reduce(smap[par_state])

        if npar > 1:
            all_scen = bm.get_buddmat(ss)
            cur_scen = [j for j in all_scen[par_state] if j != 0]
            weights = np.zeros(len(cur_scen))
            #l_sub = SUB_RATE 
            l_sym = 1 - l_sub
            weights[-1] = l_sym

            for w in range(len(weights)-1):
                weights[w] = l_sub / (len(weights) - 1)
        elif npar == 1:
            cur_scen = [j for j in range(len(smap)) if np.add.reduce(smap[j]) == 1]
            weights = np.zeros(len(cur_scen))
            par_ind = [j for j in range(len(cur_scen)) if cur_scen[j] == par_state][0]
            #l_jump = JUMP_RATE
            l_sym = 1 - l_jump
            weights[par_ind] = l_sym
            for j in range(len(weights)):
                if j != par_ind:
                    weights[j] = l_jump / float(len(weights) - 1)

        last_state = choices(cur_scen,k=1,weights=weights)[0]
        clado_diffs += (smap[par_state] != smap[last_state]).sum()

        dt = node.lower - node.upper
        pmat = qmats.calc_single_p_mat(dt,ss)
        trans_probs = pmat[last_state]
        trans_probs[0] = 0.0
        trans_probs = [k / sum(trans_probs) for k in trans_probs]
        states = [k for k in range(len(trans_probs))]
        sim_state = choices(states,k=1,weights=trans_probs)[0]

        ana_diffs += ((smap[last_state] == 1) & (smap[sim_state] == 0)).sum() 
        # this version below is if you want to count total _differences_ not reductions in num of states
        #ana_diffs += (smap[last_state] != smap[sim_state]).sum()

        cur_traits = [0.0] * ( 2 ** ss )
        cur_traits[sim_state] = 1.0

        node.budd_marginals[0][i]=np.array(cur_traits)
        node.disc_traits[i] = np.array(cur_traits)
    #print(clado_diffs,ana_diffs)
    return clado_diffs, ana_diffs  


def sim_traits_across_tree2(tree,qmats,num_traits,ss):
    tree_utils.sort_children_by_age(tree) 
    ana   = 0
    clado = 0
    for n in tree.iternodes(0):
        trait_probs = np.zeros((num_traits,2 ** ss),dtype=np.double)
        n.disc_traits = trait_probs
        
        marginals = []
        if len(n.children) > 0:
            for i in range(len(n.children)):
                marginals.append(trait_probs)

            n.budd_marginals = np.array(marginals)
        if n.istip:
            clado_diffs, ana_diffs = sim_along_branch2(n, qmats, ss, None, SUB_RATE, JUMP_RATE)
        else:
            clado_diffs, ana_diffs = sim_along_hyp_anc(n, qmats, ss, None, SUB_RATE, JUMP_RATE)
        ana += ana_diffs
        clado += clado_diffs
    return clado, ana

    resim = get_pars_uninf(tree, num_traits)
    while True:
        if len(resim) == 0:
            break
        for n in tree.iternodes(0):
            if n.istip:
                clado_diffs, ana_diffs = sim_along_branch2(n, qmats, ss, resim, SUB_RATE, JUMP_RATE)
            else:
                clado_diffs, ana_diffs = sim_along_hyp_anc(n, qmats, ss, resim, SUB_RATE, JUMP_RATE)
        resim = get_pars_uninf(tree, num_traits)


if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("usage:" + sys.argv[0] + " <tree> <stratigraphic data> <num_traits> <max_num_states> <gain_rate> <loss_rate> <sub_proportion>")
        sys.exit()

    nwk = open(sys.argv[1],"r").readline()
    tree = tree_reader.read_tree_string(nwk)
    tree_utils.map_strat_to_tree(tree,sys.argv[2]) 
    num_traits = int(sys.argv[3])
    ss = int(sys.argv[4])

    """
    SUB_RATE = 0.06094244 * 10.0 #0.629
    JUMP_RATE = 0.0
    qmats = qmat.Qmat(0.035, 0.075)
    qmats = qmat.Qmat(3.95897301, 0.03811306) # trilobite rates
    SUB_RATE = 0.05373935 * 10 # trilobite
    #SUB_RATE = 0.09425363 * 10    # notharctus
    #qmats = qmat.Qmat(4.98, 0.08) # notharctus
    qmats = qmat.Qmat(0.01435794, 0.02261108) # barycrinus
    SUB_RATE = 0.08076945 * 10  # barycrinus
    """
    qmats = qmat.Qmat(float(sys.argv[5]), float(sys.argv[6]))
    SUB_RATE = float(sys.argv[7]) * 10.0    
    for row in qmats.get_qmat(2):
        print(list(row))
    stratlike.calibrate_brlens_strat(tree, 0.3)
    tree_utils.sort_children_by_age(tree) 


    ratios = []    
    for _ in range(500):
        c, a = sim_traits_across_tree2(tree,qmats,num_traits,ss)
        #print("Total clado diffs:", c)
        #print("Total ana diffs:", a)
        ratios.append((c / float(a+c)) * 100.)
    median_val = np.median(ratios)

    plt.figure(figsize=(8, 5))
    plt.hist(ratios, bins=30, color='skyblue', edgecolor='black')
    plt.axvline(median_val, color='red', linestyle='dashed', linewidth=2, label=f'Median = {median_val:.2f}')
    plt.xlabel('Percentage of Cladogenetic Changes')
    plt.ylabel('Frequency')
    plt.title('Histogram of Ratios')
    plt.legend()
    plt.tight_layout()
    plt.show()

    #for i in trait_dict:
    #    print(">"+i)
    #    print(" ".join(trait_dict[i]))

