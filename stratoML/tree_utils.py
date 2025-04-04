import sys
import numpy as np
import random
import stratlike
import node
from scipy.optimize import minimize
from scipy.optimize import basinhopping 
import mfc
import qmat
import biparts as bp
#import read_fasta, tree_reader


def map_tree_disc_traits(tree,traits,ss):
    ntrait = len(list(traits.values())[0])
    for n in tree.iternodes(1):
        if n.istip:
            try:
                curtr = traits[n.label]
                n.add_disc_traits(curtr,ss)
                #for i in n.disc_traits:
                #    print(list(i))
            except:
                print(n.label, "is present in the tree, but has no traits in the fasta")
                sys.exit()
            #print(curtr)
            #print(list(n.disc_starts))
            #print(list(n.disc_states))
        else:
            n.disc_traits = np.zeros((ntrait,128))

#def get_all_reattachment_nodes(tree,pluck_node):
#    for n in tree.iternodes(order=0):

       

def prune_subtree(pluck_node):
    if pluck_node.parent.parent == None:
        print("ERROR: cannot prune branch from the root. exiting search")
        sys.exit()
    if pluck_node.parent.istip:
        pluck_node.parent.children.remove(pluck_node)
        par = pluck_node.parent
        pluck_node.parent = None
        sib = None
    else:
        #print("HERE2")
        par = pluck_node.parent.parent
        sib = pluck_node.get_sib()
        par.children.remove(pluck_node.parent)
        par.children.append(sib)
        pluck_node.parent.children.remove(sib)
        sib.parent = par
        pluck_node.parent.parent = None
    return par, sib


def regraft_subtree(pluck_node,regraft_node):
    if regraft_node.istip and pluck_node.lower < regraft_node.lower:   # if the pruned subtree root is younger than the regraft point, make the regraft a direct ancestor
        pluck_node.parent = regraft_node
        regraft_node.children.append(pluck_node)
    else:
        #print("HERE")
        if pluck_node.parent == None:
            make_new_hyp_anc(pluck_node)
        pluck_node.parent.parent = regraft_node.parent
        regraft_node.parent.children.append(pluck_node.parent)
        regraft_node.parent.children.remove(regraft_node)
        regraft_node.parent = pluck_node.parent
        pluck_node.parent.children.append(regraft_node)

def reattach_to_orig_parent(pluck_node,orig_parent,sib):
    if pluck_node.parent.istip:
        pluck_node.parent.children.remove(pluck_node)
        if sib != None:
            make_new_hyp_anc(pluck_node)
    else:
        par = pluck_node.parent.parent
        rm_sib = pluck_node.get_sib()
        par.children.append(rm_sib)        
        par.children.remove(pluck_node.parent)
        pluck_node.parent.children.remove(rm_sib)
        rm_sib.parent = par
    
    if sib != None:
        orig_parent.children.append(pluck_node.parent)
        orig_parent.children.remove(sib)
        pluck_node.parent.parent = orig_parent
        sib.parent = pluck_node.parent
        pluck_node.parent.children.append(sib)
    else:
        pluck_node.parent = orig_parent
        orig_parent.children.append(pluck_node)
    #return

def check_bif(n):
    if n.parent != None:
        return True
    else:
        return False

def regraft_bif_subtree(pluck_node,regraft_node):
    if pluck_node.parent == None:
        #print("problem in regraft_bif_subtree")
        #sys.exit()
        newpar = make_new_hyp_anc(pluck_node)

    pluck_node.parent.parent = regraft_node.parent
    regraft_node.parent.children.append(pluck_node.parent)
    regraft_node.parent.children.remove(regraft_node)
    regraft_node.parent = pluck_node.parent
    pluck_node.parent.children.append(regraft_node)

def regraft_AD_subtree(pluck_node,regraft_node):
        pluck_node.parent = regraft_node
        regraft_node.children.append(pluck_node)

"""
def pick_spr_prune_regraft_nodes(tree):
    nodes = [n for n in tree.iternodes() if n != tree and n.parent != tree]
    pluck_node = random.choice(nodes)
    nodes = [n for n in tree.iternodes() if n != tree and n != prev_par and n != sib and n.subtree == pluck_node.subtree]
    if len(nodes) == 0:
        pick_spr_prune_regraft_nodes(tree)
"""


def find_best_spr(tree,ss,startaic = None,tree_mod="bds"):
    if startaic == None:
        startaic,_,_ = calc_tree_ll(tree,ss,tree_mod)
    bestaic = startaic
    for i in range(5):
        nodes = [n for n in tree.iternodes() if n != tree and n.parent != tree]
        pluck_node = random.choice(nodes)
        prev_par, sib = prune_subtree(pluck_node)
        nodes = [n for n in tree.iternodes() if n != tree and n != prev_par and n != sib and n.subtree == pluck_node.subtree]
        if len(nodes) > 0:
            break
        if sib:
            regraft_bif_subtree(pluck_node,sib)
        else:
            regraft_AD_subtree(pluck_node,prev_par)
        if i == 4:
            print("couldn't find suitable node to prune in find_best_spr()")
            return None,None

    aics = []

    for regraft_node in nodes:
        regraft_bif_subtree(pluck_node,regraft_node)
        curaic,_,_ = calc_tree_ll(tree,ss,tree_mod)
        prune_subtree(pluck_node)
        aics.append((curaic,regraft_node))
    
    aics = sorted(aics, key = lambda x: x[0])
    bestrearr = aics[0]
    changed = False
    if bestrearr[0] < startaic:
        #prune_subtree(pluck_node)
        regraft_bif_subtree(pluck_node,bestrearr[1])
        bestaic = bestrearr[0]
        changed = True
    else:
        #reattach_to_orig_parent(pluck_node,prev_par,sib)
        if sib:
            regraft_bif_subtree(pluck_node,sib)
        else:
            regraft_AD_subtree(pluck_node,prev_par)
        #sys.exit()
    return changed, bestaic

def find_all_possible_descendant_nodes(tree):
    descnodes = [n for n in tree.iternodes() if n != tree and n.parent != tree and n.istip]
    goodnodes = []
    for n in descnodes:
        ancnodes = get_possible_ancestors(n,tree)#[nn for nn in tree.iternodes() if nn != n and nn != n.parent and nn != n.get_sib() and nn.istip and nn.strat[0] >= n.strat[0] and nn.subtree == n.subtree]
        #print("TEST",n.label,len(ancnodes),[n.label for n in ancnodes])
        if len(ancnodes) > 0:
            goodnodes.append(n)
    return goodnodes

"""
def pick_candidate_ancestor(tree):
    descnodes = [n for n in tree.iternodes() if n != tree and n.parent != tree and n.istip]
    desc = random.choice(descnodes)
    ancnodes = [n for n in tree.iternodes() if n !=  and n.istip and n.strat[0] >= pluck_node.strat[0]]
    if len(ancnodes) == 0:
        pick_candidate_ancestor(tree)
"""

def choose_descnode(descnodes):
    all_root = True
    for n in descnodes:
        if n.parent.parent != None:
            all_root = False
    if all_root:
        print("ERROR: all candidate descendants are attached to root")
        sys.exit()
    pluck_node = random.choice(descnodes)
    if pluck_node.parent.parent == None:
        pluck_node = choose_descnode(descnodes)
    return pluck_node

def get_possible_ancestors(pluck_node,tree):
    nodes = [nn for nn in tree.iternodes() if nn != pluck_node and nn != pluck_node.parent and nn.istip and nn.strat[0] >= pluck_node.strat[0] and nn.subtree == pluck_node.subtree] #and nn != pluck_node.get_sib()] 
    return nodes

def find_new_ancestor(tree,ss,startaic=None,tree_mod="bds"):
    descnodes = find_all_possible_descendant_nodes(tree)
    if startaic == None:
        startaic,_,_ = calc_tree_ll(tree,ss,tree_mod)
    bestaic = startaic
    pluck_node = choose_descnode(descnodes)    # random.choice(descnodes)
    prev_par, sib = prune_subtree(pluck_node)
    if sib:
        spare_hyp_anc = pluck_node.parent
    nodes = get_possible_ancestors(pluck_node,tree) #[nn for nn in tree.iternodes() if nn != pluck_node and nn != pluck_node.parent and nn != pluck_node.get_sib() and nn.istip and nn.strat[0] >= pluck_node.strat[0] and nn.subtree == pluck_node.subtree]
    if len(nodes) == 0:
        #print("FOUND NO ANCESTORS",pluck_node.label,pluck_node.parent,pluck_node.parent.label,tree)
        print("FOUND NO ANCESTORS")
        if sib:
            pluck_node.parent = spare_hyp_anc
            regraft_bif_subtree(pluck_node,sib)
        else:
            regraft_AD_subtree(pluck_node,prev_par)
        return False,bestaic

    aics = []
    #allseen = []
    for regraft_node in nodes:
        regraft_AD_subtree(pluck_node,regraft_node)
        curaic,_,_ = calc_tree_ll(tree,ss,tree_mod)
        #allseen.append((curaic,tree.get_newick_repr()))
        prune_subtree(pluck_node)
        aics.append((curaic,regraft_node))
    aics = sorted(aics, key = lambda x: x[0])
    bestrearr = aics[0]
    changed = False
    if bestrearr[0] < bestaic:
        regraft_AD_subtree(pluck_node,bestrearr[1])
        bestaic = bestrearr[0]
        changed = True

    else:
        if sib:
            pluck_node.parent = spare_hyp_anc
            regraft_bif_subtree(pluck_node,sib)
        else:
            regraft_AD_subtree(pluck_node,prev_par)
    return changed,bestaic


def random_spr(tree):
    nodes = [n for n in tree.iternodes() if n != tree and n.parent != tree]
    try:
        pluck_node = random.choice(nodes)
    except:
        print("choosing node to prune failed")
        print(len(nodes))
        print(nodes)
        sys.exit()
    prev_par, sib = prune_subtree(pluck_node)
    nodes = [n for n in tree.iternodes() if n != tree and n != prev_par]
    regraft_node = random.choice(nodes)
    regraft_subtree(pluck_node,regraft_node)
    return pluck_node,prev_par,sib

def calc_tree_ll(tree,ss,tree_model="bds"):
    stratlike.calibrate_brlens_strat(tree,0.2)
    qmats = qmat.Qmat(0.01,0.05)
    #for n in tree.iternodes():
    #    n.update_pmat(qmats,max(ss))

    if tree_model == "bds":
        pqr_start = np.array([0.5,0.5,1.0])
        res = basinhopping(stratlike.bds_neg_ll,x0=pqr_start,niter=5,minimizer_kwargs={"method":"L-BFGS-B","args":(tree),"bounds":((0.00001,10),(0.00001,10),(0.00001,20))})
        p = res.x[0]
        q = res.x[1]
        r = res.x[2]
        stratlike.bds_dates(p,q,r,tree)
        bdsll = stratlike.bds_loglike(p,q,r,tree)
        res_tr = minimize(mfc.evaluate_m_l,x0=np.array([0.1,0.0002]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.2),(0.00001,0.2)))
        traitll = -res_tr.fun
        tree_ll = traitll + bdsll
        nparam = 3.0 + 2.0
    elif tree_model == "hr97":
        res_st = minimize(stratlike.poisson_neg_ll,x0=np.array([1.0]),args=(tree),method="Nelder-Mead")
        bdsll = -res_st.fun
        res_tr = minimize(mfc.evaluate_m_l,x0=np.array([0.1,0.0002]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.2),(0.00001,0.2)))
        print(res_tr.x)
        traitll = -res_tr.fun
        tree_ll = traitll + bdsll
        nparam = 1.0 + 2.0
        #print(res_st.fun,bdsll,res_tr.fun)
        #print(stratlike.poisson_loglike(1.0,tree))
        #print(res_st)
        #sys.exit()
    else:
        print("stratigraphic range model not recognized. please type either \"bds\" or \"hr97\"")
        sys.exit()

    nparam += float(len([n for n in tree.iternodes()]) - 1)
   
    aic = (2. * nparam) - (2. * tree_ll) 

    return aic,traitll,bdsll


# this function tries to insert hypothetical ancestors into all possible AD pairs
# TODO need to fix so that it works for ancestors with MULTIPLE descendants
def search_bifurcating(tree,ss,startaic = None,tree_mod="bds"):
    if startaic == None:
        startaic,_,_ = calc_tree_ll(tree,ss,tree_mod)

    bestaic = startaic
    testnodes = []
    for n in tree.iternodes():
        if n == tree:
            continue
        if n.istip and n.parent.istip:
            testnodes.append(n)
    aics = []
    changed = False
    if len(testnodes) > 0:
        for n in testnodes:
            make_bifurcating(n)
            curaic,_,_ = calc_tree_ll(tree,ss,tree_mod)
            aics.append((curaic,n))
            make_ancestor(n.parent)

        aics = sorted(aics, key = lambda x: x[0])
        seenpar = {}
        for tup in aics:
            curbif = tup[1]
            try: 
                seenpar[curbif]
                continue
            except:
                seenpar[curbif] = True

            make_bifurcating(curbif)
            curaic,morph,strat = calc_tree_ll(tree,ss,tree_mod)

            if curaic < bestaic:
                bestaic = curaic
                changed = True
            else:
                make_ancestor(curbif.parent)
                break
    return changed, bestaic



def search_ancestors(tree,ss,startaic = None,tree_mod="bds"):
    if startaic == None:
        startaic,_,_ = calc_tree_ll(tree,ss,tree_mod)

    testnodes = []
    for n in tree.iternodes():
        if n == tree:
            continue
        if n.istip == False:
            n_real_ch = 0
            real = None
            for chi in range(len(n.children)):
                if n.children[chi].istip:
                    n_real_ch += 1
                    real = chi
            if n_real_ch == 0:
                continue
            #elif n_real_ch == 1 and n.children[real].strat[0] < n.children[real^1].lower:
            elif n_real_ch == 1 and n.children[real].strat[0] < max([nn.strat[0] for nn in n.children[real^1].iternodes()]):
                #print("HERE")
                continue
            testnodes.append(n)
    #print(testnodes)

    aics = []
    changed = False
    bestaic = startaic
    if len(testnodes) > 0:
        for n in testnodes:
            desc = make_ancestor(n)
            curaic,_,_ = calc_tree_ll(tree,ss,tree_mod)

            aics.append((curaic,n))
            make_bifurcating(desc,n)

        aics = sorted(aics, key = lambda x: x[0])

        for tup in aics:
            curanc = tup[1]
            make_ancestor(curanc)
            curaic,morph,strat = calc_tree_ll(tree,ss,tree_mod)
            if curaic < bestaic:
                bestaic = curaic
                changed = True
            else:
                make_bifurcating(desc,curanc)
                break
        #print(bestaic)
        #print(tree.get_newick_repr())
    return changed, bestaic

def tree_search(tree, ss):
    #sl.calibrate_brlens_strat(tree)
    bestll,_,_ = calc_tree_ll(tree,ss)
    besttree = tree.get_newick_repr()
    for i in range(100):
        print(i,bestll)
        pluck_node, prev_par, sib = random_spr(tree)
        curll,morphll,bdsll = calc_tree_ll(tree,ss)
        if curll > bestll:
            reattach_to_orig_parent(pluck_node,prev_par,sib)
        else:
            bestll = curll
            besttree = tree.get_newick_repr()
    print(bestll)
    print(besttree)

def single_tree_aic(tree,ss,tree_mod="bds",anc_start=False):
    stratlike.calibrate_brlens_strat(tree,0.3)
    bestaic,traitll,bdsll = calc_tree_ll(tree,ss,tree_mod)
    #print("STARTING AIC:",bestaic)
    #print("traitll stratll treell")
    #print(traitll,bdsll,traitll+bdsll)
    return bestaic

def tree_search2(tree,ss,tree_mod="bds",anc_start=False):
    stratlike.calibrate_brlens_strat(tree,0.3)
    bestaic,_,_ = calc_tree_ll(tree,ss,tree_mod)
    besttree = tree.get_newick_repr()
    tipdic = bp.get_tip_indices(tree)
    curbp = bp.decompose_tree(tree,tipdic)
    seen = set([curbp])
    nums = [0,1,2,3] # 0 = spr; 1 = find_new_ancestor; 2 = search_ancestors; 3 = search_bifurcating
    #weights = [0.4,0.4,0.1,0.1]
    #weights = [0.45,0.45,0.1,0.0]
    weights = [0.0,1.0,0.0,0.0]
    #weights = [0.,1.,0.,0.]
    outfl = open("stratoML.outtrees","w")
    outfl.write(str(bestaic)+" "+besttree+"\n")
    lastchange = 0
    if anc_start == True:
        changed, curaic = search_ancestors(tree,ss,bestaic,tree_mod)
    else:
        changed = False
        curaic = bestaic
    for i in range(200):
        if i-lastchange >=100:
            break
        move = random.choices(population=nums,weights=weights,k=1)[0] 
        if move == 0:
            changed, curaic = find_best_spr(tree,ss,bestaic,tree_mod)
        elif move == 1:
            changed, curaic = find_new_ancestor(tree,ss,bestaic,tree_mod)
        elif move == 2: 
            changed, curaic = search_ancestors(tree,ss,bestaic,tree_mod)
        elif move == 3:
            changed, curaic = search_bifurcating(tree,ss,bestaic,tree_mod)
        else:
            print("need to specify a valid move")
            sys.exit()
        print("ITERATION",i,changed,curaic,move)
        if curaic < bestaic:
            bestaic = curaic
            besttree = tree.get_newick_repr()
            curbp = bp.decompose_tree(tree,tipdic)
            lastchange = i
            print("CURRENT:",bestaic,tree.get_newick_repr())
            if curbp not in seen:
                seen.add(curbp)
                outfl.write(str(bestaic)+" "+besttree+"\n")
    print("BEST",bestaic,besttree)
    outfl.close()
    return besttree, bestaic

def make_new_hyp_anc(descendant):
    hyp_anc = node.Node()
    descendant.parent = hyp_anc
    descendant.parent.disc_traits = np.zeros((len(descendant.disc_traits),128))
    hyp_anc.children.append(descendant)
    return hyp_anc

def make_bifurcating(descendant,hyp_anc = None):
    if descendant.parent == None:
        print("trying to transform AD pair into bifurcating but there is no ancestor!")
        print(descendant.label)
        sys.exit()
    elif descendant.parent.istip == False:
        print("trying to transform AD pair into bifurcating but the ancestor is already hypothetical!")
        print(descendant.label)
        sys.exit()

    if hyp_anc == None:
        hyp_anc = node.Node()
        hyp_anc.disc_traits = np.zeros((len(descendant.disc_traits),128))

    curpar = descendant.parent
    if curpar.parent != None:
        curpar.parent.children.append(hyp_anc)
        #print(curpar.label,curpar.parent.label)
        curpar.parent.children.remove(curpar)
        hyp_anc.parent = curpar.parent

    curpar.children.remove(descendant)
    curpar.parent = hyp_anc
    descendant.parent = hyp_anc
    hyp_anc.children.append(curpar) 
    hyp_anc.children.append(descendant) 

    oldest = max([i.lower for i in hyp_anc.children])
    hyp_anc.upper = oldest 
    hyp_anc.lower = oldest + 0.1
    for n in hyp_anc.children:
        n.lower = oldest

def make_ancestor(hyp_anc):
    if hyp_anc.istip:
        print("cannot make bifurcating pair into anc-desc from already sampled ancestor")
        sys.exit()
    elif len(hyp_anc.children) != 2:
        print("cannot make AD pair from hyp_anc with more than 2 children")
        sys.exit()
    elif hyp_anc.children[0].istip == False and hyp_anc.children[1].istip == False:
        print("cannot make two bifurcating nodes AD pair")
        sys.exit()

    if hyp_anc.children[0].istip == False:
        anc = hyp_anc.children[1]
    elif hyp_anc.children[1].istip == False:
        anc = hyp_anc.children[0]
    elif hyp_anc.children[0].lower == hyp_anc.children[1].lower:
        anc = random.choice(hyp_anc.children)
    else:
        oldest = 0.0
        anc = None
        for ch in hyp_anc.children:
            if ch.lower > oldest:
                anc = ch
                oldest = ch.lower

    desc = anc.get_sib()
    desc.parent = anc
    anc.children.append(desc)
    hyp_anc.children = []
    if hyp_anc.parent != None:
        hyp_anc.parent.children.append(anc)
        anc.parent = hyp_anc.parent
        hyp_anc.parent.children.remove(hyp_anc)
    return desc    


def random_nni(tree):
    keepnodes = []
    #for n in tree.iternodes():
    #    if le 
    # 

def map_strat_to_tree(tree, flnm):
    fl = open(flnm, "r")
    fl.readline()
    ranges = {}

    for line in fl:
        spls = line.strip().split()
        spnm = spls[0]
        fad = float(spls[1])
        lad = float(spls[2])
        ranges[spnm] = np.array([fad,lad])

    for n in tree.iternodes():
        if n.istip:
            try:
                strat = ranges[n.label]
                n.strat = strat
            except:
                print(n.label," is in the tree but was not found in the stratigraphic range data")
                sys.exit()

    stratlike.calibrate_brlens_strat(tree)
