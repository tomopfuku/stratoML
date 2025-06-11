import sys
sys.path.append("../../")
import node,tree_reader,tree_utils
import stratlike
import numpy as np
import math
import biparts as bp
from scipy.optimize import basinhopping



def check_all_bps(tree_bp,cur_bp):
    if len(cur_bp) != len(tree_bp):
        return False
    elif len(cur_bp.symmetric_difference(tree_bp)) > 0:
        return False
    elif len(cur_bp.symmetric_difference(tree_bp)) == 0:
        return True
 
def assign_codes(tree):
    nn = 0
    codes = {}    
    for n in tree.iternodes():
        codes[n.label] = str(nn)
        nn+=1
    return codes

def label_internal_nodes(tree):
    count = 0
    for n in tree.iternodes():
        if n.istip == False:
            lab = "n"+str(count)
            n.label = lab
            count += 1


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage:",sys.argv[0]," <tree> <stratigraphic ranges> <top_trees>")
        sys.exit()

    nwk = open(sys.argv[1],"r").readline().strip()
    tree = tree_reader.read_tree_string(nwk)
    tree_utils.map_strat_to_tree(tree,sys.argv[2])
    label_internal_nodes(tree)
    stratlike.calibrate_brlens_strat(tree,0.4)

    #for n in tree.iternodes():
    #    print(n.label,n.lower,n.upper,[(ch.label,ch.lower) for ch in n.children])
    #sys.exit()
    pqr_start = np.array([0.5,0.5,1.0])
    res = basinhopping(stratlike.bds_neg_ll,x0=pqr_start,niter=20,minimizer_kwargs={"method":"L-BFGS-B","args":(tree),"bounds":((0.00001,10),(0.00001,10),(0.00001,20))})
    p = res.x[0]
    q = res.x[1]
    r = res.x[2]

    #stratlike.bds_dates(p,q,r,tree)

    print("speciation:",p,"\nextinction:",q,"\npreservation",r,"\n")
    #for n in tree.iternodes():
    #    print(n.label,n.lower,n.upper,n.strat[0],n.strat[1])

    codes = assign_codes(tree)

    node_cis = stratlike.bds_CIs(p,q,r,tree)

    nl_cis = {}
    for n in node_cis:
        nl_cis[n.label] = node_cis[n]

    node_cis = nl_cis


    ############################
    ## NOW MAP SUPPORT VALUES ##
    ############################

    fl = open(sys.argv[3],"r")
    trees = []
    for line in fl:
        spls = line.strip().split()
        aic = float(spls[0])
        nwk = spls[1]
        curtree = tree_reader.read_tree_string(nwk)
        trees.append((aic,curtree))

    best = 100000000000.
    besttree = None
    for i in trees:
        if i[0] < best:
            best = i[0]
            besttree = i[1]


    tipdic = bp.get_tip_indices(tree)
    curbp = bp.decompose_tree(tree,tipdic)
    origbp = bp.decompose_tree(tree,tipdic)
    bpsupp = {}
    for i in curbp:
        bpsupp[i]=0

    bpall = []
    sumrellike = 0.
    origAIC = None
    for i in trees:
        aic = i[0]
        rellike = math.exp(0.5*(best-aic))
        sumrellike+=rellike
        curtree = i[1]
        curbp = bp.decompose_tree(curtree,tipdic)
        same_tree = check_all_bps(origbp,curbp)
        if same_tree:
            origAIC = rellike
        bpall.append((rellike,curbp))

    for best_bp in bpsupp:
        cumaic = 0.
        for i in bpall:
            rellike = i[0]
            decomp_tree = i[1]
            for curbp in decomp_tree:
                #print(best_bp,curbp)
                if best_bp == curbp:
                    weight = rellike / sumrellike
                    bpsupp[best_bp] += weight


    tree_supp = origAIC / sumrellike
    print("TREE_SUPPORT:",tree_supp)

    node_supports = {}
    for n in tree.iternodes():
        if n == tree or len(n.children) == 0:
            node_supports[n.label] = "NA" 
            continue

        bs = bp.make_bipart(besttree,n,tipdic)
        for i in bpsupp:
            if i == bs:
                node_supports[n.label] = bpsupp[i]
                #print(n.label,bs,bpsupp[i])



    #print(node_supports)
    outfl = open(".".join(sys.argv[1].strip().split(".")[0:-1])+".tree_table","w")
    outfl.write("Name,Code,Start,End,FAD,LAD,2.5HPD,97.5HPD,Parent,Support\n")
    for n in tree.inorder():
        curcode = codes[n.label]
        if n.parent != None:
            parcode = codes[n.parent.label]
        else:
            parcode = "NA"
        if n.istip:
            outfl.write(n.label+","+curcode+","+str(n.lower)+","+str(n.upper)+","+str(n.strat[0])+","+str(n.strat[1])+","+str(node_cis[n.label][0])+","+str(node_cis[n.label][1])+","+parcode+","+str(node_supports[n.label])+"\n")
        else:
            #outfl.write(n.label+","+curcode+","+str(n.lower)+","+str(n.upper)+","+str(n.lower)+","+str(n.upper)+","+str(node_cis[n.label][0])+","+str(node_cis[n.label][1])+","+parcode+","+str(node_supports[n.label])+"\n")
            outfl.write(""+","+curcode+","+str(n.lower)+","+str(n.upper)+","+"NA"+","+"NA"+","+"NA"+","+"NA"+","+parcode+","+str(node_supports[n.label])+"\n")
    outfl.close()


