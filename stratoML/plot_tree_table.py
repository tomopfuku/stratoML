import node,tree_reader,tree_utils
import stratlike
import sys
import numpy as np
from scipy.optimize import basinhopping 

def assign_codes(tree):
    nn = 0
    codes = {}    
    for n in tree.iternodes():
        codes[n] = str(nn)
        nn+=1
    return codes

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage:",sys.argv[0]," <tree> <stratigraphic ranges>")
        sys.exit()

    nwk = open(sys.argv[1],"r").readline().strip()
    tree = tree_reader.read_tree_string(nwk)
    tree_utils.map_strat_to_tree(tree,sys.argv[2])
    stratlike.calibrate_brlens_strat(tree,0.4)
    pqr_start = np.array([0.5,0.5,1.0])
    res = basinhopping(stratlike.bds_neg_ll,x0=pqr_start,niter=20,minimizer_kwargs={"method":"L-BFGS-B","args":(tree),"bounds":((0.00001,10),(0.00001,10),(0.00001,20))})
    p = res.x[0]
    q = res.x[1]
    r = res.x[2]

    stratlike.bds_dates(p,q,r,tree)

    print("speciation:",p,"\nextinction:",q,"\npreservation",r,"\n")
    #for n in tree.iternodes():
    #    print(n.label,n.lower,n.upper,n.strat[0],n.strat[1])

    codes = assign_codes(tree)

    node_cis = stratlike.bds_CIs(p,q,r,tree)

    outfl = open(".".join(sys.argv[1].strip().split(".")[0:-1])+".tree_table","w")
    outfl.write("Name,Code,Start,End,FAD,LAD,2.5HPD,97.5HPD,Parent\n")
    for n in tree.inorder():
        curcode = codes[n]
        if n.parent != None:
            parcode = codes[n.parent]
        else:
            parcode = "NA"
        if n.istip:
            outfl.write(n.label+","+curcode+","+str(n.lower)+","+str(n.upper)+","+str(n.strat[0])+","+str(n.strat[1])+","+str(node_cis[n][0])+","+str(node_cis[n][1])+","+parcode+"\n")
        else:
            outfl.write(n.label+","+curcode+","+str(n.lower)+","+str(n.upper)+","+str(n.lower)+","+str(n.upper)+","+str(node_cis[n][0])+","+str(node_cis[n][1])+","+parcode+"\n")
    outfl.close()