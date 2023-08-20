import node,tree_reader,tree_utils
import stratlike
import sys
import numpy as np
from scipy.optimize import basinhopping 
import bd
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def calc_all_desc_probs(p,q,r,n_chld,obs_range,duration):
    pvec = []
    for n_chld in range(50):
        #r_like = np.log(bd.calc_prob_range(r,duration,obs_range))
        #p_like = np.log(bd.prob_n_obs_desc(p,q,r,n_chld,duration))
        #q_like = np.log(bd.prob_extinction_t(q,duration))
        #ll = sum([r_like,p_like,q_like])
        ll = np.log(bd.prob_n_obs_desc(p,q,r,n_chld,duration))
        pvec.append(ll)

        
    best = max(pvec)
    rel = []
    for i in pvec:
        rll = np.exp(-.5*(best-i))
        rel.append(rll)

    sumrel = np.sum(rel)

    probs = []
    for n_chld,val in enumerate(rel):
        curp = val/sumrel
        #probs.append(str(round(val/sumrel,3)))
        probs.append(curp)

    cump = []
    for n_chld,val in enumerate(probs):
        if n_chld == 0:
            curcum = val
        else:
            curcum = np.sum(probs[n_chld:])

        #cump.append(str(round(curcum,2)))
        cump.append(curcum)
        if n_chld > 3:
            break
    return cump


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

alltax   = []
allprobs = []
for n in tree.iternodes():
    if n.istip:
        n_chld = len(n.children)
        obs_range = n.strat[0] - n.strat[1]
        duration = n.length
        probs = calc_all_desc_probs(p,q,r,n_chld,obs_range,duration)
        alltax.append(n.label.strip().split("_")[-1])
        allprobs.append(np.array(probs)) 
        print(n.label," ".join([str(round(i,2)) for i in probs]))

allprobs = np.array(allprobs)
col = []
for i in range(len(allprobs[0])):
    if i == 0:
        col.append("P("+str(i)+")")
    else:
        col.append("P("+str(i)+"+)")
allpdf = pd.DataFrame(allprobs,columns = col,index = alltax)
print(allpdf)
#seaborn.heatmap(allprobs,show=True)
#plt = seaborn.heatmap(allprobs)
with sns.axes_style("white"):
    ax=sns.heatmap(allpdf,annot=True,cmap="Greys",cbar=False)
    ax.xaxis.tick_top()
    ax.tick_params(left=False, top=False)

    #ax.set_title(title)
    #plt.savefig(title+".mod.svg",format="svg")
plt.show()