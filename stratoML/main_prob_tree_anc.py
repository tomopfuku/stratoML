import node,tree_reader,tree_utils
import stratlike
import sys
import numpy as np
from scipy.optimize import basinhopping
import bd
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def calc_desc_probs(p,q,r,n_chld,obs_range,duration):
    pvec = []
    for n_chld in range(50):
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
        probs.append(curp)
    return(probs)

def expect_desc(p,q,r,n_chld,obs_range,duration):
    probs = calc_desc_probs(p,q,r,n_chld,obs_range,duration)
    expect = [i for i,j in enumerate(probs) if j == max(probs)][0]
    weight_expect = 0.0
    for i, p in enumerate(probs):
        weight_expect += (i + 1) * p
    return [expect, weight_expect - 1.0]

def calc_cum_desc_probs(p,q,r,n_chld,obs_range,duration):
    probs = calc_desc_probs(p,q,r,n_chld,obs_range,duration)
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

print("speciation rate:",p)
print("extinction rate:",q)
print("preservation rate:",r)
print("completeness:",r / (r + q))

stratlike.bds_dates(p,q,r,tree)

alltax   = []
allprobs = []
tax_expect = {}
for n in tree.iternodes():
    if n.istip:
        n_chld = len(n.children)
        obs_range = n.strat[0] - n.strat[1]
        duration = n.length
        probs = calc_cum_desc_probs(p,q,r,n_chld,obs_range,duration)
        alltax.append(n.label.strip().split("_")[-1])
        allprobs.append(np.array(probs))
        tax_expect[n.label] = [n_chld] + expect_desc(p,q,r,n_chld,obs_range,duration)
        #print(n.label," ".join([str(round(i,2)) for i in probs]))


"""
expect_count_mat = {0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[]}

for exp in range(8):
    expect_count_mat[exp] = [0 for i in range(len(expect_count_mat))]


for i in tax_expect:
    ex = tax_expect[i]
    obs = ex[0]
    exp = ex[1]

    expect_count_mat[exp][obs] += 1

print(expect_count_mat)
expect = pd.DataFrame(expect_count_mat).transpose()#,columns=["E"+str(i) for i in range(8)],index = ["O"+str(i) for i in range(8)])
print(expect)
with sns.axes_style("white"):
    ax=sns.heatmap(expect,annot=True,cmap="Greys",cbar=False)
    ax.xaxis.tick_top()
    ax.tick_params(left=False, top=False)
plt.show()
"""

expect_df = pd.DataFrame(tax_expect).transpose().rename(columns={0:"n_obs",1:"expect_mle",2:"expect_weighted"})
plt.plot(expect_df["expect_weighted"],expect_df["n_obs"],"o")
plt.ylabel("num observed desc")
plt.xlabel("weighted expected num desc")
plt.show()


plt.plot(np.linspace(-0.5,7.5,20),np.linspace(-0.5,7.5,20),color="grey")
sns.boxplot(data=expect_df,y="n_obs",x="expect_weighted",orient="h",color="grey")
#sns.boxplot(data=expect_df,x="n_obs",y="expect_weighted",color="grey")
plt.xlim((-0.5,7.5))
plt.ylim((-0.5,7.5))
plt.show()

expect_df = pd.DataFrame(tax_expect).transpose().rename(columns={0:"n_obs",1:"expect_mle",2:"expect_weighted"})
plt.plot(expect_df["expect_mle"],expect_df["n_obs"],"o")
plt.ylabel("num observed desc")
plt.xlabel("mle expected num desc")
plt.show()


allprobs = np.array(allprobs)
col = []
for i in range(len(allprobs[0])):
    if i == 0:
        col.append("P("+str(i)+")")
    else:
        col.append("P("+str(i)+"+)")
allpdf = pd.DataFrame(allprobs,columns = col,index = alltax)
#seaborn.heatmap(allprobs,show=True)
#plt = seaborn.heatmap(allprobs)
with sns.axes_style("white"):
    ax=sns.heatmap(allpdf,annot=True,cmap="Greys",cbar=False)
    ax.xaxis.tick_top()
    ax.tick_params(left=False, top=False)

    #ax.set_title(title)
    #plt.savefig(title+".mod.svg",format="svg")
plt.show()
