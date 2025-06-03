import node,tree_reader,tree_utils
import stratlike
import sys
import numpy as np
from scipy.optimize import basinhopping
#import bd
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from main_prob_tree_anc import *


if len(sys.argv) < 3:
    print("usage:",sys.argv[0]," <tree> <stratigraphic ranges> <OPTIONAL: COMPLETENESS LEVEL TO CALCULATE>")
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
completeness = r / (r + q)
print("completeness:",completeness)

stratlike.bds_dates(p,q,r,tree)

if len(sys.argv) == 4:
    completeness = float(sys.argv[3])
    r = abs ( (completeness * q) / (completeness - 1.) )
    print("SIMULATED preservation rate (based on "+str(completeness)+" completeness): " +str(r)) 

alltax   = []
allprobs = []
allpointprobs = []
tax_expect = {}
for n in tree.iternodes():
    if n.istip:
        n_chld = len(n.children)
        obs_range = n.strat[0] - n.strat[1]
        duration = n.length
        probs = calc_cum_desc_probs(p,q,r,n_chld,obs_range,duration)
        alltax.append(n.label.strip().split("_")[-1])
        allprobs.append(np.array(probs))
        #print(n.label)
        probs = calc_desc_probs(p,q,r,n_chld,obs_range,duration)
        allpointprobs.append(np.array(probs[0:9]))
        tax_expect[n.label] = [n_chld] + expect_desc(p,q,r,n_chld,obs_range,duration) + [duration]

expect_df = pd.DataFrame(tax_expect).transpose().rename(columns={0:"n_obs",1:"expect_mle",2:"expect_weighted",3:"duration"})
expect_df["n_below_exp"] = expect_df["n_obs"] - expect_df["expect_weighted"]


flpref = ".".join(sys.argv[1].strip().split(".")[0:-1]) + "." + str(round(completeness * 100))
expect_df.to_csv(flpref+".expect_desc",sep=",")

plt.scatter(expect_df["duration"],expect_df["expect_weighted"],color="black")
plt.scatter(expect_df["duration"],expect_df["n_obs"],color="grey")
#plt.xscale("log")
#plt.yscale("log")
plt.show()


print(expect_df)
allpointprobs = np.array(allpointprobs)
col = []
for i in range(len(allpointprobs[0])):
    if i == 0:
        col.append("P("+str(i)+")")
    else:
        col.append("P("+str(i)+")")

allpdf = pd.DataFrame(allpointprobs,columns = col,index = alltax)
n_sig_above = 0
n_sig_below = 0
for n in tree.iternodes():
    if n.istip:
        expect = list(expect_df.loc[[n.label]]["expect_weighted"])[0]
        sp = n.label.strip().split("_")[-1]
        nchld=len(n.children)
        row = allpdf.loc[[sp]]
        p=row.iloc[0,nchld]
        if p < 0.1:
            if nchld > expect:
                print(n.label,nchld,expect)
                n_sig_above += 1
            else:
                n_sig_below += 1

print("clade n_sig_above n_sig_below")
print(sys.argv[1].strip().split(".")[0],sys.argv[1].strip().split(".")[1], n_sig_above, n_sig_below)

