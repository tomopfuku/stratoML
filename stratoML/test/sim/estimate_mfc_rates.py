import sys
sys.path.append("../../")
import node
import tree_reader,read_fasta,tree_utils
import numpy as np
import qmat, mfc
from scipy.optimize import minimize
from random import random
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("usage: "+ sys.argv[0]+ " <newick> <simulated trait fastas> <stratigraphic data> <stratigraphic model>")
        sys.exit()

    nwk = open(sys.argv[1],"r").readline()
    tree = tree_reader.read_tree_string(nwk)
    tree_utils.map_strat_to_tree(tree,sys.argv[3])    
    pathlist = Path(sys.argv[2]).rglob('*.fa')
    #true_m = 0.02
    #true_f = 0.02
    m_inf = []
    f_inf = []

    true_m_ls = []
    true_f_ls = []

    m_cont = []
    f_cont = []
    #print("true_m est_m true_f est_f")
    for fl in pathlist:
        flstr = str(fl).strip().split("/")[1]
        pref = flstr.strip().replace(".fa","")
        spls = pref.strip().split("_")
        true_m = float(spls[1])
        true_m_ls.append(true_m)
        true_f = float(spls[2])
        true_f_ls.append(true_f)
        traits,ss = read_fasta.read_fasta(fl)
        retraits  = read_fasta.recode_poly_traits(traits,ss)
        tree_utils.map_tree_disc_traits(tree,retraits,ss)
        qmats = qmat.Qmat(random()*.1,random()*.1)
        for n in tree.iternodes():
            n.update_pmat(qmats,max(ss))

        res_tr = minimize(mfc.evaluate_m_l,x0=np.array([0.05,0.05]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.2),(0.00001,0.2)))
        m = res_tr.x[0]
        f = res_tr.x[1]
        #print("0.02 " + str(m) + " 0.02 " + str(l))
        m_cont.append(m - true_m)
        f_cont.append(f - true_f)
        m_inf.append(m)
        f_inf.append(f)
        #aic = tree_utils.single_tree_aic(tree,ss,sys.argv[4])
        #print(aic,tree.get_newick_repr()+";")


    cont =  pd.DataFrame(
            {"m_rate":m_cont,
             "f_rate":f_cont}
            )
    
    #print(dat.melt())
    #cont.plot(kind="box")
    fig = sns.boxplot(data = cont.melt(), x = 'variable', y = 'value', color = "grey")
    fig.set(xlabel="",ylabel="true and est. rate difference")
    plt.show()

    rates = pd.DataFrame(
            {"true_m":true_m_ls,
             "infer_m":m_inf,
             "true_f":true_f_ls,
             "infer_f":f_inf
             }
            )

    biggest = max(true_m_ls + m_inf) + 0.01 #rates.to_numpy().max() + 0.01
    x=np.linspace(0.,biggest,101) 

    plt.plot(rates["true_m"],rates["infer_m"],"o",color="grey")
    plt.xlabel("true mut. rate")
    plt.ylabel("inferred mut. rate")
    plt.xlim(0.,biggest)
    plt.ylim(0.,biggest)
    plt.plot(x,x,'k-')
    plt.show()

    biggest = max(true_f_ls + f_inf) + 0.01 #rates.to_numpy().max() + 0.01
    x=np.linspace(0.,biggest,101) 

    plt.plot(rates["true_f"],rates["infer_f"],"o",color="grey")
    plt.xlabel("true fix. rate")
    plt.ylabel("inferred fix. rate")
    plt.xlim(0.,biggest)
    plt.ylim(0.,biggest)
    plt.plot(x,x,'k-')
    plt.show()

   

