import sys
import node
import tree_reader,read_fasta,tree_utils,stratlike,mfc, glc_bd
import numpy as np
import qmat
import lam_mat
from scipy.optimize import minimize
from scipy.optimize import basinhopping
import matplotlib.pyplot as plt
import time


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("usage: "+ sys.argv[0]+ " <newick> <trait fasta file> <stratigraphic data> <stratigraphic model> <morphologic model>")
        sys.exit()
    
    traits,ss = read_fasta.read_fasta(sys.argv[2])
    retraits  = read_fasta.recode_poly_traits(traits,ss)

    times = []
    for line in open(sys.argv[1],"r"):
        if line.strip() == "":
            continue

        nwk = line.strip().split()[-1]
        tree = tree_reader.read_tree_string(nwk)

        tree_utils.map_strat_to_tree(tree,sys.argv[3])    
        #stratlike.calibrate_brlens_strat(tree,0.3)
        tree_utils.map_tree_disc_traits(tree,retraits,ss)
        tree_utils.sort_children_by_age(tree)
        tree_utils.init_budd_marginals(tree,len(ss), ss)

        #mono_prob = 0.6
        #tree_utils.fix_obs_lv(tree, True, True, ss, mono_prob) 
        tree_utils.fix_obs_lv(tree) 
        qmats = qmat.Qmat(0.01,0.05)

        for n in tree.iternodes():
            n.update_pmat(qmats,max(ss), "mid")

        #lm2 = lam_mats.get_ratemat(2)
        #lm3 = lam_mats.get_ratemat(3)
        #for row in lm3:
        #    print(list(row))
        pqr_start = [0.1, 0.1, 1.1]
        res = basinhopping(stratlike.bds_neg_ll,x0=pqr_start,niter=20,minimizer_kwargs={"method":"L-BFGS-B","args":(tree),"bounds":((0.00001,10),(0.00001,10),(0.00001,20))})
        p = res.x[0]
        q = res.x[1]
        r = res.x[2]
        pqr_start = [p,q,r]
        print("BDS rates:",pqr_start)
        X = np.linspace(0.005, 0.2, 10)
        Y = np.linspace(0.005, 0.2, 10)
        print(X,Y)
        Xgrid, Ygrid = np.meshgrid(X, Y)
        z = []
        highest = -10000000000.0
        highest_vals = []
        count = 0
        for yi in Y: 
            y_ll = []
            for xi in X:
                print(count)
                count += 1
                lsub = p * xi 
                lsym = p - lsub 
                lam_mats = lam_mat.lam_mat(p, lsub)

                #start_rates = np.array([yi, 0.09, xi]) # gain rates
                start_rates = np.array([xi, yi, 0.5 / 10.])  # loss rates
                ll = -glc_bd.evaluate_m_l3(start_rates, tree, qmats, lam_mats, ss, np.array(pqr_start))
                if ll > highest:
                    highest = ll
                    highest_vals = [xi, yi]
                y_ll.append(ll)
                
            z.append(y_ll)
        
        print(highest_vals)
        levels = [highest - 20, highest - 10, highest - 5, highest - 1, highest - .1] 
        fig, ax = plt.subplots()
        CS = ax.contour(Xgrid,Ygrid,z,levels=levels)
        ax.clabel(CS, inline=True, fontsize=10)
        plt.xlabel("sub clado fraction")
        plt.ylabel("loss rate")
        plt.show()

        sys.exit()
        t1 = time.time()
        aic,traitll,bdsll = tree_utils.calc_tree_ll2(tree,qmats,ss,"hr97")
        t2 = time.time()
        print(line.strip().split()[0:-1], aic)
        #print(aic,nwk)
        #print("TIME OPTIMIZING",t2-t1)
        times.append(t2-t1)

    print("MEAN TIME:",   np.mean(times))
    print("MEDIAN TIME:", np.median(times))
    print("TOTAL TIME:", sum(times))
