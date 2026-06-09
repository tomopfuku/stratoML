import sys
import node
import tree_reader,read_fasta,tree_utils,stratlike,mfc, glc_bd
import numpy as np
import qmat
import lam_mat
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy.optimize import differential_evolution        
import time
import math

def evaluate_m_l3_wrapper(params, *args):
    return glc_bd.evaluate_m_l3(params, *args)

def evaluate_m_l_jump_wrapper(params, *args):
    return glc_bd.evaluate_m_l_jump(params, *args)

def get_tree_length(tree):
    length = 0.0
    for n in tree.iternodes():
        length += (n.lower - n.upper)
    return length

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("usage: "+ sys.argv[0]+ " <newick> <trait fasta file> <stratigraphic data> <stratigraphic model> <morphologic model>")
        sys.exit()
    
    traits,ss = read_fasta.read_fasta(sys.argv[2])
    ntraits = float(len(list(traits.values())[0]) - 1)
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

        #mono_prob = 0.5
        #tree_utils.fix_obs_lv(tree, True, True, ss, mono_prob) 
        tree_utils.fix_obs_lv(tree) 
        qmats = qmat.Qmat(0.01,0.05)

        for n in tree.iternodes():
            n.update_pmat(qmats,max(ss), "mid")


        pqr_start = [0.1, 0.1, 1.1]
        res = basinhopping(stratlike.bds_neg_ll,x0=pqr_start,niter=20,minimizer_kwargs={"method":"L-BFGS-B","args":(tree),"bounds":((0.00001,10),(0.00001,10),(0.00001,20))})
        p = res.x[0]
        q = res.x[1]
        r = res.x[2]
        #p = q = 0.1
        pqr_start = [p,q,r]
        print("BDS rates:",pqr_start)
        sub_frac = 0.5
        jump_frac = 0.1
        lsub = p * sub_frac 
        #lsym = p - lsub 
        ljump = p * jump_frac 
        lam_mats = lam_mat.lam_mat(p, lsub, 0.0)#ljump)

        #for n in tree.iternodes():
        #    print(n.label, n.midpoint)
        #exit()
        
        #for row in qmats.get_qmat(3):
        #    print(list(row))
        #print("LAM:", sub_frac, lsub, p - lsub)
        #for row in list(lam_mats.get_ratemat(4)):
        #    print(sum(row))
        #    print(list(row))
        #exit()
        N = float(len([n for n in tree.iternodes() if n.istip == True])) * ntraits 
        #print(res.fun)
        start_rates = np.array([0.05, 0.05, sub_frac / 10.])
        k = len(start_rates)
        #res_tr = minimize(glc_bd.evaluate_m_l3,x0=np.array(start_rates),args=(tree,qmats,lam_mats,ss, np.array(pqr_start)),method="L-BFGS-B",bounds=((0.00001,5.0),(0.00001,5.0),(0.00001,1.0)))
        #res_tr = minimize(glc_bd.evaluate_m_l3,x0=np.array(start_rates),args=(tree,qmats,lam_mats,ss, np.array(pqr_start)),method="Nelder-Mead")
        """
        res_tr = minimize(
            glc_bd.evaluate_m_l3,
            x0=np.array(start_rates),
            args=(tree, qmats, lam_mats, ss, np.array(pqr_start)),
            method="L-BFGS-B",
            bounds=((1e-5, 5.0), (1e-5, 5.0), (1e-5, 1.0)),
            options={
                'eps': 1e-4,       # Larger step for gradient estimation
                'ftol': 1e-7,      # Slightly more relaxed function tolerance
                'gtol': 1e-5,      # Gradient tolerance
                'maxls': 50,       # Increase max line search steps
            }
        )

        print("no jump params:", res_tr.x)
        exit()
        minimizer_kwargs = {
            "method": "L-BFGS-B",
            "args": (tree, qmats, lam_mats, ss, np.array(pqr_start)),
            "bounds": ((1e-5, 5.0), (1e-5, 5.0), (1e-5, 1.0)),
            "options": {"eps": 1e-4}
        }

        res_tr = basinhopping(
            glc_bd.evaluate_m_l3,
            x0=np.array(start_rates),
            minimizer_kwargs=minimizer_kwargs,
            niter=20,     # Try 20 different starting locations
            stepsize=0.5  # Magnitude of the "jump" between iterations
        )
        """
        res_tr = differential_evolution(
            glc_bd.evaluate_m_l3,
            bounds=((1e-5, 5.0), (1e-5, 5.0), (1e-5, 1.0)),
            args=(tree, qmats, lam_mats, ss, np.array(pqr_start)),
            strategy='best1bin', # Standard reliable strategy
            popsize=15,          # Total candidates = popsize * parameters
            mutation=(0.5, 1),   # Exploration vs exploitation
            recombination=0.7,
            tol=0.01,            # Relative tolerance for convergence
            workers=1           # USES ALL CORES (This is the speed booster!)
        )
        #res_tr = minimize(glc_bd.evaluate_m_l3,x0=np.array(start_rates),args=(tree,qmats,lam_mats,ss, np.array(pqr_start)),method="Powell")
        m1_ll = -res_tr.fun
        m1_aic = (2. * k) - (2.  * m1_ll)
        m1_bic = (k * math.log(N)) - (2. * m1_ll)
        #start_rates = np.array([0.05, 0.05, sub_frac, jump_frac])
        #nll = glc_bd.evaluate_m_l_jump(start_rates, tree, qmats, lam_mats, ss, np.array(pqr_start))
        
        print("no jump params:", res_tr.x)

        gainr = res_tr.x[0]
        lossr = res_tr.x[1]
        subfrac = res_tr.x[2] * 10.0 
        lsub = p * subfrac
        print("GAIN:", gainr)
        print("LOSS:", lossr)
        print("LAMBDA SUB:", lsub)
        print("no jump AIC:",m1_aic)
        #tr_len = get_tree_length(tree)
        exit()
        
        t1 = time.time()
        """
        #res_tr = minimize(glc_bd.evaluate_m_l3,x0=np.array(start_rates),args=(tree,qmats,lam_mats,ss, np.array(pqr_start)),method="L-BFGS-B",bounds=((0.00001,5.0),(0.00001,5.0),(0.00001,1.0)))
        start_rates = np.array([0.05, 0.05, sub_frac / 10., jump_frac / 10.])
        #res_tr = minimize(glc_bd.evaluate_m_l_jump,x0=np.array(start_rates),args=(tree,qmats,lam_mats,ss, np.array(pqr_start)),method="L-BFGS-B",bounds=((0.00001,5.0),(0.00001,5.0),(0.00001,1.0),(0.00001,1.0)))
        minimizer_kwargs = {
            "method": "L-BFGS-B",
            "args": (tree, qmats, lam_mats, ss, np.array(pqr_start)),
            "bounds": ((1e-5, 5.0), (1e-5, 5.0), (1e-5, 1.0), (1e-5, 1.0)),
            "options": {"eps": 1e-4}
        }

        res_tr = basinhopping(
            glc_bd.evaluate_m_l_jump,
            x0=np.array(start_rates),
            minimizer_kwargs=minimizer_kwargs,
            niter=20,     # Try 20 different starting locations
            stepsize=0.5  # Magnitude of the "jump" between iterations
        )
        """
        res_tr = differential_evolution(
            glc_bd.evaluate_m_l_jump,
            bounds=((1e-5, 5.0), (1e-5, 5.0), (1e-5, 1.0/10.), (1e-5, 1.0/10.)),
            args=(tree, qmats, lam_mats, ss, np.array(pqr_start)),
            strategy='best1bin', # Standard reliable strategy
            popsize=15,          # Total candidates = popsize * parameters
            mutation=(0.5, 1),   # Exploration vs exploitation
            recombination=0.7,
            tol=0.01,            # Relative tolerance for convergence
            workers=1           # USES ALL CORES (This is the speed booster!)
        )


        t2 = time.time()
        m2_ll = -res_tr.fun
        k += 1
        m2_aic = (2. * k) - (2.  * m2_ll)
        m2_bic = (k * math.log(N)) - (2. * m2_ll)
        print("no jump AIC:",m1_aic)
        print("no jump BIC:",m1_bic)
        print("\n\nJUMP MODEL\n\n")
        print(res_tr)

        print("jump AIC:",m2_aic)
        print("jump BIC:",m2_bic)
        print("jump params:",res_tr.x)
        sys.exit()
        
        ### NEED TO IMPLEMENT FULL TREE LIKE (INCL STRAT MODEL) 
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
