import sys
import os
import multiprocessing as mp
from contextlib import contextmanager
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import node
import tree_reader,read_fasta,tree_utils,stratlike,mfc, glc_bd
import numpy as np
import qmat
import lam_mat
from scipy.optimize import minimize
from scipy.optimize import basinhopping
import time
import math
import emcee
import corner
import matplotlib.pyplot as plt
import pandas as pd

_worker_tree = None
_worker_qmats = None
_worker_lam_mats = None
_worker_ss = None
_worker_pqr_start = None
MCMC_SAVE_DISCARD = 100
MCMC_SAVE_THIN = 15


def init_worker(tree, qmats, lam_mats, ss, pqr_start):
    global _worker_tree, _worker_qmats, _worker_lam_mats, _worker_ss, _worker_pqr_start
    _worker_tree = tree
    _worker_qmats = qmats
    _worker_lam_mats = lam_mats
    _worker_ss = ss
    _worker_pqr_start = np.array(pqr_start, dtype=np.float64)


def log_probability(params, tree, qmats, lam_mats, ss, pqr_start):
    params_array = np.array(params, dtype=np.float64)
    if np.any(params_array <= 1e-5) or np.any(params_array > 10.0):
        return -np.inf
    try:
        ll = -glc_bd.evaluate_m_l3(
            params_array, 
            tree, 
            qmats, 
            lam_mats, 
            ss, 
            np.array(pqr_start, dtype=np.float64)
        )
        return ll
    except Exception:
        return -np.inf


def log_probability_worker(params):
    return log_probability(
        params,
        _worker_tree,
        _worker_qmats,
        _worker_lam_mats,
        _worker_ss,
        _worker_pqr_start,
    )


def get_max_mcmc_workers(nwalkers):
    default_workers = min(nwalkers, os.cpu_count() or 1)
    max_workers = os.environ.get("MCMC_MAX_THREADS", os.environ.get("MCMC_NPROCS"))

    if max_workers is None:
        return default_workers

    try:
        max_workers = int(max_workers)
    except ValueError:
        print(f"Ignoring invalid MCMC_MAX_THREADS/MCMC_NPROCS value: {max_workers!r}")
        return default_workers

    return max(1, min(nwalkers, max_workers))


@contextmanager
def mcmc_pool(tree, qmats, lam_mats, ss, pqr_start, nwalkers):
    nproc = get_max_mcmc_workers(nwalkers)

    if nproc == 1:
        yield None
        return

    if "fork" not in mp.get_all_start_methods():
        print("Fork multiprocessing is unavailable; running MCMC serially to avoid pickling Cython Node objects.")
        yield None
        return

    ctx = mp.get_context("fork")
    pool = ctx.Pool(
        processes=nproc,
        initializer=init_worker,
        initargs=(tree, qmats, lam_mats, ss, pqr_start),
    )
    try:
        print(
            f"Running MCMC with {nproc} forked worker processes "
            f"for {nwalkers} walkers."
        )
        yield pool
    finally:
        pool.close()
        pool.join()


def load_restart_positions(chain_file, nwalkers, ndim):
    if not os.path.exists(chain_file):
        raise FileNotFoundError(f"Restart chain file does not exist: {chain_file}")

    chain_df = pd.read_csv(chain_file)
    column_names = [f"param_{i}" for i in range(ndim)]
    missing_columns = [col for col in column_names if col not in chain_df.columns]
    if missing_columns:
        raise ValueError(
            "Restart chain is missing expected columns: "
            + ", ".join(missing_columns)
        )
    if len(chain_df) < nwalkers:
        raise ValueError(
            f"Restart chain has {len(chain_df)} samples, but {nwalkers} walkers are needed."
        )

    positions = chain_df[column_names].tail(nwalkers).to_numpy(dtype=np.float64)
    if positions.shape != (nwalkers, ndim):
        raise ValueError(
            f"Restart positions have shape {positions.shape}; expected {(nwalkers, ndim)}."
        )
    if np.any(~np.isfinite(positions)):
        raise ValueError("Restart chain contains non-finite parameter values.")
    if np.any(positions <= 1e-5) or np.any(positions > 10.0):
        raise ValueError("Restart chain contains parameter values outside the MCMC prior bounds.")

    return positions, chain_df[column_names].copy()


def estimate_saved_generations(nsamples, nwalkers):
    if nsamples <= 0:
        return 0
    saved_steps = int(math.ceil(float(nsamples) / float(nwalkers)))
    return MCMC_SAVE_DISCARD + saved_steps * MCMC_SAVE_THIN


if __name__ == "__main__":
    if len(sys.argv) not in (7, 8):
        print("usage: "+ sys.argv[0]+ " <newick> <trait fasta file> <stratigraphic data> <stratigraphic model> <morphologic model> <num_gen> [previous mcmc samples csv]")
        sys.exit()

    try:
        num_gen = int(sys.argv[6])
    except ValueError:
        print(f"num_gen must be an integer; received {sys.argv[6]!r}")
        sys.exit()
    if num_gen <= 0:
        print(f"num_gen must be positive; received {num_gen}")
        sys.exit()

    restart_chain_file = sys.argv[7] if len(sys.argv) == 8 else None
    
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

        #mono_prob = 0.6
        #tree_utils.fix_obs_lv(tree, True, True, ss, mono_prob) 
        tree_utils.fix_obs_lv(tree) 
        qmats = qmat.Qmat(0.05,0.05)

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
        #p = q = 0.1
        pqr_start = [p,q,r]
        print("BDS rates:",pqr_start)
        sub_frac = 0.5
        jump_frac = 0.0
        lsub = p * sub_frac 
        #lsym = p - lsub 
        ljump = p * jump_frac 
        lam_mats = lam_mat.lam_mat(p, lsub)#, ljump)

        N = float(len([n for n in tree.iternodes() if n.istip == True])) * ntraits 
        start_rates = np.array([0.05, 0.05, sub_frac / 10.])

        ndim = 3           # Number of parameters
        nwalkers = 9      # Number of parallel "walkers" (usually 2x or 4x ndim)
        column_names = [f"param_{i}" for i in range(ndim)]
        previous_samples = None
        if restart_chain_file is None:
            p0 = abs(np.array(start_rates) + 1.5e-2 * np.random.randn(nwalkers, ndim)) # Initial spread
        else:
            p0, previous_samples = load_restart_positions(restart_chain_file, nwalkers, ndim)
            print(
                f"Restarting MCMC from the last {nwalkers} samples in "
                f"'{restart_chain_file}'."
            )
        total_num_gen = num_gen
        if previous_samples is not None:
            total_num_gen += estimate_saved_generations(len(previous_samples), nwalkers)
        
        # 3. Run the MCMC
        print("Running MCMC...")
        with mcmc_pool(tree, qmats, lam_mats, ss, pqr_start, nwalkers) as pool:
            if pool is None:
                sampler = emcee.EnsembleSampler(
                    nwalkers,
                    ndim,
                    log_probability,
                    args=(tree, qmats, lam_mats, ss, pqr_start),
                )
            else:
                sampler = emcee.EnsembleSampler(
                    nwalkers,
                    ndim,
                    log_probability_worker,
                    pool=pool,
                )
            sampler.run_mcmc(
                p0,
                num_gen,
                progress=True,
                skip_initial_state_check=previous_samples is not None,
            )
        
        # --- 1. Convergence Analysis ---
        print("\n--- Convergence Diagnostics ---")

        # Integrated Autocorrelation Time
        try:
            tau = sampler.get_autocorr_time()
            print(f"Autocorrelation time for each parameter: {tau}")
            
            # Rule of thumb: chain length should be >> 50 * tau
            for i, t in enumerate(tau):
                if total_num_gen > 50 * t:
                    print(f"  Parameter {i}: Chain length looks sufficient (> 50 * tau).")
                else:
                    print(f"  Parameter {i}: Warning! Chain might be too short ({total_num_gen} steps < 50 * tau).")
        except emcee.autocorr.AutocorrError as e:
            print("Autocorrelation time warning:", e)
            print("The chain is likely too short or hasn't converged enough to reliably estimate autocorrelation.")

        # Acceptance Fraction
        # Ideally should be between 0.2 and 0.5 for emcee
        acc_fraction = sampler.acceptance_fraction
        print(f"Mean acceptance fraction: {np.mean(acc_fraction):.3f}")
        if np.any(acc_fraction < 0.1) or np.any(acc_fraction > 0.6):
            print("Warning: Some walkers have low/high acceptance fractions. Check your initial spread or model.")


        # --- 2. Extract and Save Flat Samples to CSV ---
        # Your original discard and thin parameters
        save_discard = 0 if previous_samples is not None else MCMC_SAVE_DISCARD
        flat_samples = sampler.get_chain(discard=save_discard, thin=MCMC_SAVE_THIN, flat=True)

        # Convert to a Pandas DataFrame for easy CSV exporting
        # (Adjust column names to match your actual parameters if desired)
        new_df = pd.DataFrame(flat_samples, columns=column_names)
        if previous_samples is None:
            df = new_df
        else:
            df = pd.concat([previous_samples, new_df], ignore_index=True)

        # Construct CSV filename from sys.argv[1]
        output_base = sys.argv[2] 
        csv_filename = f"{output_base}_mcmc_samples.csv"
        df.to_csv(csv_filename, index=False)

        print(f"\nSuccessfully saved {len(df)} flat samples to '{csv_filename}'.")
        if previous_samples is not None:
            print(f"  Appended {len(new_df)} new samples to {len(previous_samples)} restart samples.")

        #flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)

        # Get the "Best Fit" (Median) and 16th/84th percentiles
        flat_samples = df.to_numpy(dtype=np.float64)
        for i in range(ndim):
            mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
            print(f"Param {i}: {mcmc[1]:.4f} (+{mcmc[2]-mcmc[1]:.4f} / -{mcmc[1]-mcmc[0]:.4f})")
        
        # 2. Corner Plot: The 4D "Diagnosis"
        labels = ["gain rate", "loss rate", "sub ratio"]
        fig = corner.corner(
            flat_samples, 
            labels=labels,
            quantiles=[0.16, 0.5, 0.84],
            show_titles=True,
            title_fmt=".4f",
            color="forestgreen"
        )

        #plt.suptitle("3-Parameter Model Posterior", fontsize=16)
        plt.show()

        exit()


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
        """
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
        #res_tr = minimize(glc_bd.evaluate_m_l3,x0=np.array(start_rates),args=(tree,qmats,lam_mats,ss, np.array(pqr_start)),method="Powell")
        for row in list(lam_mats.get_ratemat(2)):
            print(list(row))
        print(res_tr)
        m1_ll = -res_tr.fun
        m1_aic = (2. * k) - (2.  * m1_ll)
        m1_bic = (k * math.log(N)) - (2. * m1_ll)
        #start_rates = np.array([0.05, 0.05, sub_frac, jump_frac])
        #nll = glc_bd.evaluate_m_l_jump(start_rates, tree, qmats, lam_mats, ss, np.array(pqr_start))
        
        print("no jump params:", res_tr.x)
        t1 = time.time()
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
