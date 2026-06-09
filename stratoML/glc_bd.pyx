# cython: language_level=3
# cython: initializedcheck=True

import node
cimport node
import numpy as np
#cimport numpy as cnp
cimport numpy as np
np.import_array()
cimport cython
import sys
from scipy.linalg import expm
import smaps
import qmat
cimport qmat
import lam_mat
cimport lam_mat
import bd
import mfc
import spltmat as sm
import buddmat as bm
import stratlike
import tree_utils
import math


cpdef double[:] max_scale_timeslice(double[:] anc_marg):
    cdef double[:] scaled_anc_marg
    cdef int i

    #print("IN MAX SCSALE TIMESLICE", list(anc_marg))
    scaled_anc_marg = np.zeros(len(anc_marg))
    for i in range(len(anc_marg)):
        scaled_anc_marg[i] = anc_marg[i] / max(anc_marg)

    return scaled_anc_marg

cdef calc_M_matrix(double dt, qmat.Qmat qmats, lam_mat.lam_mat lam_mats, double[:] bds_rates):
    return

cdef calc_prob_surv(double dt, double[:] bds_rates):
    cdef double E_t, p_zero_desc, p_surv, surv_scalar 
    #p_zero_desc = bd.prob_n_obs_desc(bds_rates[0], bds_rates[1], bds_rates[2], 0, dt)
    E_t = bd.calc_extinction_prob_eq(bds_rates[0], bds_rates[1], bds_rates[2])
    p_zero_desc = math.exp(-(bds_rates[0] * ( 1.0 - E_t ) * dt))
    p_surv = math.exp(-bds_rates[1] * dt)
    surv_scalar = p_zero_desc * p_surv
    return surv_scalar

def segment_branch(double[:] times, double[:] strat_range):
    cdef int MAX_SEGMENTS = 5 
    cdef np.ndarray[double, ndim=1] durations = np.empty(MAX_SEGMENTS, dtype=np.float64)
    cdef np.ndarray[int, ndim=1] types = np.empty(MAX_SEGMENTS, dtype=np.int32)
    cdef double cur_t, start, end, fad, lad, duration, cutoff
    cdef bint is_obs
    cdef int model_type, count = 0

    start = times[0]
    end = times[1]
    fad = strat_range[0]
    lad = strat_range[1]

    #NOTE NEED TO CHANGE THIS BELOW lad - cutoff
    cutoff = 1.0
    cur_t = end 
    while cur_t < start:
        is_obs = (cur_t >= lad - 0.0) and (cur_t < fad + 0.0)
        if is_obs:
            #model_type = "Q"
            model_type = 0
            duration = min(start, fad) - cur_t
        else:
            #model_type = "M"
            model_type = 1
            if cur_t < lad:
                duration = lad - end
            elif cur_t >= fad:
                duration = start - fad
        #curseg = Segment(duration, model_type)
        durations[count] = duration
        types[count] = model_type
        cur_t += duration
        count+=1
        if count > MAX_SEGMENTS:
            print("ERROR IN glc_bd.segment_branch()")
            exit()
    #print(list([float(i) for i in durations[:count]]))
    #print(list(types[:count]))
    return durations[:count], types[:count]

cdef calc_like_to_base(double[:, :] last_tr, qmat.Qmat qmats, lam_mat.lam_mat lam_mats, long[:] ss, double[:] bds_rates, double[:] seg_durations, int[:] seg_models):
    cdef Py_ssize_t i, chari, cur_k, nstates, j, k
    cdef double cur_dt, p_surv, val
    cdef int model_type
    cdef Py_ssize_t nseg = seg_durations.shape[0]
    cdef Py_ssize_t nch = ss.shape[0]
    cdef Py_ssize_t ncol = last_tr.shape[1]
    cdef double[:, :] working_tr = np.empty((nch, ncol), dtype=np.float64)
    cdef double[:] D_anc_buffer = np.empty(ncol, dtype=np.float64)
    cdef double[:, :] curp
    cdef object pmats

    # Copy last_tr to working_tr
    for chari in range(nch):
        for j in range(ncol):
            working_tr[chari, j] = last_tr[chari, j]

    for i in range(nseg):
        cur_dt = seg_durations[i]
        model_type = seg_models[i]

        if model_type == 0:
            pmats = qmats.calc_p_mats(cur_dt)
            p_surv = calc_prob_surv(cur_dt, bds_rates)
        elif model_type == 1:
            pmats = qmats.calc_m_mats(cur_dt, lam_mats, bds_rates)

        for chari in range(nch):
            cur_k = ss[chari]
            nstates = 1 << cur_k  # Faster than 2 ** cur_k

            if cur_k == 1:
                continue

            curp = pmats[cur_k - 2]

            for j in range(nstates):
                val = 0.0
                for k in range(nstates):
                    val += curp[j, k] * working_tr[chari, k]
                if model_type == 0:
                    val *= p_surv
                D_anc_buffer[j] = val

            # Use slice assignment for speed
            for j in range(ncol):
                if j < nstates:
                    working_tr[chari, j] = D_anc_buffer[j]
                else:
                    working_tr[chari, j] = 0.0

    return np.asarray(working_tr)

cdef calc_like_to_base_OLD(double[:, :] last_tr, qmat.Qmat qmats, lam_mat.lam_mat lam_mats, long[:] ss, double[:] bds_rates, double[:] seg_durations, int[:] seg_models):
    for i, cur_dt in enumerate(seg_durations):
        model_type = seg_models[i]
        
        if model_type == 0:
            pmats = qmats.calc_p_mats(cur_dt)
        elif model_type == 1:
            pmats = qmats.calc_m_mats(cur_dt, lam_mats, bds_rates)

        for chari in range(0, len(ss)):
            cur_k = ss[chari]
            curp = pmats[cur_k - 2]
            nstates = int(2 ** cur_k)
            if cur_k == 1:
                continue

            prev_tr = last_tr[chari]
            p_surv = calc_prob_surv(cur_dt, bds_rates)
            D_anc = np.array(curp[:nstates, :nstates]) @ np.array(prev_tr[:nstates])
            print(list(seg_durations))
            for j in range(nstates):
                print(j,list(curp[j]))
            print(D_anc)
            print("PREVTRAIT",list(prev_tr)) 
            print("PSURV", p_surv)
            if model_type == 0:
                D_anc = p_surv * D_anc
            for j in range(len(last_tr[chari])):
                if j < len(D_anc):
                    last_tr[chari][j] = D_anc[j]
                else:
                    last_tr[chari][j] = 0.0
    return last_tr

cpdef double[:, :] missing_trait_vec_all(long[:] ss):
    cdef double[:, :] tr = np.zeros([len(ss), 128])
    cdef int j, k
    
    for j, cur_k in enumerate(ss):
        ana_tr = mfc.missing_trait_vec(cur_k)
        for k in range(len(ana_tr)):
            tr[j][k] = ana_tr[k]
    return tr

def budd_node_join(double[:,:] anc_ll, 
                   double[:, :] base_ll_ch, 
                   lam_mat.lam_mat lam_mats, 
                   long[:] ss, 
                   bint rescale = False):
    
    cdef int n_chars = ss.shape[0]
    cdef double[:,:] budd_ll_ch = np.zeros((n_chars, base_ll_ch.shape[1]), dtype=np.float64)
    
    cdef int chari, i, j, cur_k, nstates
    cdef double[:,:] mat_view
    cdef double dot_product_val
    
    for chari in range(n_chars):
        cur_k = ss[chari]
        if cur_k == 1:
            continue
            
        nstates = 1 << cur_k  # Efficient way to do 2**cur_k
        
        # Get the matrix as a memoryview to avoid Python overhead
        # Assuming get_ratemat returns something that can be viewed as double[:,:]
        mat_view = lam_mats.get_ratemat(cur_k)
        
        # Manually perform Matrix-Vector multiplication and element-wise join
        # This replaces: (Mat @ Base) * Anc
        for i in range(nstates):
            dot_product_val = 0.0
            for j in range(nstates):
                dot_product_val += mat_view[i, j] * base_ll_ch[chari, j]
            
            # Combine the result with anc_ll and store directly in the output
            budd_ll_ch[chari, i] = dot_product_val * anc_ll[chari, i]

    return budd_ll_ch

def budd_node_join_UNOPTIMIZED(double[:,:] anc_ll, double[:, :] base_ll_ch, lam_mat.lam_mat lam_mats, long[:] ss, bint rescale = False):
    cdef double[:,:] budd_ll_ch = np.zeros_like(base_ll_ch) 

    for chari in range(len(ss)):
        cur_k = ss[chari]
        nstates = int(2 ** cur_k)
        if cur_k == 1:
            continue

        bud_likes = np.array(lam_mats.get_ratemat(cur_k)) @ np.array(base_ll_ch[chari][:nstates])
        join_ll = bud_likes * anc_ll[chari][:nstates]

        for i in range(len(bud_likes)):
            budd_ll_ch[chari][i] = join_ll[i]

    return budd_ll_ch

def split_like_marg(node.Node n, 
                    qmat.Qmat qmats, 
                    lam_mat.lam_mat lam_mats, 
                    long[:] ss, 
                    double[:] bds_rates):
    
    cdef int chari, i, j, k, cur_k, nstates
    cdef node.Node ch
    cdef double[:] seg_durations
    cdef int[:] seg_models
    cdef double[:, :] base_ll_ch, mat_view
    cdef double dot_product
    cdef int n_chars = ss.shape[0]
    
    # Pre-allocate or access existing views
    # Assuming n.children[0].timeslice_lv[-1] provides the shape (n_chars, max_states)
    cdef double[:, :] join_ll = np.zeros((n_chars, n.children[0].timeslice_lv[-1].shape[1]), dtype=np.float64)
    cdef double[:] scaled_anc_marg_i

    if len(n.children) != 2:
        print("hypothetical ancestors can only have two children in the model.")
        sys.exit()

    for i in range(len(n.children)):
        ch = n.children[i]
        
        if ch.istip:
            ch_times = mfc.get_child_dt(ch, True)
            seg_durations, seg_models = segment_branch(np.array(ch_times), ch.strat)
        else:
            seg_durations = np.array([mfc.get_child_dt(ch, False)], dtype=np.float64)
            seg_models = np.array([1], dtype=np.int32)

        base_ll_ch = calc_like_to_base(ch.timeslice_lv[-1], qmats, lam_mats, ss, bds_rates, seg_durations, seg_models)

        for chari in range(n_chars):
            cur_k = ss[chari]
            if cur_k == 1:
                continue
            
            nstates = 1 << cur_k
            mat_view = lam_mats.get_ratemat(cur_k)

            for j in range(nstates):
                dot_product = 0.0
                for k in range(nstates):
                    dot_product += mat_view[j, k] * base_ll_ch[chari, k]
                
                if i == 0:
                    join_ll[chari, j] = dot_product
                else:
                    join_ll[chari, j] *= dot_product

        for chari in range(n_chars):
            # Scaling factors update
            # Assuming n.scaling_factors is an object/array accessible via indices
            n.scaling_factors[ch.parent_lv_index][chari] = max(join_ll[chari])
            
            # This is likely a bottleneck if it returns a new list/array
            scaled_anc_marg_i = mfc.max_scale_timeslice(join_ll[chari])

            for j in range(len(scaled_anc_marg_i)):
                n.timeslice_lv[ch.parent_lv_index][chari][j] = scaled_anc_marg_i[j]


cdef double split_like_marg_UNOPTIMIZED(node.Node n, qmat.Qmat qmats, lam_mat.lam_mat lam_mats, long[:] ss, double[:] bds_rates):
    cdef double nodelike, charll, dt
    cdef int chari, cur_k #, ti, tj, anc1, anc2, traitprob_i, nscenarios
    cdef node.Node ch
    cdef double[:, :, :] pmats0, pmats1
    cdef double[:, :] join_ll = np.zeros_like(n.children[0].timeslice_lv[-1])
    cdef int[:] seg_models

    if len(n.children) != 2:
        print("hypothetical ancestors can only have two children in the model.")
        sys.exit() 
    for i in range(len(n.children)):
        ch = n.children[i]
        if ch.istip:
            ch_times = mfc.get_child_dt(ch, True)
            seg_durations, seg_models = segment_branch(np.array(ch_times), ch.strat)
        else:
            seg_durations = np.array([mfc.get_child_dt(ch, False)])
            seg_models = np.array([1], dtype=np.int32)
        #print(list(seg_durations), list(seg_models)) 
        ch_tr = ch.timeslice_lv[-1]
        base_ll_ch = calc_like_to_base(ch_tr, qmats, lam_mats, ss, bds_rates, seg_durations, seg_models)
        for chari in range(len(ss)):
            cur_k = ss[chari]
            nstates = int(2 ** cur_k)
            if cur_k == 1:
                continue
            bud_likes = np.array(lam_mats.get_ratemat(cur_k)) @ np.array(base_ll_ch[chari][:nstates])
            for j in range(len(bud_likes)):
                if join_ll[chari][j] == 0.0:
                    join_ll[chari][j] = bud_likes[j]
                else:
                    join_ll[chari][j] *= bud_likes[j]
        #print(i, list(join_ll[2]))
        for chari in range(len(ss)):
            n.scaling_factors[ch.parent_lv_index][chari] = max(join_ll[chari]) 
            scaled_anc_marg_i = mfc.max_scale_timeslice(join_ll[chari])

            for j in range(len(scaled_anc_marg_i)):
                n.timeslice_lv[ch.parent_lv_index][chari][j] = scaled_anc_marg_i[j] 


cdef budd_like_marg(node.Node n, qmat.Qmat qmats, lam_mat.lam_mat lam_mats, long[:] ss, double[:] bds_rates):
    cdef int chari, cur_k, chd_i, ch_ind, num_ch = len(n.children) #, i ,j,nscen, ancst, k, maxstates
    cdef double prev_time, dt1, dt2 
    cdef node.Node ch
    cdef bint past_mid
    cdef double[:, :, :] pmats1, pmats2
    cdef double[:, :] p1, p2, ana_tr, ch_tr, base_ll_ch, base_ll_par_br, anc_marg 
    cdef double[:] seg_durations, ch_times, scaled_anc_marg_i
    cdef int[:] seg_models
    cdef list times

    times = [ch.lower for ch in n.children]
    times.append(n.midpoint)
    times.sort()
    prev_time = n.upper
    past_mid = False
    for i in range(len(times)):
        cur_time = times[i]
        #print("PREENTRY PAR", n.label, cur_time, prev_time, list(n.strat))
        seg_durations, seg_models = segment_branch(np.array([cur_time,prev_time]), n.strat)
        if i == 0:
            ana_tr = missing_trait_vec_all(ss)
        else:
            ana_tr = n.timeslice_lv[i - 1]

        base_ll_par_br = calc_like_to_base(ana_tr, qmats, lam_mats, ss, bds_rates, seg_durations, seg_models)
        if i == n.midpoint_lv_index:
            anc_marg = np.zeros_like(base_ll_par_br)
            for chari in range(len(ss)):
                cur_tr = np.asarray(n.timeslice_lv[n.midpoint_lv_index][chari]) * np.asarray(base_ll_par_br[chari])
                for j in range(len(cur_tr)):
                    anc_marg[chari][j] = cur_tr[j] 
                n.scaling_factors[n.midpoint_lv_index][chari] = max(anc_marg[chari]) 
                scaled_anc_marg_i = mfc.max_scale_timeslice(anc_marg[chari])

                for j in range(len(scaled_anc_marg_i)):
                    n.timeslice_lv[n.midpoint_lv_index][chari][j] = scaled_anc_marg_i[j] 
            past_mid = True
        else:
            if past_mid == False:
                ch_ind = num_ch - 1 - i
                ch = n.children[ch_ind] ## NEED TO FIX THIS INDEXING IMMEDIATELY. N.CHILDREN IS SORTED OLDEST TO NEWEST WHICH IS OPPOSITE OF times
                #ch = rev_ch[i]
            else: # need to do the node join
                ch_ind = num_ch - i
                ch = n.children[ch_ind]
                #ch = rev_ch[i - 1]
            if ch.istip:
                ch_times = mfc.get_child_dt(ch, True)
                seg_durations, seg_models = segment_branch(np.array(ch_times), ch.strat)
            else:
                seg_durations = np.array([mfc.get_child_dt(ch, False)])
                seg_models = np.array([1], dtype=np.int32)
                #print("LINE278")
                #print(list(seg_durations))
                #exit()
            ch_tr = ch.timeslice_lv[-1]
            base_ll_ch = calc_like_to_base(ch_tr, qmats, lam_mats, ss, bds_rates, seg_durations, seg_models)
            anc_marg = budd_node_join(base_ll_par_br, base_ll_ch, lam_mats, ss, True)
            #print([c.label for c in n.children])
            #print(n.label, ch.label, ch_ind, i, list(anc_marg[0]))
            for chari in range(len(ss)):
                n.scaling_factors[ch.parent_lv_index][chari] = max(anc_marg[chari]) 
                #print(n.label, ch.label, chari, list(anc_marg[chari]))
                scaled_anc_marg_i = mfc.max_scale_timeslice(anc_marg[chari])

                for j in range(len(scaled_anc_marg_i)):
                    n.timeslice_lv[ch.parent_lv_index][chari][j] = scaled_anc_marg_i[j] 
        prev_time = cur_time                

def calc_bud_root_ll_single_trait(node.Node tree, qmat.Qmat qmats, int cur_k, int chari):
    cdef double scen_cond_like,scenario_like,  weight, marg_prob, dt, traitprob
    cdef double[:,:] pmat
    cdef double[:] chd_tr, anc_marg
    cdef long[:,:] cur_scen
    cdef long[:] inher, ancsts
    cdef int  nscen, j, ancst, startst, k

    cur_scen = bm.get_buddmat(cur_k)

    if tree.midpoint_lv_index != len(tree.children):
        dt = tree.lower - tree.children[0].lower
    elif tree.midpoint_lv_index == len(tree.children) or len(tree.children) == 0:
        dt = tree.lower - tree.midpoint

    #print("HERE")
    pmat = qmats.calc_single_p_mat(dt, cur_k) # P matrix for segment along the child leading to the first likelihood vector

    chd_tr = tree.timeslice_lv[-1][chari] 
    #print(list(chd_tr))
    ancsts = np.array([i for i in range(1,len(cur_scen))])
    anc_marg = np.zeros(len(ancsts)+1)
    for ancst in ancsts:
        inher = cur_scen[ancst]
        #print(ancst,list(inher))
        nscen = mfc.count_nscenarios(inher)
        weight = 1.0 / ( float(nscen) * float(len(cur_scen)-1) )
        marg_prob = 0.0
        for j in range(len(inher)): # calc cond like for each allowed inheritance scen
            startst = inher[j]
            if startst == 0:
                break

            #scenario_like = 0.0 # marginal likelihood of one scenario, marginalized over all combinations of anagenetic transitions and budded descendant states

            scen_cond_like = 0.0
            for k in range(len(chd_tr)):
                traitprob = chd_tr[k]
                if k == int(2) ** cur_k: # NEED TO FIX
                    break
                if traitprob == 0.0:
                    continue
                scen_cond_like += pmat[startst][k] * traitprob
                #print("CLAD PROB FROM ",startst,"TO",k,scen_cond_like)
            scenario_like = scen_cond_like 
            
            #print("UNWEIGHTED SCENARIO LIKE",scenario_like)
            scenario_like *= weight
            #print("WEIGHTED SCENARIO LIKE",scenario_like)
            marg_prob += scenario_like
        anc_marg[ancst] = marg_prob
    #print(list(anc_marg))
    return anc_marg

def calc_bud_root_ll(node.Node tree, qmat.Qmat qmats, long[:] ss):
    cdef double[:,:] root_partials #last_lv 
    cdef int chari, cur_k, i
    cdef double[:] anc_marg

    root_partials = np.zeros((len(ss),128))
    #last_lv = tree.timeslice_lv[-1]
    for chari in range(len(ss)):
        cur_k = ss[chari]
        
        if cur_k == 1:
            continue

        anc_marg = calc_bud_root_ll_single_trait(tree, qmats, cur_k, chari)
        #print(chari,list(anc_marg))
        for i in range(len(anc_marg)):
            root_partials[chari][i] = anc_marg[i]
        #print(chari,list(root_partials[chari]))
    return root_partials

def mfc3_treell(node.Node tree, qmat.Qmat qmats, lam_mat.lam_mat lam_mats, long[:] ss, double[:] bds_rates, bint asc = True):
    cdef double treell, asc_treell, plike, flat_prior, invarll, sum_plikes, sublike
    cdef double[:,:] root_marg_likes
    cdef double[:] plikes, sum_log_sf
    cdef int i, j, k, count
    cdef node.Node n

    treell = 0.0
    for n in tree.iternodes(1):
        if len(n.children) == 0:
            continue
        if n.istip == False:
            split_like_marg(n, qmats, lam_mats, ss, bds_rates)
        elif n.istip:
            budd_like_marg(n, qmats, lam_mats, ss, bds_rates)

    if tree.istip == False:
        root_marg_likes = tree.timeslice_lv[-1]
    else:
        root_marg_likes = calc_bud_root_ll(tree, qmats, ss)

    sum_log_sf = mfc.calc_logsum_scaling_factors(tree, ss) 

    treell = 0.0
    for i in range(1, len(root_marg_likes)):
        #for j in i:
        if ss[i] == 1: # invariant traits do not contribute to the tree likelihood bc we use ascertainment bias correction
            continue
        plikes = root_marg_likes[i]
        sum_plikes = sum(plikes)
        sublike = 0.0
        count = 0
        
        for j in range(len(plikes)):
            if j > 0.0:
                count+=1
        
        #flat_prior = 1.0 / float(count)
        for j in range(len(plikes)):
            plike = plikes[j]
            sublike += plike * (plike / sum_plikes) #flat_prior 
        #print(i, np.log(sublike), np.log(sublike) + sum_log_sf[i], np.exp(np.log(sublike) + sum_log_sf[i]))
        sublike = np.log(sublike) + sum_log_sf[i]
        treell += sublike

    if asc == True:
        invarll = calc_invar_ll_marg(tree,qmats)
        #invarll = np.log(0.5)
        treell = treell - np.log(1.0-np.exp(invarll))
    return treell


def log1mexp(x):
    if x > -0.69314718056: # ln(2)
        return np.log(-np.expm1(x))
    else:
        return np.log1p(-np.exp(x))

cdef double calc_invar_ll_marg(node.Node tree, qmat.Qmat qmats):
    cdef double desc_weight, nodell, sum_log_sf, sublike, sum_plikes, plike
    cdef double[:] plikes
    cdef node.Node n
    cdef int i, j, count
    cdef double [:, :, :] pmats1, pmats2

    """
    for n in tree.iternodes(1):
        if len(n.children) == 0:
            continue

        if n.istip == False:
            dt = get_child_dt(n.children[0])
            p1 = qmats.calc_single_p_mat(dt, 2)
            dt = get_child_dt(n.children[1])
            p2 = qmats.calc_single_p_mat(dt, 2)
            split_loglike_single_trait_marg(n, p1, p2, 2, 0)
        elif n.istip:
            budd_loglike_single_trait_marg(n, qmats, 2, 0)
    """
    if tree.istip:
        plikes = calc_bud_root_ll_single_trait(tree, qmats, 2, 0)
    else:
        plikes = tree.timeslice_lv[-1][0]

    sum_log_sf = 0.0 
    for n in tree.iternodes(1):
        if n == tree and n.istip == False: 
            continue
        if len(n.children) == 0:
            continue
        for i in range(len(n.scaling_factors)):
            sum_log_sf += np.log(n.scaling_factors[i][0])
        
    sum_plikes = sum(plikes)
    sublike = 0.0
    count = 0
    for j in range(len(plikes)):
        if plikes[j] > 0.0:
            count+=1

    for j in range(len(plikes)):
        plike = plikes[j]
        sublike += plike * (plike / sum_plikes) #1.0 / count 

    sublike = np.log(sublike) + sum_log_sf
    return sublike 

def evaluate_m_l3_mono_weight(double[:] params, node.Node tree, qmat.Qmat qmats,long[:] ss):
    cdef node.Node n
    cdef double treell, weight
    cdef int maxstates = max(ss)
    if params[0] < 0.0001 or params[1] < 0.0001 or params[2] < 0.0001 or params[2] > 1.0:
        return 10000000000

    weight = params[2]

    tree_utils.fix_obs_lv(tree, True, True, ss, weight) 
    qmats.update_all_qmats(params[0],params[1]) 

    tree_utils.sort_children_by_age(tree)
    treell = mfc3_treell(tree, qmats, ss)
    return -treell

def evaluate_m_l_jump(double[:] params, node.Node tree, qmat.Qmat qmats, lam_mat.lam_mat lam_mats, long[:] ss, double[:] bds_rates):
    cdef node.Node n
    cdef double treell
    cdef int maxstates = max(ss)

    for i in range(len(params)):
        if params[i] < 0.00001:
            return 10000000000

    gainr = params[0]
    lossr = params[1]
    qmats.update_all_qmats(gainr, lossr) 

    lam = bds_rates[0]
    mu = bds_rates[1]
    psi = bds_rates[2]

    lsub_w = params[2] * 10.
    ljump_w = params[3] * 10.
    #print(gainr, lossr)
    if lsub_w > 0.9999 or ljump_w > 0.9999:
        return 10000000000

    lsub = lam * lsub_w
    ljump = lam * ljump_w
    lam_mats.update_all_mats(lam, lsub, ljump)

    tree_utils.sort_children_by_age(tree)
    #treell = mfc.mfc2_treell(tree, qmats, ss)
    treell = mfc3_treell(tree, qmats, lam_mats, ss, bds_rates)
    return -treell



def evaluate_m_l3(double[:] params, node.Node tree, qmat.Qmat qmats, lam_mat.lam_mat lam_mats, long[:] ss, double[:] bds_rates):
    cdef node.Node n
    cdef double treell
    cdef int maxstates = max(ss)

    for i in range(len(params)):
        if params[i] < 0.00001:
            print("HERE")
            return 10000000000

    gainr = params[0]
    lossr = params[1]
    qmats.update_all_qmats(gainr, lossr) 

    lam = bds_rates[0]
    mu = bds_rates[1]
    psi = bds_rates[2]

    lsub_w = params[2] * 10.
    if lsub_w > 0.99:
        return 10000000000

    lsub = lam * lsub_w
    lam_mats.update_all_mats(lam, lsub)

    tree_utils.sort_children_by_age(tree)
    #treell = mfc.mfc2_treell(tree, qmats, ss)
    treell = mfc3_treell(tree, qmats, lam_mats, ss, bds_rates)

    #print(lam, lsub_w, gainr, lossr, -treell)
    return -treell

