# cython: language_level=3
# cython: initializedcheck=True

import node
cimport node
import numpy as np
cimport numpy as np
np.import_array()
cimport cython
import sys
from scipy.linalg import expm
import smaps
import qmat
cimport qmat
import spltmat as sm
import buddmat as bm
import stratlike
import tree_utils
#from scipy.optimize import minimize
#from cython.parallel import prange
#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.locals(i=cython.int, x=cython.double)
# cython: declare_variables=True

#@cython.wraparound(False)
#@cython.boundscheck(False)

def mat_mult(double[:,:] mat, double m):
    cdef int i,j,nrow
    nrow = len(mat)
    cdef double[:,:] newmat = np.empty((nrow,nrow),dtype = np.double)
    for i in range(nrow):
        for j in range(nrow):
            newmat[i][j] = mat[i][j] * m
    return newmat

#def create_calc_p_matrix():
#    qmat.Qmat ratemats

def calc_p_matrix(double[:,:] ratemat, double t):
    cdef double[:,:] pmat
    cdef double val
    cdef int i

    pmat = expm(mat_mult(ratemat,t))
    for i in range(len(pmat[0])):
        val = pmat[0][i]
        if np.isnan(val):
            print("NAN")
            sys.exit()
    return pmat

cdef double split_loglike_single_trait_marg(node.Node n, double[:, :] p1, double[:, :] p2, int cur_k, int chari):
    cdef double traitprob, curp1, curp2, anclike, weight, dt
    cdef long[:] spltanc 
    cdef int ti, tj, anc1, anc2, traitprob_i, nscenarios, count
    cdef double[:] ch1_tr,ch2_tr
    cdef long[:,:,:] cur_scen
    cdef node.Node ch

    cur_scen = sm.get_spltmat(cur_k)
    ch = n.children[0]
    ch1_tr = ch.timeslice_lv[-1][chari]

    ch = n.children[1]
    ch2_tr = ch.timeslice_lv[-1][chari]
  
    count = 0
    for ti in range(len(cur_scen)):
        if ti == 0:
            continue
        anclike = 0.0
        nscenarios = 0
        for tj in range(len(cur_scen[ti])):
            spltanc = cur_scen[ti][tj]
            if spltanc[0] == 0: #np.add.reduce(spltanc) == 0:
                break
            nscenarios += 1

        weight = 1.0 / (float(nscenarios) * float(len(cur_scen)-1))
        for tj in range(len(cur_scen[ti])):
            spltanc = cur_scen[ti][tj]
            if spltanc[0] == 0: #np.add.reduce(spltanc) == 0:
                break
            count += 1
            anc1 = spltanc[0]
            curp1 = 0.0
            for traitprob_i in range(len(ch1_tr)):
                traitprob = ch1_tr[traitprob_i]
                if traitprob == 0.0:
                    continue
                curp1 += p1[anc1][traitprob_i] * traitprob

            anc2 = spltanc[1]
            curp2 = 0.0
            for traitprob_i in range(len(ch2_tr)):
                traitprob = ch2_tr[traitprob_i]
                if traitprob == 0.0:
                    continue
                curp2 += p2[anc2][traitprob_i] * traitprob
            anclike += (curp1 * curp2) * weight
        n.timeslice_lv[0][chari][ti] = anclike
    if n.parent != None:
        n.scaling_factors[0][chari] = max(n.timeslice_lv[0][chari]) # update scaling factor vector
        try:
            n.timeslice_lv[0][chari] = max_scale_timeslice(n.timeslice_lv[0][chari])
        except:
            print("ERROR SCALING LIKELIHOOD VECTOR in `mfc.split_loglike_single_trait_marg()`")
            sys.exit()

cdef double split_like_marg(node.Node n, qmat.Qmat qmats, long[:] ss):
    cdef double nodelike, charll, dt
    cdef int chari, cur_k #, ti, tj, anc1, anc2, traitprob_i, nscenarios
    cdef node.Node ch
    cdef double[:, :, :] pmats0, pmats1
    cdef double[:, :] p0, p1

    if len(n.children) != 2:
        print("hypothetical ancestors can only have two children in the model.")
        sys.exit() 

    dt = get_child_dt(n.children[0])
    pmats1 = qmats.calc_p_mats(dt)
    dt = get_child_dt(n.children[1])
    pmats2 = qmats.calc_p_mats(dt)

    for chari in range(0,len(ss)):
        cur_k = ss[chari]
        if cur_k == 1:
            continue

        p1 = pmats1[cur_k - 2]
        p2 = pmats2[cur_k - 2]
        split_loglike_single_trait_marg(n, p1, p2, cur_k, chari)

def max_scale_lv(node.Node n):
    cdef int i, j, k
    cdef double[:,:] cur_lv
    cdef double[:] curtr_lv
    cdef double max_val

    for i in range(len(n.timeslice_lv)):
        cur_lv = n.timeslice_lv[i]
        for j in range(1,len(cur_lv)):
            curtr_lv = cur_lv[j]
            max_val = max(curtr_lv)
            for k in range(len(curtr_lv)):
                n.timeslice_lv[i][j][k] = curtr_lv[k] / max_val 
            n.scaling_factors[i][j] = max_val
 
cdef double[:] max_scale_timeslice(double[:] anc_marg):
    cdef double[:] scaled_anc_marg
    cdef int i

    #print("IN MAX SCSALE TIMESLICE", list(anc_marg))
    scaled_anc_marg = np.zeros(len(anc_marg))
    for i in range(len(anc_marg)):
        scaled_anc_marg[i] = anc_marg[i] / max(anc_marg)

    return scaled_anc_marg

cdef bint check_mis(double[:] tr_vec):
    cdef bint mis
    cdef int i

    mis = False
    if np.add.reduce(tr_vec) == 0.0:
        mis = True
        return mis
    
    for i in range(len(tr_vec)):
        if tr_vec[i] < 1.0 and tr_vec[i] > 0.0:
            mis = True
            return mis
    return mis

cdef double[:] missing_trait_vec(int cur_k):
    cdef double[:] tr = np.zeros(128)
    cdef int j, nstate

    nstate = int(2) ** cur_k
    for j in range(1, len(tr)):
        if j == nstate:
            break
        tr[j] = 1.0
    return tr

"""
def calc_midpoint_ll(node.Node n, qmat.Qmat qmats, double dt, int chari, int cur_k):
    cdef double[:,:] pmat
    cdef double[:] last_like#, all_marg
    cdef int i, j
    cdef double tr_prob, last_like_val, cond_prob, marg_prob 
    pmat = qmats.calc_single_p_mat(dt, cur_k)
    if n.midpoint_lv_index == 0:
        last_like = missing_trait_vec(cur_k)
    else:
        last_like = n.timeslice_lv[n.midpoint_lv_index-1][chari]

    for i in range(len(n.disc_traits[chari])):
        tr_prob = n.disc_traits[chari][i]
        if tr_prob == 0.0:
            continue
        marg_prob = 0.0
        for j in range(len(last_like)):
            last_like_val = last_like[j]
            if last_like_val == 0.0:
                continue
            cond_prob = pmat[i][j]
            cond_prob *= last_like_val
            marg_prob += cond_prob
        n.timeslice_lv[n.midpoint_lv_index][chari][i] = marg_prob
    n.scaling_factors[n.midpoint_lv_index][chari] = max(n.timeslice_lv[n.midpoint_lv_index][chari]) # update scaling factor vector
    n.timeslice_lv[n.midpoint_lv_index][chari] = max_scale_timeslice(n.timeslice_lv[n.midpoint_lv_index][chari])

cdef budd_loglike_single_trait_marg(node.Node n, qmat.Qmat qmats, int cur_k, int chari): #, double desc_weight):
    cdef double scen_cond_like,scenario_like, ana_cond_like, weight, ana_prob, marg_prob, dt, traitll, stateprob, traitprob, prev_time, allstprob = 0.0
    cdef double[:,:] p1, p2
    cdef double[:] chd_tr, par_tr, ana_tr, anc_marg, miss_tr
    cdef long[:,:] cur_scen
    cdef long[:] inher, ancsts
    cdef double[:,:] lv 
    cdef bint mis, past_mid#, last_ts
    cdef int lv_i, chd_i, i, nscen, j, ancst, startst, k, l, maxstates
    cdef node.Node ch

    cur_scen = bm.get_buddmat(cur_k)
    past_mid = False
    prev_time = n.upper
    miss_tr = missing_trait_vec(cur_k)
    #prev_ls = [prev_time]
    #prev_ch = ["tip"]
    if n.midpoint_lv_index == 0:
            dt = n.midpoint - prev_time
            calc_midpoint_ll(n, qmats, dt, chari, cur_k) # need to compute the likelihood at the midpoint if we've hit the first child beyond the midpoint
            past_mid = True
            prev_time = n.midpoint
            #prev_ls.append(prev_time)
            #prev_ch.append("mid")

    for chd_i in reversed( range( len(n.children) ) ): # start with child farthest from the root
        ch = n.children[chd_i]

        #prev_ch.append(ch.label)
        #prev_ls.append(prev_time)
        if len(ch.children) > 0 and ch.midpoint_lv_index != len(ch.children):
            dt = ch.lower - ch.children[0].lower
        elif len(ch.children) > 0 and ch.midpoint_lv_index == len(ch.children) or len(ch.children) == 0:
            dt = ch.lower - ch.midpoint
        
        if dt < 0.0000001:
            print("SPOT1",n.label, ch.label, dt, "ERROR: ZERO BLEN ENCOUNTERED")
            print(n.label,ch.label,ch.lower,ch.upper,ch.midpoint,ch.children[0].lower,ch.children[0].label)
            sys.exit()
        
        p1 = qmats.calc_single_p_mat(dt, cur_k) # P matrix for segment along the child leading to the first likelihood vector

        dt = ch.lower - prev_time
        p2 = qmats.calc_single_p_mat(dt, cur_k)
        if chd_i == len(n.children) - 1 and past_mid == False: # if we are on most distal desc and that occurs before midpoint
            ana_tr = miss_tr 
        else:
            ana_tr = n.timeslice_lv[ch.parent_lv_index - 1][chari]
 
        if dt < 0.0000001:
            print("SPOT2",n.label, ch.label, dt, "ERROR: ZERO BLEN ENCOUNTERED")
            print(n.label,ch.label,ch.lower,ch.upper,ch.midpoint,ch.children[0].lower,ch.children[0].label)
            sys.exit()

        chd_tr = ch.timeslice_lv[-1][chari] #.disc_traits[chari]            
        #print(n.label,ch.label,list(chd_tr), dt0, dt)
        allstprob = 0.0
        ancsts = np.array([i for i in range(1,len(cur_scen))])

        anc_marg = np.zeros(len(ancsts)+1)
        for ancst in ancsts:
            inher = cur_scen[ancst]
            nscen = count_nscenarios(inher)
            weight = 1.0 / ( float(nscen) * float(len(cur_scen)-1) )
            marg_prob = 0.0
            ana_cond_like = 0.0
            for l in range(len(ana_tr)):
                ana_prob = ana_tr[l]
                if ana_prob == 0.0:
                    continue
                ana_cond_like += p2[ancst][l] * ana_prob

            for j in range(len(inher)): # calc cond like for each allowed inheritance scen
                startst = inher[j]
                if startst == 0:
                    break

                scen_cond_like = 0.0
                for k in range(len(chd_tr)):
                    traitprob = chd_tr[k]
                    if k == int(2) ** cur_k: # NEED TO FIX
                        break
                    if traitprob == 0.0:
                        continue
                    scen_cond_like += p1[startst][k] * traitprob
                scenario_like = ana_cond_like * scen_cond_like 
                
                scenario_like *= weight
                marg_prob += scenario_like
            anc_marg[ancst] = marg_prob 
        n.scaling_factors[ch.parent_lv_index][chari] = max(anc_marg) # update scaling factor vector
        
        try:
            anc_marg = max_scale_timeslice(anc_marg)
        except:
            print("ERROR RESCALING LIKELIHOOD VECTOR `budd_loglike_single_trait_marg()")
            sys.exit()
    
        for j in range(len(anc_marg)):
            n.timeslice_lv[ch.parent_lv_index][chari][j] = anc_marg[j] 

        prev_time = ch.lower
        if ch.parent_lv_index == n.midpoint_lv_index - 1:
            dt = n.midpoint - prev_time
            calc_midpoint_ll(n, qmats, dt, chari, cur_k) # need to compute the likelihood at the midpoint if we've hit the first child beyond the midpoint
            past_mid = True
            prev_time = n.midpoint
            #prev_ls.append(prev_time)
            #prev_ch.append("mid_post")




#def update_scaling_factors(double[:] anc_marg, ):
#    for 

cdef budd_like_marg(node.Node n, qmat.Qmat qmats, long[:] ss):
    cdef int chari, cur_k #, chd_i, i ,j,nscen, ancst, k, maxstates
    cdef double nodell, traitll#, desc_weight

    #desc_weight = 1.0 / float(len(n.children))
    
    nodell = 0.0
    for chari in range(1, len(ss)):
        cur_k = ss[chari]
        
        if cur_k == 1:
            continue

        budd_loglike_single_trait_marg(n, qmats, cur_k, chari)#, desc_weight)
        
"""
def calc_midpoint_ll(node.Node n, double[:, :, :] pmats, long[:] ss):
    cdef int cur_k
    cdef double [:, :] pmat

    for chari in range(0,len(n.timeslice_lv[0])):
        cur_k = ss[chari]
        if cur_k == 1:
            continue
        pmat = pmats[cur_k - 2] 
        calc_midpoint_ll_single_trait(n, pmat, cur_k, chari)

def calc_midpoint_ll_single_trait(node.Node n, double[:, :] pmat, int cur_k, int chari):
    #cdef double[:,:] pmat
    cdef double[:] last_like#, all_marg
    cdef int i, j #, chari, curk
    cdef double tr_prob, last_like_val, cond_prob, marg_prob 
    if n.midpoint_lv_index == 0:
        last_like = missing_trait_vec(cur_k)
    else:
        last_like = n.timeslice_lv[n.midpoint_lv_index-1][chari]
    
    for i in range(len(n.disc_traits[chari])):
        tr_prob = n.disc_traits[chari][i]
        if tr_prob == 0.0:
            continue
        marg_prob = 0.0
        for j in range(len(last_like)):
            last_like_val = last_like[j]
            if last_like_val == 0.0:
                continue
            cond_prob = pmat[i][j]
            cond_prob *= last_like_val
            marg_prob += cond_prob

        n.timeslice_lv[n.midpoint_lv_index][chari][i] = marg_prob
    n.scaling_factors[n.midpoint_lv_index][chari] = max(n.timeslice_lv[n.midpoint_lv_index][chari]) # update scaling factor vector
    n.timeslice_lv[n.midpoint_lv_index][chari] = max_scale_timeslice(n.timeslice_lv[n.midpoint_lv_index][chari])

cdef budd_loglike_single_trait_marg(node.Node n, node.Node ch, double[:,:] p1, double[:,:] p2, int cur_k, int chari): #, double desc_weight):
    cdef double scen_cond_like,scenario_like, ana_cond_like, weight, ana_prob, marg_prob, dt, traitll, stateprob, traitprob, prev_time, allstprob = 0.0
    #cdef double[:,:] p1, p2
    cdef double[:] chd_tr, par_tr, ana_tr, anc_marg, miss_tr
    cdef long[:,:] cur_scen
    cdef long[:] inher, ancsts
    cdef int lv_i, chd_i, i, nscen, j, ancst, startst, k, l, maxstates
    
    if ch.parent_lv_index == 0:
        ana_tr = missing_trait_vec(cur_k)
    else:
        ana_tr = n.timeslice_lv[ch.parent_lv_index - 1][chari]

    chd_tr = ch.timeslice_lv[-1][chari] 
    cur_scen = bm.get_buddmat(cur_k)
 
    allstprob = 0.0
    ancsts = np.array([i for i in range(1,len(cur_scen))])

    anc_marg = np.zeros(len(ancsts)+1)
    for ancst in ancsts:
        inher = cur_scen[ancst]
        nscen = count_nscenarios(inher)
        weight = 1.0 / ( float(nscen) * float(len(cur_scen)-1) )
        marg_prob = 0.0
        ana_cond_like = 0.0
        for l in range(len(ana_tr)):
            ana_prob = ana_tr[l]
            if ana_prob == 0.0:
                continue
            ana_cond_like += p2[ancst][l] * ana_prob

        for j in range(len(inher)): # calc cond like for each allowed inheritance scen
            startst = inher[j]
            if startst == 0:
                break

            scen_cond_like = 0.0
            for k in range(len(chd_tr)):
                traitprob = chd_tr[k]
                #print("TRAITPROB",traitprob)

                #print("PMATVAL",p1[startst][k])
                if k == int(2) ** cur_k: # NEED TO FIX
                    break
                if traitprob == 0.0:
                    continue
                scen_cond_like += p1[startst][k] * traitprob
            #print("SUBLIKES:", ana_cond_like, scen_cond_like)
            scenario_like = ana_cond_like * scen_cond_like 
            
            scenario_like *= weight
            marg_prob += scenario_like
        anc_marg[ancst] = marg_prob 
    n.scaling_factors[ch.parent_lv_index][chari] = max(anc_marg) # update scaling factor vector
    
    try:
        anc_marg = max_scale_timeslice(anc_marg)
    except:
        print("ERROR RESCALING LIKELIHOOD VECTOR `budd_loglike_single_trait_marg()")
        sys.exit()

    for j in range(len(anc_marg)):
        n.timeslice_lv[ch.parent_lv_index][chari][j] = anc_marg[j] 


cdef budd_like_marg(node.Node n, qmat.Qmat qmats, long[:] ss):
    cdef int chari, cur_k, chd_i #, i ,j,nscen, ancst, k, maxstates
    cdef double prev_time, dt 
    cdef node.Node ch
    cdef bint past_mid
    cdef double[:, :, :] pmats1, pmats2
    cdef double[:, :] p1, p2
    #desc_weight = 1.0 / float(len(n.children))

    past_mid = False
    prev_time = n.upper
    
    if n.midpoint_lv_index == 0:
            dt = n.midpoint - prev_time
            pmats1 = qmats.calc_p_mats(dt)
            calc_midpoint_ll(n, pmats1, ss) # need to compute the likelihood at the midpoint if we've hit the first child beyond the midpoint
            """
            if n.label == "Copelemur_australotutus":
                #print(n.label,"disc traits")
                #print(list(n.disc_traits[chari]))
                print(n.label, "MIDPOINT",n.midpoint_lv_index)
                print("CURMID:",n.midpoint)
                print(n.children[0].label, "DESCLOWER:",n.children[0].lower)
                print(n.children[0].label, n.children[0].parent_lv_index)
                print(list(n.timeslice_lv[n.midpoint_lv_index][chari]))
                exit()
            """
            past_mid = True
            prev_time = n.midpoint

    for chd_i in reversed( range( len(n.children) ) ): # start with child farthest from the root
        ch = n.children[chd_i]
        
        ## NEED TO COME BACK AND ADD get_child_dt() after debugging
        #if len(ch.children) > 0 and ch.midpoint_lv_index != len(ch.children):
        #    dt = ch.lower - ch.children[0].lower
        #elif len(ch.children) > 0 and ch.midpoint_lv_index == len(ch.children) or len(ch.children) == 0:
        #    dt = ch.lower - ch.midpoint

        dt = get_child_dt(ch)
        pmats1 = qmats.calc_p_mats(dt)
        dt = ch.lower - prev_time
        pmats2 = qmats.calc_p_mats(dt)
        
        for chari in range(0, len(ss)):
            cur_k = ss[chari]
            if cur_k == 1:
                continue
            
            p1 = pmats1[cur_k - 2]
            p2 = pmats2[cur_k - 2]
            
            budd_loglike_single_trait_marg(n, ch, p1, p2, cur_k, chari)#, desc_weight)
 
        prev_time = ch.lower
        if ch.parent_lv_index == n.midpoint_lv_index - 1:
            dt = n.midpoint - prev_time
            pmats1 = qmats.calc_p_mats(dt)
            calc_midpoint_ll(n, pmats1, ss) # need to compute the likelihood at the midpoint if we've hit the first child beyond the midpoint

            past_mid = True
            prev_time = n.midpoint

cdef int count_nscenarios(long[:] inher):
    cdef int nscen, j
    nscen = 0
    for j in range(len(inher)):
        if inher[j] == 0:
            break
        nscen += 1
    return nscen

# this loops over all possible inherited states after speciation and computes the probability 
# of each possible anagenetic transition and returns their sum 
cdef double p_over_desc(long[:] inher, double[:] chd_tr, double[:,:] p1, double weight):
    cdef int j, k, ancst
    cdef double marg_prob,stateprob,traitprob
    marg_prob = 0.0
    print(list(inher))
    for j in range(len(inher)):
        stateprob = 0.0
        ancst = inher[j]
        if ancst == 0:
            break
        for k in range(len(chd_tr)):
            traitprob = chd_tr[k]
            if traitprob == 0.0:
                continue
                #print(ancst,k,traitprob)

            stateprob += p1[ancst][k] * traitprob
            stateprob *= weight
        marg_prob += stateprob
    return marg_prob

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
        nscen = count_nscenarios(inher)
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

def calc_logsum_scaling_factors(node.Node tree, long[:] ss):
    sum_log_sf = np.zeros(len(ss)) 
    for n in tree.iternodes(1):
        if n == tree and n.istip == False: 
            continue
        if len(n.children) == 0:
            continue
        for i in range(len(n.scaling_factors)):
            for j in range(1,len(n.scaling_factors[i])):
                if ss[j] > 1:
                    sum_log_sf[j] += np.log(n.scaling_factors[i][j])
    return sum_log_sf

def mfc2_treell(node.Node tree, qmat.Qmat qmats, long[:] ss, bint asc = True):
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
            split_like_marg(n, qmats, ss)
        elif n.istip:
            budd_like_marg(n, qmats, ss)

    if tree.istip == False:
        root_marg_likes = tree.timeslice_lv[-1]
    else:
        root_marg_likes = calc_bud_root_ll(tree, qmats, ss)

    sum_log_sf = calc_logsum_scaling_factors(tree, ss) 

    treell = 0.0
    for i in range(1, len(root_marg_likes)):
        #for j in i:
        if ss[i] == 1: # invariant traits do not contribute to the tree likelihood bc we use ascertainment bias correction
            continue
        plikes = root_marg_likes[i]
        sum_plikes = sum(plikes)
        sublike = 0.0
        #print("PLIKES",list(plikes)[0:8])
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
    #print("MFC2 TREELL:",treell)
    if asc == True:
        invarll = calc_invar_ll_marg(tree,qmats)
        treell = treell - np.log(1.0-np.exp(invarll))
    
    #print("CORRECTED MFC2 TREELL:",treell)
    return treell

def calc_ASR_down_budd_node(node.Node n, qmat.Qmat qmats, long[:] ss):#, double[:, :] prev_marg):
    cdef node.Node ch
    cdef double[:,:] p1
    cdef double dt, last_time
    cdef int i, cur_k, ch_i
    last_time = n.lower
    past_mid = False
    mid_first = False
    if n.midpoint_lv_index == len(n.timeslice_lv) - 1: # check if midpoint is oldest timeslice along branch
            n.timeslice_lv[n.midpoint_lv_index] = n.disc_traits # no need to calculate midpoint ASRs since we condition on taxon's observed state
            past_mid = True
            last_time = n.midpoint
            mid_first = True

    for ch_i in range(len(n.children)):
        ch = n.children[ch_i]
        if ch_i > 0 or mid_first == True:
            dt = last_time - ch.lower
            prev_marg = n.timeslice_lv[ch.parent_lv_index + 1] # previous marg vec is the index one higher than cur (moving forward in time)
            budd_anagenesis_one_step(n, ch, dt, qmats, ss, prev_marg)

        dt = get_child_dt(ch)
        budd_cladogenesis_one_step(ch, qmats, ss, dt)
        last_time = ch.lower
        if ch.parent_lv_index == n.midpoint_lv_index + 1 and past_mid == False:
            n.timeslice_lv[n.midpoint_lv_index] = n.disc_traits # no need to calculate midpoint ASRs since we condition on taxon's observed state
            past_mid = True
            last_time = n.midpoint

def budd_cladogenesis_one_step(node.Node ch, qmat.Qmat qmats, long[:] ss, double dt):
    cdef int i, j, k, cur_k, nscen, ancst, traitprob_i, anc1
    cdef double cumprob, curprob, weight, traitprob
    cdef double[:, :] p1, desc_lv, anc_marg
    cdef double[:, :, :] pmats
    cdef long[:] ancsts

    if ch.parent == None:
        print("ERROR: cannot use cladogenetic process on branch with no parent!!")
        sys.exit()

    dt = get_child_dt(ch)
    pmats = qmats.calc_p_mats(dt)
    
    desc_lv = ch.timeslice_lv[-1]
    anc_marg = ch.parent.timeslice_lv[ch.parent_lv_index]
    
    for i in range(1, len(ss)):
        forward_probs1 = np.zeros(len(desc_lv[i])) 
        cur_k = ss[i]
        if cur_k == 1:
            continue
        #print(cur_k, list(anc_marg[i][0:10]))
        cur_scen = bm.get_buddmat(cur_k)
        p1 = pmats[cur_k - 2] 
        for traitprob_i in range(len(desc_lv)):
            traitprob = desc_lv[i][traitprob_i]
            if traitprob == 0.0:
                continue
            for ti in range(1,len(cur_scen)):
                anc_prob = anc_marg[i][ti] 
                nscen = 0
                inher = cur_scen[ti]
                nscen = count_nscenarios(inher)
                weight = 1.0 / ( float(nscen) * float(len(cur_scen)-1) )
                
                cum_anc_prob = 0.0
                for tj in range(len(cur_scen[ti])):
                    anc1 = cur_scen[ti][tj]
                    if anc1 == 0: #np.add.reduce(spltanc) == 0:
                        break
                    curp1 = p1[anc1][traitprob_i]
                    cum_anc_prob += curp1 * weight
                cum_anc_prob *= anc_prob
                forward_probs1[traitprob_i] += cum_anc_prob
        #print(ch.parent.label, ch.label, list(forward_probs1)[0:5])
        forward_probs1 = forward_probs1 * desc_lv[i]
        cumprob = sum(forward_probs1)
        for j in range(len(forward_probs1)):
            forward_probs1[j] = forward_probs1[j] / cumprob
            ch.timeslice_lv[-1][i][j] = forward_probs1[j] 
        #print(ch.parent.label, ch.label, list(forward_probs1)[0:5])
        #ch.timeslice_lv[-1][i] = forward_probs1

def budd_anagenesis_one_step(node.Node n, node.Node ch, double dt, qmat.Qmat qmats, long[:] ss, double[:,:] prev_marg):
    cdef int i, j, k, cur_k
    cdef double cumprob, curprob, sum_all
    cdef double[:, :] p1, cur_lv#, prev_marg
    cdef double[:, :, :] pmats

    pmats = qmats.calc_p_mats(dt)
    if ch.parent_lv_index == len(n.timeslice_lv) - 1:
        print("ERROR: cannot calc ASRs based on anagenetic process for first timeslice along branch. need to do cladogenetic process")
        sys.exit()

    cur_lv = n.timeslice_lv[ch.parent_lv_index]
    #prev_marg = n.timeslice_lv[ch.parent_lv_index + 1] # previous marg vec is the index one higher than cur (moving forward in time)
    for i in range(len(ss)):
        cur_k = ss[i] 
        p1 =  pmats[cur_k - 2] 

        forward_probs1 = np.zeros(len(cur_lv[i])) 
        for j in range(len(cur_lv[i])):
            if cur_lv[i][j] == 0.0:
                continue

            cumprob = 0.0 
            for k in range(len(prev_marg[i])):
                if prev_marg[i][k] == 0.0:
                    continue
                curprob = p1[k][j]
                cumprob += curprob * prev_marg[i][k]
            forward_probs1[j] = cumprob

        forward_probs1 = forward_probs1 * cur_lv[i]
        sum_all = sum(forward_probs1)

        for j in range(len(forward_probs1)):
            forward_probs1[j] = forward_probs1[j] / sum_all

        for j in range(len(forward_probs1)):
            n.timeslice_lv[ch.parent_lv_index][i][j] = forward_probs1[j]

def compute_mfc2_ASRs(node.Node tree, qmat.Qmat qmats, long[:] ss):
    cdef node.Node n, ch
    cdef double[:,:] root_marg_likes, curlike, p1, prev_marg
    cdef double[:] plikes, sum_log_sf, curvec
    cdef double lasttime, dt, cond, anc, marg
    cdef int i, j, count, cur_k

    mfc2_treell(tree,qmats,ss,True) # propogate backward conditional likelihoods
    if tree.istip == False:
        if len(tree.timeslice_lv) > 1:
            print("root is bifurcating but there is more than one likelihood vector. something is wrong here.")
            print("error in `mfc.compute_mfc2_ASRs()`")
            sys.exit(0)

        root_marg_likes = tree.timeslice_lv[0]
    else:
        root_marg_likes = calc_bud_root_ll(tree, qmats, ss)

    root_marg_likes = normalize_root_marg_likes(root_marg_likes, ss, "flat")
    
    if tree.istip == False:
        tree.timeslice_lv[0] = root_marg_likes 
        splitting_forward_probs(tree, qmats, ss)
    else: # TODO: need to come back and fix case where root is not bifurcating but is a sampled anc
        print("CAUTION: CASE WHERE ROOT IS SAMPLED ANC WITH BUDD DESCENDANTS IS NOT WORKING RIGHT YET")
        #if tree.children[0].lower > tree.midpoint:
            #calc_ASR_down_budd_node(tree, qmats, ss, root_marg_likes)

    for n in tree.iternodes(0):
        if n == tree:
            continue

        if n.istip == False:
            splitting_forward_probs(n, qmats, ss)
        else:
            #prev_marg = n.parent.timeslice_lv[n.parent_lv_index]
            calc_ASR_down_budd_node(n, qmats, ss)


def splitting_forward_probs(node.Node n, qmat.Qmat qmats,  long[:] ss):
    cdef int i
    cdef double dt
    cdef double[:,:,:] pmats1, pmats2
    cdef node.Node ch

    ch = n.children[0]
    dt = get_child_dt(ch)
    pmats1 = qmats.calc_p_mats(dt)

    ch = n.children[1]
    dt = get_child_dt(ch)
    pmats2 = qmats.calc_p_mats(dt)


    for i in range(1,len(n.timeslice_lv[0])): # i is a trait index
        if ss[i] == 1:
            continue

        splitting_forward_probs_single_trait(n, pmats1, pmats2, ss[i], i)


def get_child_dt(node.Node ch):
    cdef double dt

    if ch.istip:
        if len(ch.children) > 0 and ch.midpoint_lv_index != len(ch.children):
            dt = ch.lower - ch.children[0].lower
        elif len(ch.children) > 0 and ch.midpoint_lv_index == len(ch.children) or len(ch.children) == 0:
            dt = ch.lower - ch.midpoint
    else:
        dt = ch.length

    if dt < 0.0000001:
        print(ch.label, dt)
        print("ERROR: found zero length in `mfc.get_child_dt()`")
        sys.exit()
    return dt

def splitting_forward_probs_single_trait(node.Node n, double[:,:,:] pmats1, double[:,:,:] pmats2, int cur_k, int chari):
    cdef double traitprob, curp1, curp2, weight, dt, anc_prob, cum_anc_prob
    cdef long[:] spltanc 
    cdef int ti, tj, anc1, anc2, traitprob_i, nscenarios
    cdef double[:,:] p1, p2
    cdef double[:] ch1_tr,ch2_tr
    cdef long[:,:,:] cur_scen
    cdef node.Node ch

    if n.istip or len(n.children) != 2:
        print("trying to calc forward splitting cladogenetic probs on non-splitting node")
        print("problem in `splitting_forward_probs_single_trait()`")
        sys.exit()

    cur_scen = sm.get_spltmat(cur_k)
    
    ch = n.children[0]
    ch1_tr = ch.timeslice_lv[-1][chari]
    p1 = pmats1[cur_k - 2]
 
    ch = n.children[1]
    ch2_tr = ch.timeslice_lv[-1][chari]
    p2 = pmats2[cur_k - 2]

    forward_probs1 = np.zeros(len(ch1_tr)) 
    forward_probs2 = np.zeros(len(ch2_tr)) 
    for traitprob_i in range(len(ch1_tr)):
        traitprob = ch1_tr[traitprob_i]
        for ti in range(len(cur_scen)):  # each ti is a state present in the ancestor BEFORE cladogenesis
            if ti == 0:
                continue
            anc_prob = n.timeslice_lv[0][chari][ti]
            nscenarios = 0
            for tj in range(len(cur_scen[ti])):
                spltanc = cur_scen[ti][tj]
                if spltanc[0] == 0: #np.add.reduce(spltanc) == 0:
                    break
                nscenarios += 1


            cum_anc_prob = 0.0
            weight = 1.0 / (float(nscenarios) * float(len(cur_scen)-1))

            for tj in range(len(cur_scen[ti])):
                spltanc = cur_scen[ti][tj]
                if spltanc[0] == 0: #np.add.reduce(spltanc) == 0:
                    break
                anc1 = spltanc[0]
                
                if traitprob == 0.0:
                    continue
                
                curp1 = p1[anc1][traitprob_i] #* traitprob
                cum_anc_prob += curp1 * weight
            cum_anc_prob *= n.timeslice_lv[0][chari][ti]
            forward_probs1[traitprob_i] += cum_anc_prob

    for traitprob_i in range(len(ch2_tr)):
        traitprob = ch2_tr[traitprob_i]
        if traitprob == 0.0:
            continue

        for ti in range(len(cur_scen)):  # each ti is a state present in the ancestor BEFORE cladogenesis
            anc_prob = n.timeslice_lv[0][chari][ti]
            if ti == 0:
                continue
            nscenarios = 0
            for tj in range(len(cur_scen[ti])):
                spltanc = cur_scen[ti][tj]
                if spltanc[0] == 0: #np.add.reduce(spltanc) == 0:
                    break
                nscenarios += 1

            cum_anc_prob = 0.0
            weight = 1.0 / (float(nscenarios) * float(len(cur_scen)-1))

            for tj in range(len(cur_scen[ti])):
                spltanc = cur_scen[ti][tj]
                if spltanc[1] == 0: #np.add.reduce(spltanc) == 0:
                    break
             
                anc2 = spltanc[1]
                curp2 = p2[anc2][traitprob_i] #* traitprob
                cum_anc_prob += curp2 * weight
            cum_anc_prob *= anc_prob 
            forward_probs2[traitprob_i] += cum_anc_prob

    forward_probs1 = forward_probs1 * ch1_tr
    forward_probs2 = forward_probs2 * ch2_tr
    curp1 = sum(forward_probs1)
    curp2 = sum(forward_probs2) # just reusing the variable curp2 to store sum of each vec to normalize 
    for ti in range(len(forward_probs1)):
        forward_probs1[ti] = forward_probs1[ti] / curp1
        forward_probs2[ti] = forward_probs2[ti] / curp2

    n.children[0].timeslice_lv[-1][chari] = forward_probs1
    n.children[1].timeslice_lv[-1][chari] = forward_probs2

def normalize_root_marg_likes(double[:,:] root_marg_likes, long[:] ss, prior = "flat"):
    cdef int i, j, count
    cdef double flat_prior, sumvec
    cdef double[:] curvec
    for i in range(len(root_marg_likes)):
        if ss[i] == 1:
            continue

        curvec = root_marg_likes[i]
        sumvec = sum(curvec) 
        #print("BEFOREBEFORE",list(curvec))
        if prior == "flat":
            count = 0
            for j in range(len(curvec)):
                if curvec[j] > 0.0:
                    count+=1
            
            flat_prior = 1.0 / float(count)
            for j in range(len(curvec)):
                curvec[j] = curvec[j] * flat_prior


        elif prior == "estimate":
            for j in range(len(curvec)):
                plike = curvec[j]
                curvec[j] = plike * (plike / sumvec) #flat_prior 


        sumvec = sum(curvec) 
        #print("BEFORE",list(curvec))
        for j in range(len(curvec)):
            if curvec[j] == 0.0:
                continue
            #print(curvec[j], curvec[j] / sumvec)
            curvec[j] = curvec[j] / sumvec
        #if i == 1:
        #    print("AFTER",list(curvec), sum(curvec))
        root_marg_likes[i] = curvec
    return root_marg_likes

def index_of_max(double[:] la):
    cdef int i
    cdef double highest = -10000000000.
    cdef int highest_ind = -1
    for i in range(len(la)):
        if la[i] > highest:
            highest = la[i]
            highest_ind = i
    return highest_ind

def log_sum(double[:] la):
    cdef int i, idx = index_of_max(la)
    cdef double sum_exp = 0.0

    for i in range(len(la)):
        if i == idx:
            continue
        sum_exp += np.exp(la[i] - la[idx])

    return la[idx] + np.log1p(sum_exp)

def compute_anagenetic_marg(double[:,:] pmat, double[:] anc_marg, double[:, :] all_lv):
    cdef double[:] cur_lv, logmarg_vec 
    cdef int i, j, k
    cdef double cond, anc, logmarg, logcond, sum_logmarg, marg, log_frac

    for i in range(1,len(all_lv)): # loop over trait indicies
        cur_lv = all_lv[i]
        logmarg_vec = np.zeros(len(cur_lv))
        for j in range(len(cur_lv)): # loop over individual state partial likelihoods for each trait index
            cond = cur_lv[j] 
            if cond == 0.0:
                continue

            logcond = np.log(cond) 
            
            for k in range(len(anc_marg)):
                if anc_marg[k] == 0.0:
                    continue
                anc += pmat[k][j] * anc_marg[k]
            logmarg = np.log(anc) + logcond
            logmarg_vec[j] = logmarg
        
        sum_logmarg = log_sum(logmarg_vec)
        for j in range(len(logmarg_vec)):
            if logmarg_vec[j] == 0.0:
                continue
            log_frac = logmarg_vec[j] - sum_logmarg
            marg = np.exp(log_frac)
            all_lv[i][j] = marg
    return all_lv


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
    
def evaluate_m_l2(double[:] params, node.Node tree, qmat.Qmat qmats,long[:] ss):
    cdef node.Node n
    cdef double treell
    cdef int maxstates = max(ss)

    if params[0] < 0.0001 or params[1] < 0.0001:
        return 10000000000

    qmats.update_all_qmats(params[0],params[1]) 

    #for n in tree.iternodes():
    #    n.update_pmat(qmats, maxstates)


    #print("ABOUT TO SORT")
    tree_utils.sort_children_by_age(tree)
    #print("SORTED")
    #sys.exit()
    treell = mfc2_treell(tree, qmats, ss)
    return -treell



"""
def evaluate_m_l(double[:] params, node.Node tree, qmat.Qmat qmats,long[:] ss):
    cdef node.Node n
    cdef double treell
    cdef int maxstates = max(ss)

    if params[0] < 0.0001 or params[1] < 0.0001:
        return 10000000000

    qmats.update_all_qmats(params[0],params[1]) 

    for n in tree.iternodes():
        n.update_pmat(qmats, maxstates)

    treell = mfc_treell(tree, ss)
    return -treell
"""
"""
def create_q_matrix(int nstates, double m, double l):
    cdef double[:,:] qmat
    cdef long[:,:] smap
    cdef int matsize
    smap = smaps.get_smap(nstates)
    matsize = len(smap) #2**nstates
    qmat = np.zeros((matsize,matsize),dtype=np.double)
    update_mut_loss(qmat,smap,m,l)
    #for row in qmat:
    #    print(list(row))                    
    return qmat, smap

def update_mut_loss(double[:,:] ratemat, long[:,:] smap, double m, double l):
    cdef int i, j, ii, ndiff, nrow, nsti, nstj, nstates 
    nstates = len(smap[0])
    nrow = len(ratemat)
    for i in range(nrow):
        if i == 0:
            continue
        for j in range(nrow):
            if i != j:
                ndiff = 0
                for ii in range(nstates):
                    if smap[i][ii] != smap[j][ii]:
                        ndiff += 1
                nsti = np.add.reduce(smap[i])
                nstj = np.add.reduce(smap[j])
                if nstj - nsti < 0:
                    if ndiff == 1:
                        ratemat[i][j] = l
                elif nstj - nsti > 0:
                    if ndiff == 1:
                        ratemat[i][j] = m
    for i in range(nrow):
        ratemat[i][i] = -np.add.reduce(ratemat[i])
"""

"""
cdef double budd_like_single_inheritance1(long[:] inher, double[:] chd_tr, double[:] ana_tr, double[:,:] p1, double[:,:] p2, double weight):
    cdef int j, k, ancst
    cdef double marg_prob,stateprob,traitprob
    marg_prob = 0.0
    print(list(inher))
    for j in range(len(inher)):
        stateprob = 0.0
        ancst = inher[j]
        if ancst == 0:
            break
        for k in range(len(chd_tr)):
            traitprob = chd_tr[k]
            if traitprob == 0.0:
                continue
                #print(ancst,k,traitprob)

            stateprob += p1[ancst][k] * traitprob
            print("PROB FROM ",ancst,"TO",k,stateprob)
            stateprob *= weight
            print("WEIGHTED PROB FROM ",ancst,"TO",k,stateprob)
        marg_prob += stateprob
    return marg_prob

cdef double budd_like_single_inheritance(long[:] inher, double[:] chd_tr, double[:] ana_tr, double[:,:] p1, double[:,:] p2, double weight):
    cdef int j, k, ancst
    cdef double marg_prob,stateprob,traitprob
    marg_prob = 0.0
    print(list(inher))
    for j in range(len(inher)):
        stateprob = 0.0
        ancst = inher[j]
        if ancst == 0:
            break
        for k in range(len(chd_tr)):
            traitprob = chd_tr[k]
            if traitprob == 0.0:
                continue
                #print(ancst,k,traitprob)

            stateprob += p1[ancst][k] * traitprob
            print("PROB FROM ",ancst,"TO",k,stateprob)
            stateprob *= weight
            print("WEIGHTED PROB FROM ",ancst,"TO",k,stateprob)
        marg_prob += stateprob
    return marg_prob


cdef budd_loglike_single_trait_marg1(node.Node n, qmat.Qmat qmats, int cur_k, int chari, double desc_weight):
    cdef double scen_cond_like,scenario_like, ana_cond_like, weight, ana_prob, marg_prob, dt, traitll, stateprob, traitprob, prev_time, allstprob = 0.0
    cdef double[:,:] p1, p2
    cdef double[:] chd_tr, par_tr, ana_tr, anc_marg, miss_tr
    cdef long[:,:] cur_scen
    cdef long[:] inher, ancsts
    cdef double[:,:] lv 
    cdef bint mis, past_mid, last_ts
    cdef int lv_i, chd_i, i, nscen, j, ancst, startst, k, l, maxstates
    cdef node.Node ch

    traitll = 0.0
    cur_scen = bm.get_buddmat(cur_k)
    #for inher in cur_scen:
    #    print(list(inher))
    
    #par_tr = n.disc_traits[chari]
    past_mid = False
    last_ts = True
    prev_time = n.upper
    miss_tr = missing_trait_vec(cur_k)
    for chd_i in reversed( range( len(n.children) ) ): # start with child farthest from the root
        ch = n.children[chd_i]
        if ch.lower > n.midpoint:
            if past_mid == False:
                calc_midpoint_ll(n, qmats, n.midpoint - prev_time, chari, cur_k) # need to compute the likelihood at the midpoint if we've hit the first child beyond the midpoint
                past_mid = True
                prev_time = n.midpoint
                if last_ts:
                    last_ts = False


        # NEED TO TEST THESE ONCE I CAN GET THROUGH FIRST BRANCH
        if len(ch.children) > 0 and ch.midpoint_lv_index != len(ch.children):
            dt = ch.lower - ch.children[-1].lower
        elif len(ch.children) > 0 and ch.midpoint_lv_index == len(ch.children) or len(ch.children) == 0:
            dt = ch.lower - ch.midpoint

        p1 = qmats.calc_single_p_mat(dt, cur_k) # P matrix for segment along the child leading to the first likelihood vector

        if chd_i == len(n.children) - 1 and past_mid == False: # if we are on most distal desc and that occurs before midpoint
            dt = ch.lower - prev_time
            p2 = qmats.calc_single_p_mat(dt, cur_k)
            ana_tr = miss_tr 
        else:
            dt = ch.lower - prev_time
            p2 = qmats.calc_single_p_mat(dt, cur_k)
            ana_tr = n.timeslice_lv[ch.parent_lv_index - 1][chari]
        
        chd_tr = ch.timeslice_lv[-1][chari] #.disc_traits[chari]            
        allstprob = 0.0
        ancsts = np.array([i for i in range(1,len(cur_scen))])
        anc_marg = np.zeros(len(ancsts)+1)
        for ancst in ancsts:
            inher = cur_scen[ancst]
            nscen = count_nscenarios(inher)
            weight = 1.0 / ( float(nscen) * float(len(cur_scen)-1) )
            marg_prob = 0.0
            for j in range(len(inher)): # calc cond like for each allowed inheritance scen
                startst = inher[j]
                if startst == 0:
                    break
                scenario_like = 0.0 # marginal likelihood of one scenario, marginalized over all combinations of anagenetic transitions and budded descendant states
                for k in range(len(chd_tr)):
                    traitprob = chd_tr[k]
                    if k == int(2) ** cur_k: # NEED TO FIX
                        break
                    if traitprob == 0.0:
                        continue
                    scen_cond_like = p1[startst][k] * traitprob
                    for l in range(len(ana_tr)):
                        ana_prob = ana_tr[l]
                        if ana_prob == 0.0:
                            continue
                        ana_cond_like = p2[ancst][l] * ana_prob
                        print("ANA PROB FROM", ancst, "TO",l, ana_cond_like)
                        scenario_like += ana_cond_like * scen_cond_like 

                    #print("PROB FROM ",startst,"TO",k,scen_cond_like)
                #print("UNWEIGHTED SCENARIO LIKE",scenario_like)
                scenario_like *= weight
                #print("WEIGHTED SCENARIO LIKE",scenario_like)
                marg_prob += scenario_like
                anc_marg[ancst] = marg_prob 
                #print("MARG PROB FOR ANC STATE: ", ancst, marg_prob)
        #print("MARGINAL PROBABILITY VECTOR: ",list(anc_marg))
        for j in range(len(anc_marg)):
            n.timeslice_lv[ch.parent_lv_index][chari][j] = anc_marg[j] 

        #print("MARGINAL PROBABILITY VECTOR: ",list(n.timeslice_lv[ch.parent_lv_index][chari]))
    sys.exit()

 
"""

cdef double split_loglike_single_trait(node.Node n, int cur_k, int chari):
    cdef double traitprob, curp1, curp2, anclike, weight, charlike
    cdef long[:] spltanc 
    cdef int ti, tj, anc1, anc2, traitprob_i, nscenarios
    cdef double[:,:] p1, p2
    cdef double[:] ch1_tr,ch2_tr
    cdef long[:,:,:] cur_scen


    cur_scen = sm.get_spltmat(cur_k)
    ch1_tr = n.children[0].disc_traits[chari]
    p1 = n.children[0].pmats[cur_k-2] 
    ch2_tr = n.children[1].disc_traits[chari]
    p2 = n.children[1].pmats[cur_k-2] 
    charlike = 0.0
    for ti in range(len(cur_scen)):
        if ti == 0:
            continue
        anclike = 0.0
        nscenarios = 0
        for tj in range(len(cur_scen[ti])):
            spltanc = cur_scen[ti][tj]
            if spltanc[0] == 0: #np.add.reduce(spltanc) == 0:
                break
            nscenarios += 1

        weight = 1.0 / (float(nscenarios) * float(len(cur_scen)-1))
        for tj in range(len(cur_scen[ti])):
            spltanc = cur_scen[ti][tj]

            if spltanc[0] == 0: #np.add.reduce(spltanc) == 0:
                break

            anc1 = spltanc[0]
            curp1 = 0.0
            for traitprob_i in range(len(ch1_tr)):
                traitprob = ch1_tr[traitprob_i]
                if traitprob == 0.0:
                    continue
                curp1 += p1[anc1][traitprob_i] * traitprob

            anc2 = spltanc[1]
            curp2 = 0.0
            for traitprob_i in range(len(ch2_tr)):
                traitprob = ch2_tr[traitprob_i]
                if traitprob == 0.0:
                    continue
                curp2 += p2[anc2][traitprob_i] * traitprob
                #print(curp1,curp2,curp1*curp2)
            anclike += (curp1 * curp2) * weight
        n.disc_traits[chari][ti] = anclike
        charlike += anclike
    return np.log(charlike)

cdef double split_like(node.Node n, long[:] ss):
    cdef double nodelike, charll
    cdef int chari, cur_k #, ti, tj, anc1, anc2, traitprob_i, nscenarios
 
    if len(n.children) != 2:
        print("hypothetical ancestors can only have two children in the model.")
        sys.exit() 

    nodelike = 0.0
    for chari in range(len(ss)):
        cur_k = ss[chari]
        if cur_k == 1:
            continue
        try:
            charll = split_loglike_single_trait(n, cur_k, chari)
        except:
            print("HERE")
            sys.exit()
        nodelike += charll
        #print(n.label,nodelike)
    return nodelike


cdef double budd_loglike_single_trait(node.Node n, int cur_k, int chari, double desc_weight):
    cdef double weight, stateprob,  traitll,allstprob = 0.0
    cdef double[:,:] p1
    cdef double[:] chd_tr, par_tr
    cdef long[:,:] cur_scen
    cdef long[:] inher
    cdef bint mis
    cdef int chd_i, i ,nscen #,j, ancst, k, maxstates

    traitll = 0.0
    cur_scen = bm.get_buddmat(cur_k)
    par_tr = n.disc_traits[chari]
    mis = check_mis(par_tr)
    for chd_i in range(len(n.children)):
        p1 = n.children[chd_i].pmats[cur_k-2] #calc_p_matrix(cur_q,n.children[0].length)
        chd_tr = n.children[chd_i].disc_traits[chari]            
        
        if mis == True: # parent is missing trait
            allstprob = 0.0
            for i in range(1,len(cur_scen)):
                inher = cur_scen[i]
                nscen = count_nscenarios(inher)
                weight = 1.0 / ( float(nscen) * float(len(cur_scen)-1) )

                stateprob = p_over_desc(inher, chd_tr, p1, weight)
                allstprob += stateprob
        elif mis == False:
            for i in range(len(par_tr)): # i == character state in parent state vector
                if par_tr[i] == 1.0:
                    inher = cur_scen[i]
                    nscen = count_nscenarios(inher)
                    weight = 1.0 / float(nscen)
                    allstprob = p_over_desc(inher, chd_tr, p1, weight)
        traitll += np.log(allstprob)
       
    return traitll


cdef double budd_like(node.Node n, long[:] ss):
    cdef int chari, cur_k #, chd_i, i ,j,nscen, ancst, k, maxstates
    cdef double nodell, traitll, desc_weight

    desc_weight = 1.0 / float(len(n.children))
 
    nodell = 0.0
    for chari in range(len(ss)):
    #for chari in prange(len(ss),nogil=True):
        cur_k = ss[chari]
        if cur_k == 1:
            continue
        traitll = budd_loglike_single_trait(n, cur_k, chari, desc_weight)
        #if n.label == "Micraster_quebrada":
        #    print(n.label,traitll)
        nodell += traitll
    return nodell


def mfc_treell(node.Node tree, long[:] ss, bint asc = True):
    cdef double nodell, treell, asc_treell, invarll
    cdef node.Node n
    treell = 0.0
    for n in tree.iternodes(1):
        if len(n.children) == 0:
            continue
        if n.istip == False:
            nodell = split_like(n, ss)
        elif n.istip:
            nodell = budd_like(n, ss)
        if nodell == 0.0:
            print("error in calculating morph ll")
            sys.exit()
        treell += nodell
        #print(n.label,nodell)

    #print(treell)
    if asc == True:
        invarll = calc_invar_ll(tree)
        #print(invarll)
        #print(1.0-np.exp(invarll))
        treell = treell - np.log(1.0-np.exp(invarll))
    #print(treell)
    #sys.exit()
    return treell


cdef double calc_invar_ll(node.Node tree):
    cdef double desc_weight, nodell, treell
    cdef node.Node n 

    treell = 0.0
    for n in tree.iternodes(1):
        if len(n.children) == 0:
            continue
        if n.istip == False:
            nodell = split_loglike_single_trait(n, 2, 0)
            #print(nodelike)
        elif n.istip:
            desc_weight = 1.0 / float(len(n.children))
            nodell = budd_loglike_single_trait(n, 2, 0, desc_weight)
        treell += nodell
    #print(treell)    
    return treell

