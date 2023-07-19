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
#from scipy.optimize import minimize
#from cython.parallel import prange

@cython.wraparound(False)
@cython.boundscheck(False)

def mat_mult(double[:,:] mat, double m):
    cdef int i,j,nrow
    nrow = len(mat)
    cdef double[:,:] newmat = np.empty((nrow,nrow),dtype = np.double)
    for i in range(nrow):
        for j in range(nrow):
            newmat[i][j] = mat[i][j] * m
    return newmat

def calc_p_matrix(double[:,:] ratemat, double t):
    cdef double[:,:] pmat
    pmat = expm(mat_mult(ratemat,t))
    return pmat

cdef double split_loglike_single_trait(node.Node n, int cur_k, int chari):
    cdef double traitprob, curp1, curp2, anclike, weight
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
        charll = split_loglike_single_trait(n, cur_k, chari)
        nodelike += charll
    return nodelike

cdef bint check_mis(double[:] tr_vec):
    cdef bint mis

    mis = False
    if np.add.reduce(tr_vec) == 0.0:
        mis = True
        return mis
    
    for i in range(len(tr_vec)):
        if tr_vec[i] < 1.0 and tr_vec[i] > 0.0:
            mis = True
            return mis
    return mis

cdef double budd_loglike_single_trait(node.Node n, int cur_k, int chari, double desc_weight):
    cdef double weight, stateprob,  traitll,allstprob = 0.0
    cdef double[:,:] p1
    cdef double[:] chd_tr, par_tr
    cdef long[:,:] cur_scen
    cdef long[:] inher
    cdef bint mis
    #cdef int cur_k = ss[chari]
    cdef int chd_i, i ,nscen #,j, ancst, k, maxstates

    traitll = 0.0
    cur_scen = bm.get_buddmat(cur_k)
    par_tr = n.disc_traits[chari]
    mis = check_mis(par_tr)
    if mis:
        for i in range(len(par_tr)):
            par_tr[i] = 0.0
    for chd_i in range(len(n.children)):
        p1 = n.children[chd_i].pmats[cur_k-2] #calc_p_matrix(cur_q,n.children[0].length)
        chd_tr = n.children[chd_i].disc_traits[chari]            
        if mis == True: # parent is missing trait
            allstprob = 0.0
            for i in range(len(cur_scen)):
                if i == 0:
                    continue
                inher = cur_scen[i]
                nscen = count_nscenarios(inher)
                weight = 1.0 / ( float(nscen) * float(len(cur_scen)-1) )
                stateprob = p_over_desc(inher, chd_tr, p1, weight)
                par_tr[i] += stateprob * desc_weight
                #print(i, stateprob,weight)
                allstprob += stateprob
        elif mis == False:
            for i in range(len(par_tr)): # i == character state in parent state vector
                if par_tr[i] == 1.0:
                    inher = cur_scen[i]
                    nscen = count_nscenarios(inher)
                    weight = 1.0 / float(nscen)
                    allstprob = p_over_desc(inher, chd_tr, p1, weight)
                    #traitlike += stateprob
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
        nodell += traitll
    #print(nodell)
    return nodell

cdef int count_nscenarios(long[:] inher):
    cdef int nscen, j
    nscen = 0
    for j in range(len(inher)):
        if inher[j] == 0:
            break
        nscen += 1
    return nscen

cdef double p_over_desc(long[:] inher, double[:] chd_tr, double[:,:] p1, double weight):
    cdef int j, k, ancst
    cdef double stateprob,traitprob
    stateprob = 0.0
    for j in range(len(inher)):
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
    return stateprob

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

def mfc_treell(node.Node tree, long[:] ss, bint asc = True):
    cdef double nodell, treell, asc_treell
    cdef node.Node n
    treell = 0.0
    for n in tree.iternodes(1):
        if len(n.children) == 0:
            continue
        if n.istip == False:
            nodell = split_like(n, ss)
        elif n.istip:
            nodell = budd_like(n, ss)
        treell += nodell

    if asc == True:
        invarll = calc_invar_ll(tree)
        treell = treell - np.log(1.0-np.exp(invarll))
        
    return treell

def evaluate_m_l(double[:] params, node.Node tree, qmat.Qmat qmats,long[:] ss):
    cdef node.Node n
    cdef double treell
    cdef int maxstates = max(ss)

    if params[0] < 0.000001 or params[1] < 0.000001:
        return 10000000000

    qmats.update_all_qmats(params[0],params[1]) 

    for n in tree.iternodes():
        n.update_pmat(qmats, maxstates)

    treell = mfc_treell(tree, ss)
    return -treell

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

