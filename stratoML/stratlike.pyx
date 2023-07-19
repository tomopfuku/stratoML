import node
cimport node
import numpy as np
cimport numpy as np
import bd

cpdef double poisson_loglike(double lam, node.Node tree):
    cdef double tf,tl,f,l,samp_time,unsamp_time,lp_samp,lp_unsamp,lp0,lp1,lp01,brlik,treelik = 0.0
    cdef node.Node n

    if lam < 0.000001:
        return 100000000000

    for n in tree.iternodes():
        if n == tree and n.istip == False:
            continue
        tf = n.lower
        tl = n.upper
        if n.istip:
            f = n.strat[0]
            l = n.strat[1]
            samp_time = f-l
            unsamp_time = (tf-tl) - samp_time
            lp0 = -samp_time * lam
            lp1 = np.log(samp_time * lam) - (samp_time * lam)
            lp01 = np.logaddexp(lp0,lp1)
            lp_samp = np.log1p(-np.exp(lp01)) # 1-p(0)-p(1) log-space 
            if unsamp_time > 0.0:
                lp_unsamp = -unsamp_time * lam
            else:
                lp_unsamp = 0.0
            brlik = lp_samp + lp_unsamp
        #elif nfos == 1: #or nfos == 0:
        #    brlik = math.log(lam)+(-lam*(tf-tl))
            #print math.log(lam*(math.exp(-lam*(tf-tl))))
        else:
            brlik =  -(tf-tl) * lam
        #print(n.label,samp_time,unsamp_time,brlik)
        treelik += brlik
    return treelik

def poisson_neg_ll(double lam, node.Node tree):
    return -poisson_loglike(lam,tree)

def bds_loglike(double p, double q, double r, node.Node tree):
    cdef double small,tf,tl,f,l,range,duration,r_like,p_like,q_like,brlik,treelik = 0.0
    cdef node.Node n

    small = 0.0001
    if r < small or p < small or q < small:
        return -100000000000

    for n in tree.iternodes():
        if n == tree and n.istip == False:
            continue
        tf = n.lower
        tl = n.upper
        duration = tf-tl
        if n.istip:
            f = n.strat[0]
            l = n.strat[1]
            range = f-l
            r_like = bd.calc_prob_range(r,duration,range)
            p_like = bd.prob_n_obs_desc(p,q,r,len(n.children),duration)
            q_like = bd.prob_extinction_t(q,duration)
            brlik = np.log(r_like) + np.log(p_like) + np.log(q_like)
            #print(n.label,n.lower,n.upper,n.strat[0],n.strat[1],r_like,p_like,q_like, duration,range, brlik)

        else:
            #print(p,q,r,duration)
            brlik = np.log(bd.bds_hyp_anc_prob(p,q,r,duration))
            #print(p,q,r,duration)
        treelik += brlik
    return treelik

def bds_neg_ll(double[:] params, node.Node tree):
    cdef double p = params[0]
    cdef double q = params[1]
    cdef double r = params[2]
    cdef double treelik = -bds_loglike(p,q,r,tree) 
    #print([i for i in params],treelik)
    return treelik 

def bds_node_starts_neg_ll(double[:] starts, double p, double q, double r, node.Node tree):
    bad = adjust_node_starts(starts,tree)
    if bad:
        return 1000000000
    ll = bds_loglike(p,q,r,tree)
    return -ll

def adjust_node_starts(double[:] starts, node.Node tree):
    cdef node.Node n, sib
    cdef bint bad = False
    cdef double lower

    for n in tree.iternodes(0):
        if n == tree and n.istip == False:
            continue
        lower = starts[n.index]
        if lower < 0.0:
            return bad
        if n.istip:
            if lower <= n.strat[0]:
                #print("HERE1",n.label,lower,n.strat[0])
                bad = True
                return bad
        if n.parent != None:
            if lower >= n.parent.lower or lower < n.parent.upper:
                #print("HERE2")
                bad = True
                return bad
            if n.parent.istip == False:
                sib = n.get_sib()
                if lower != starts[sib.index]:
                    #print("HERE3")
                    return bad
                n.parent.upper = lower
        
        n.lower = lower
        n.length = n.lower - n.upper




def adjust_node_heights_strat(double[:,:] heights, node.Node tree):
    cdef node.Node n
    cdef bint bad = False
    cdef int ncount = 0
    cdef double[:] curh
    cdef double lower, upper 

    for n in tree.iternodes(0):
        curh = heights[ncount]
        lower = curh[0]
        upper = curh[1]
        if lower <= n.lower or upper >= n.upper or lower < upper:
            bad = True
            return bad
        if n.parent != None:
            if lower >= n.parent.lower:
                bad = True
                return bad

        n.lower = lower
        n.upper = upper
        ncount+=1



