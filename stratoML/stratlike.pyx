import node
cimport node
import numpy as np
cimport numpy as np
import bd
import sys
from scipy.special import logsumexp

cpdef double poisson_loglike(double lam, node.Node tree):
    cdef double tf,tl,f,l,samp_time,unsamp_time,lp_samp,lp_unsamp,lp0,lp1,lp01,brlik,treelik = 0.0
    cdef node.Node n

    if lam < 0.000001:
        return -100000000000.

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
        #    print math.log(lam*(math.exp(-lam*(tf-tl))))
        else:
            brlik =  -(tf-tl) * lam
        #print(n.label,samp_time,unsamp_time,brlik)
        treelik += brlik
    return treelik

def poisson_neg_ll(double lam, node.Node tree):
    return -poisson_loglike(lam,tree)

def single_range_ll(double gap, double p, double q, double r, double obs_range, int nch):
    cdef double r_like,p_like,q_like, duration, ll
    duration = obs_range + gap
    if gap <= 0.0 or gap > 20.0:
        return -1000000000.0

    r_like = np.log(bd.calc_prob_range(r,duration,obs_range))
    p_like = np.log(bd.prob_n_obs_desc(p,q,r,nch,duration))
    q_like = np.log(bd.prob_extinction_t(q,duration))
    ll = p_like + q_like + r_like
    return -ll 

def single_range_mle_optim(double p, double q, double r, node.Node n):
    obs_range = n.strat[0] - n.strat[1]
    duration = obs_range + 0.01
    #res = minimize(stratlike.single_range_ll,x0=0.01,args=(p,q,r,obs_range,nch),method="BFGS")#,bounds=(((0.0001,5.0))))


def single_range_hpd(double p,double q,double r, node.Node n, double step = 0.01, double quant = .95):
    cdef double duration, obs_range, r_like, p_like, q_like, loglike, total_h, lower,upper,z,lowerd,upperd, summed
    cdef bint foundl = False
    cdef list dist, durations

    #duration = n.lower - n.upper
    obs_range = n.strat[0] - n.strat[1]
    duration = obs_range + step
    dist = []
    durations = []
    for _ in range(50000):
        r_like = np.log(bd.calc_prob_range(r,duration,obs_range))
        p_like = np.log(bd.prob_n_obs_desc(p,q,r,len(n.children),duration))
        q_like = np.log(bd.prob_extinction_t(q,duration))
        loglike = r_like + p_like + q_like
        dist.append(loglike)
        durations.append(duration)
        #loglike = bds_extinct_branch_ll(p,q,r,n)
        #print(n.labejloglike,duration,duration-obs_range)
        #if loglike < np.log(0.000000001):
        if loglike < -300.:
            break
        duration += step

    total_h = logsumexp(dist)
    lowerd = total_h + np.log(.025)
    upperd = total_h + np.log(.975)
    z = dist[0]
    for i in range(1,len(dist)):
        summed = np.logaddexp(z,dist[i])
        if summed >= lowerd and foundl == False:
            lower = durations[i]
            foundl = True
        elif summed >= upperd:
            upper = durations[i]
            break
        z = summed
    #print(n.length,lower,upper)
    return lower,upper


def hyp_anc_hpd(double p,double q,double r, node.Node n, double step = 0.01, double quant = .95):
    cdef double duration, obs_range, r_like, p_like, q_like, loglike, total_h, lower,upper,z,lowerd,upperd, summed
    cdef bint foundl = False
    cdef list dist, durations

    #duration = n.lower - n.upper
    obs_range = n.strat[0] - n.strat[1]
    duration = obs_range + step
    dist = []
    durations = []
    for _ in range(50000):
        loglike = bd.bds_hyp_anc_log_prob(p,q,r,duration)
        dist.append(loglike)
        durations.append(duration)
        if loglike < -300.:
            break
        duration += step

    total_h = logsumexp(dist)
    lowerd = total_h + np.log(.025)
    upperd = total_h + np.log(.975)
    z = dist[0]
    for i in range(1,len(dist)):
        summed = np.logaddexp(z,dist[i])
        if summed >= lowerd and foundl == False:
            lower = durations[i]
            foundl = True
        elif summed >= upperd:
            upper = durations[i]
            break
        z = summed
    #print(n.length,lower,upper)
    return lower,upper


def bds_CIs(double p, double q, double r, node.Node tree):
    cdef double gap, mllen, f, l, hypanc_len, oldest , obs_range,lower,upper,gap_lower,gap_upper
    cdef node.Node n, ch
    cdef dict node_cis = {}

    #hypanc_len = hyp_anc_mle(p,q,r)
    for n in tree.iternodes(order=1):
        if n.istip:
            obs_range = n.strat[0] - n.strat[1]
            lower,upper = single_range_hpd(p,q,r,n)
            gap_lower = (lower-obs_range) / 2.0
            gap_upper = (upper-obs_range) / 2.0
            #print(gap_lower,n.lower,n.strat[0])
            node_cis[n] = (n.strat[0] + gap_lower, n.strat[0] + gap_upper)
        else:
            node_cis[n] = (n.lower,n.lower)
    return node_cis

def single_range_mle(double p,double q,double r, node.Node n):
    cdef double duration, obs_range, r_like, p_like, q_like, loglike, best, best_dur

    #duration = n.lower - n.upper
    obs_range = n.strat[0] - n.strat[1]
    duration = obs_range + 0.01
    best = -1000000000.0
    best_dur = duration
    while True:
        r_like = np.log(bd.calc_prob_range(r,duration,obs_range))
        p_like = np.log(bd.prob_n_obs_desc(p,q,r,len(n.children),duration))
        q_like = np.log(bd.prob_extinction_t(q,duration))
        loglike = r_like + p_like + q_like
        #loglike = bds_extinct_branch_ll(p,q,r,n)
        #print(n.label,loglike,duration,duration-obs_range)
        if loglike > best:
            best = loglike
            best_dur = duration
        elif loglike < best:
            break
        duration += 0.01
    return best_dur

cpdef double hyp_anc_mle(double p, double q, double r):
    cdef double mle
    mle = np.log( (p+q+r) / (q+r) ) / p
    return mle


def hyp_anc_mle_sim(double p, double q, double r):
    cdef double duration, best, loglike, bestdur
    duration = 0.01
    best = -1000000000.0
    bestdur = duration
    while True:
        loglike = bd.bds_hyp_anc_log_prob(p,q,r,duration)
        #print(loglike,duration)
        if loglike > best:
            best = loglike
            bestdur = duration
        elif loglike < best:
            break
        duration += 0.01
    return bestdur

def calibrate_brlens_strat(tree,gap=0.2):
    cdef double f, l, oldest, youngest
    cdef node.Node n

    for n in tree.iternodes(order=1):
        if n.istip:
            f = n.strat[0]
            l = n.strat[1]
            n.upper = l - (gap / 2.0)
            if n.upper < 0.0:
                n.upper = 0.0
            n.lower = f + (gap / 2.0)
        elif n.istip == False:
            #if n == tree:
            #    continue
            oldest = 0.0
            for ch in n.children:
                if ch.lower > oldest:
                    oldest = ch.lower 
            for ch in n.children:
                ch.lower = oldest
            n.upper = oldest
            n.lower = n.upper+gap
    for n in tree.iternodes(order = 1):
        if n.istip:
            youngest = 10000000.0
            oldest = 0.0
            for ch in n.children:
                if ch.lower < youngest:
                    youngest = ch.lower 
                if ch.lower > oldest:
                    oldest = ch.lower 
           
            if oldest >= n.lower - 0.09:
                n.lower = oldest + 0.1
            if youngest <= n.upper + 0.09:
                n.upper = youngest - 0.1
        elif n.istip == False:
            oldest = 0.0
            for ch in n.children:
                if ch.lower > oldest:
                    oldest = ch.lower 
            for ch in n.children:
                ch.lower = oldest
            n.upper = oldest
            n.lower = n.upper+gap
    
    for n in tree.iternodes():
        n.length = n.lower - n.upper



def bds_dates(double p, double q, double r, node.Node tree):
    cdef double gap, mllen, f, l, hypanc_len, oldest , obs_range
    cdef node.Node n, ch
    hypanc_len = hyp_anc_mle(p,q,r)
    for n in tree.iternodes(order=1):
        if n.istip:
            f = n.strat[0]
            l = n.strat[1]
            obs_range = f-l
            mllen = single_range_mle(p,q,r,n)
            gap = (mllen-obs_range) / 2.0
            n.upper = l - gap
            if n.upper < 0.0:
                n.upper = 0.0
            n.lower = f + gap
            n.length = mllen
        else:
            #if n == tree:
            #    continue
            oldest = 0.0
            for ch in n.children:
                if ch.lower > oldest:
                    oldest = ch.lower ## fix length of younger sister taxon?
            for ch in n.children:
                ch.lower = oldest
            n.upper = oldest
            n.lower = n.upper+hypanc_len
            n.length = hypanc_len
    #print("\n\n\n")
    for n in tree.iternodes(order = 1):
        if n.istip:
            for ch in n.children:
                if ch.lower >= n.lower:
                    n.lower = ch.lower + 0.1
                if ch.lower < n.upper:
                    n.upper = ch.lower
                templen = n.lower - n.upper
                if templen != n.length:
                    n.length = templen

def bds_extinct_branch_ll(double p,double q,double r,node.Node n):
    cdef double f,l,obs_range,tf,tl,duration,p_like,q_like,r_like,brlik
    
    tf = n.lower
    tl = n.upper
    duration = tf-tl

    f = n.strat[0]
    l = n.strat[1]
    obs_range = f-l
    r_like = bd.calc_prob_range(r,duration,obs_range)
    p_like = bd.prob_n_obs_desc(p,q,r,len(n.children),duration)
    q_like = bd.prob_extinction_t(q,duration)
    brlik = np.log(r_like) + np.log(p_like) + np.log(q_like)
    return brlik

def bds_hyp_anc_ll(double p,double q,double r,node.Node n):
    cdef double tf,tl,duration,brlik

    tf = n.lower
    tl = n.upper
    duration = tf-tl
    brlik = bd.bds_hyp_anc_log_prob(p,q,r,duration)
    return brlik

def bds_loglike(double p, double q, double r, node.Node tree):
    cdef double small,tf,tl,f,l,obs_range,duration,r_like,p_like,q_like,brlik,treelik = 0.0
    cdef node.Node n

    small = 0.0001
    if r < small or p < small or q < small:
        return -100000000000.0

    for n in tree.iternodes():
        if n == tree and n.istip == False:
            continue
        if n.istip:
            if n.upper > 0.0:
                brlik = bds_extinct_branch_ll(p,q,r,n)
            elif n.upper == 0.0:
                #print("WARNING: extant species not yet implemented for bds")
                brlik = bds_extinct_branch_ll(p,q,r,n)
        else:
            brlik = bds_hyp_anc_ll(p,q,r,n)
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



