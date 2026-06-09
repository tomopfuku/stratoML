#import sys
import math
import numpy as np

# prob of observing a _clade_ of unknown size from Didier 2017 via Wagner 2019
cpdef double prop_pres_taxa(double b,double d,double r): 
    cdef double top,bot
    top = (b+d+r) - math.sqrt((b+d+r)**2 - (4.0*b*d))
    bot = 2.0 * b
    return 1.0 - (top / bot)

# Pp given by Solow and Smith 1996 via Foote 1997
cpdef double prob_species_pres(double q, double r):
    return (r/(q+r))

"""
def calc_prob_unsamp(b,d,r,unsamp):
    Ps = calc_Ps(b,d,r)
    top = ((r + (Ps * b) )**3.0) * (unsamp**2.0) * (math.exp(-(r+(Ps*b))*unsamp))
    bot = 2.0
    return top / bot
"""

# probability of observing strat range given lineage duration from Foote 1997
cpdef double calc_prob_range(double r, double true_range, double obs_range):
    cdef double like

    if obs_range > 0.0:
        like = (r**2.0) * (true_range-obs_range) * math.exp(-r*(true_range-obs_range))
    elif obs_range == 0.0:
        like = r*(true_range)*math.exp(-r*true_range)
    return like

# poisson probability of not preserving over time_range
cpdef double calc_non_pres(double r, double true_range):
    return math.exp(-r*true_range)    

cpdef double bds_hyp_anc_log_prob(double p, double q, double r, double time):
    cdef double prob_pt,prob_qt,prob_rt

    prob_pt = math.log(1.0-math.exp(-p*time))
    #prob_pt = (p*time)*math.exp(-p*time)
    prob_qt = -q*time
    prob_rt = -r*time
    return (prob_pt + prob_qt + prob_rt )# / prob_qt

cpdef double bds_hyp_anc_prob(double p, double q, double r, double time):
    cdef double prob_pt,prob_qt,prob_rt

    prob_pt = (1.0-math.exp(-p*time))
    #prob_pt = (p*time)*math.exp(-p*time)
    prob_qt = math.exp(-q*time)
    prob_rt = math.exp(-r*time)
    return (prob_pt * prob_qt * prob_rt )# / prob_qt

cpdef double prob_extinction_t(double q, double time):
    cdef double prob_qt
    #prob_qt = math.exp(-q*time) - math.exp(-q*(time+0.000001))
    prob_qt = q * math.exp(-q*time)
    return prob_qt

cpdef double prob_n_desc(double p, double n,double duration):
    return (math.exp(-p*duration)*((p*duration)**n)) / math.factorial(int(n))
    
cpdef double prob_n_obs_desc(double p,double q,double r,int n_obs,double duration):
    cdef double p_pres, p_unobs, p_obs, p_n_obs, marg_p, cur_p
    cdef int n
    #cdef double float_n_obs = float(n_obs)
 
    p_pres = prop_pres_taxa(p,q,r)
    p_unobs = 1.0 - p_pres
    p_obs = p_pres ** float(n_obs)
    p_n_obs = prob_n_desc(p,float(n_obs),duration) #(math.exp(-p*duration)*((p*duration)**n_obs)) / math.factorial(n_obs)
    marg_p = p_n_obs * p_obs
    for n in range(n_obs+1,n_obs+20):
        cur_p = prob_n_desc(p,float(n),duration) 
        #n_unobs = n - n_obs
        marg_p += (cur_p * p_obs * p_unobs)
        p_unobs *= p_unobs
    return marg_p

def expect_gap(r):
    expect = 1.0 / r
    return expect



def calc_mean_extinction_prob(lam, mu, psi, t_start, t_end):
    D = np.sqrt((lam + mu + psi)**2 - 4*lam*mu)
    x1 = (lam + mu + psi - D) / (2 * lam)
    x2 = (lam + mu + psi + D) / (2 * lam)
    
    def integral_E(t):
        term1 = x1 * t
        log_term = np.log(x2 - x1 * np.exp(-D * t))
        term2 = ((x1 - x2) / D) * log_term
        return term1 + term2

    total_area = integral_E(t_end) - integral_E(t_start)
    duration = t_end - t_start
    
    return total_area / duration

def calc_extinction_prob(lam, mu, psi, t):
    D = np.sqrt((lam + mu + psi)**2 - 4*lam*mu)
    x1 = (lam + mu + psi - D) / (2 * lam)
    x2 = (lam + mu + psi + D) / (2 * lam)
    
    num = x1 * x2 * (1 - np.exp(-D * t))
    denom = x2 - (x1 * np.exp(-D * t))
    E_t = num / denom
    return E_t

def calc_extinction_prob_eq(lam, mu, psi):
    D = np.sqrt((lam + mu + psi)**2 - 4*lam*mu)
    
    x1 = (lam + mu + psi - D) / (2 * lam)
    E_t = x1
    return E_t