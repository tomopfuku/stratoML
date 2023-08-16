import sys
import math

def prop_pres_taxa(b,d,r): 
    top = (b+d+r) - math.sqrt((b+d+r)**2 - (4.0*b*d))
    bot = 2.0 * b
    return 1.0 - (top / bot)

def calc_prob_unsamp(b,d,r,unsamp):
    Ps = calc_Ps(b,d,r)
    top = ((r + (Ps * b) )**3.0) * (unsamp**2.0) * (math.exp(-(r+(Ps*b))*unsamp))
    bot = 2.0
    return top / bot

def calc_prob_range(r,true_range,obs_range):
    if obs_range > 0.0:
        like = (r**2.0) * (true_range-obs_range) * math.exp(-r*(true_range-obs_range))
    elif obs_range == 0.0:
        like = r*(true_range)*math.exp(-r*true_range)
    return like

def calc_non_pres(r,true_range):
    return math.exp(-r*true_range)    

def bds_hyp_anc_log_prob(p,q,r,time):
    prob_pt = math.log(1.0-math.exp(-p*time))
    #prob_pt = (p*time)*math.exp(-p*time)
    prob_qt = -q*time
    prob_rt = -r*time
    return (prob_pt + prob_qt + prob_rt )# / prob_qt

def bds_hyp_anc_prob(p,q,r,time):
    prob_pt = (1.0-math.exp(-p*time))
    #prob_pt = (p*time)*math.exp(-p*time)
    prob_qt = math.exp(-q*time)
    prob_rt = math.exp(-r*time)
    return (prob_pt * prob_qt * prob_rt )# / prob_qt

def prob_extinction_t(q,time):
    prob_qt = math.exp(-q*time) - math.exp(-q*(time+0.000001))
    return prob_qt

def prob_n_desc(p,n,duration):
    return (math.exp(-p*duration)*((p*duration)**n)) / math.factorial(n)
    
def prob_n_obs_desc(p,q,r,n_obs,duration):
    p_pres = prop_pres_taxa(p,q,r)
    p_unobs = 1.0 - p_pres
    p_obs = p_pres ** n_obs
    p_n_obs = prob_n_desc(p,n_obs,duration) #(math.exp(-p*duration)*((p*duration)**n_obs)) / math.factorial(n_obs)
    marg_p = p_n_obs * p_obs
    for n in range(n_obs+1,n_obs+10):
        cur_p = prob_n_desc(p,n,duration) 
        #n_unobs = n - n_obs
        marg_p += (cur_p * p_obs * p_unobs)
        p_unobs *= p_unobs
    return marg_p

#def BDS_range_prob(p,q,r,time):

    
def expect_gap(r):
    expect = 1.0 / r
    return expect