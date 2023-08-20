import bd
import numpy as np

def calc_all_desc_probs(p,q,r,n_chld,obs_range,duration):
    pvec = []
    for n_chld in range(50):
        r_like = np.log(bd.calc_prob_range(r,duration,obs_range))
        p_like = np.log(bd.prob_n_obs_desc(p,q,r,n_chld,duration))
        q_like = np.log(bd.prob_extinction_t(q,duration))
        ll = sum([r_like,p_like,q_like])
        pvec.append(ll)

        
    best = max(pvec)
    rel = []
    for i in pvec:
        rll = np.exp(-.5*(best-i))
        rel.append(rll)

    sumrel = np.sum(rel)

    for n_chld,val in enumerate(rel):
        print(n_chld,"descendants:",val/sumrel)
        if n_chld > 9:
            break



if __name__ == "__main__":
    p = 0.07
    q = 0.07
    r = 4.8

    obs_range = 40
    duration = 40.2
    calc_all_desc_probs(p,q,r,n_chld,obs_range,duration)