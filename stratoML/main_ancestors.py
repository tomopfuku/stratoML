import bd
import numpy as np
import matplotlib.pyplot as plt

def calc_all_desc_probs(p,q,r,obs_range,duration):
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

def calc_like(p,q,r,n_chld,obs_range,duration):
    r_like = np.log(bd.calc_prob_range(r,duration,obs_range))
    p_like = np.log(bd.prob_n_obs_desc(p,q,r,n_chld,duration))
    q_like = np.log(bd.prob_extinction_t(q,duration))
    ll = sum([r_like,p_like,q_like])
    return ll

if __name__ == "__main__":
    p = 0.07
    q = 0.07
    r = 4.8

    obs_range = 40
    duration = 40.2
    calc_all_desc_probs(p,q,r,obs_range,duration)

    print(calc_like(p,q,r,2,obs_range,duration))


    p = 0.6
    q = 0.4
    r = 1.3

    obs_range = 2.3
    durations=[]
    curdur = obs_range+0.01
    for i in range(1000):
        durations.append(curdur)
        curdur += 0.01



    likes = [calc_like(p,q,r,2,obs_range,i) for i in durations]
    plt.plot(durations[:500],likes[:500])
    plt.show()
