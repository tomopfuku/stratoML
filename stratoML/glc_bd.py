import numpy as np
import math 

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