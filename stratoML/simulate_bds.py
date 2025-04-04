#import node
import random

def sim_n_desc(p,q,r,duration):
    anc = []
    curtime = 0.
    while curtime < duration:
        w = random.random()
        curtime += w
        if curtime <= duration:
            chld = []

            anc.append(chld)

if __name__ == "__main__":
    sim_n_desc(0.1,0.1,2.3,10.)
    
