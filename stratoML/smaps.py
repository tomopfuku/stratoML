import numpy as np

smap2 = np.array([
[0, 0],
[1, 0],
[0, 1],
[1, 1],
],dtype=int)

smap3 = np.array([
[0, 0, 0],
[1, 0, 0],
[0, 1, 0],
[1, 1, 0],
[0, 0, 1],
[1, 0, 1],
[0, 1, 1],
[1, 1, 1],
],dtype=int)

smap4 = np.array([
[0, 0, 0, 0],
[1, 0, 0, 0],
[0, 1, 0, 0],
[1, 1, 0, 0],
[0, 0, 1, 0],
[1, 0, 1, 0],
[0, 1, 1, 0],
[1, 1, 1, 0],
[0, 0, 0, 1],
[1, 0, 0, 1],
[0, 1, 0, 1],
[1, 1, 0, 1],
[0, 0, 1, 1],
[1, 0, 1, 1],
[0, 1, 1, 1],
[1, 1, 1, 1],
],dtype=int)

smap5 = np.array([
[0, 0, 0, 0, 0],
[1, 0, 0, 0, 0],
[0, 1, 0, 0, 0],
[1, 1, 0, 0, 0],
[0, 0, 1, 0, 0],
[1, 0, 1, 0, 0],
[0, 1, 1, 0, 0],
[1, 1, 1, 0, 0],
[0, 0, 0, 1, 0],
[1, 0, 0, 1, 0],
[0, 1, 0, 1, 0],
[1, 1, 0, 1, 0],
[0, 0, 1, 1, 0],
[1, 0, 1, 1, 0],
[0, 1, 1, 1, 0],
[1, 1, 1, 1, 0],
[0, 0, 0, 0, 1],
[1, 0, 0, 0, 1],
[0, 1, 0, 0, 1],
[1, 1, 0, 0, 1],
[0, 0, 1, 0, 1],
[1, 0, 1, 0, 1],
[0, 1, 1, 0, 1],
[1, 1, 1, 0, 1],
[0, 0, 0, 1, 1],
[1, 0, 0, 1, 1],
[0, 1, 0, 1, 1],
[1, 1, 0, 1, 1],
[0, 0, 1, 1, 1],
[1, 0, 1, 1, 1],
[0, 1, 1, 1, 1],
[1, 1, 1, 1, 1],
],dtype=int)

smap6 = np.array([
[0, 0, 0, 0, 0, 0],
[1, 0, 0, 0, 0, 0],
[0, 1, 0, 0, 0, 0],
[1, 1, 0, 0, 0, 0],
[0, 0, 1, 0, 0, 0],
[1, 0, 1, 0, 0, 0],
[0, 1, 1, 0, 0, 0],
[1, 1, 1, 0, 0, 0],
[0, 0, 0, 1, 0, 0],
[1, 0, 0, 1, 0, 0],
[0, 1, 0, 1, 0, 0],
[1, 1, 0, 1, 0, 0],
[0, 0, 1, 1, 0, 0],
[1, 0, 1, 1, 0, 0],
[0, 1, 1, 1, 0, 0],
[1, 1, 1, 1, 0, 0],
[0, 0, 0, 0, 1, 0],
[1, 0, 0, 0, 1, 0],
[0, 1, 0, 0, 1, 0],
[1, 1, 0, 0, 1, 0],
[0, 0, 1, 0, 1, 0],
[1, 0, 1, 0, 1, 0],
[0, 1, 1, 0, 1, 0],
[1, 1, 1, 0, 1, 0],
[0, 0, 0, 1, 1, 0],
[1, 0, 0, 1, 1, 0],
[0, 1, 0, 1, 1, 0],
[1, 1, 0, 1, 1, 0],
[0, 0, 1, 1, 1, 0],
[1, 0, 1, 1, 1, 0],
[0, 1, 1, 1, 1, 0],
[1, 1, 1, 1, 1, 0],
[0, 0, 0, 0, 0, 1],
[1, 0, 0, 0, 0, 1],
[0, 1, 0, 0, 0, 1],
[1, 1, 0, 0, 0, 1],
[0, 0, 1, 0, 0, 1],
[1, 0, 1, 0, 0, 1],
[0, 1, 1, 0, 0, 1],
[1, 1, 1, 0, 0, 1],
[0, 0, 0, 1, 0, 1],
[1, 0, 0, 1, 0, 1],
[0, 1, 0, 1, 0, 1],
[1, 1, 0, 1, 0, 1],
[0, 0, 1, 1, 0, 1],
[1, 0, 1, 1, 0, 1],
[0, 1, 1, 1, 0, 1],
[1, 1, 1, 1, 0, 1],
[0, 0, 0, 0, 1, 1],
[1, 0, 0, 0, 1, 1],
[0, 1, 0, 0, 1, 1],
[1, 1, 0, 0, 1, 1],
[0, 0, 1, 0, 1, 1],
[1, 0, 1, 0, 1, 1],
[0, 1, 1, 0, 1, 1],
[1, 1, 1, 0, 1, 1],
[0, 0, 0, 1, 1, 1],
[1, 0, 0, 1, 1, 1],
[0, 1, 0, 1, 1, 1],
[1, 1, 0, 1, 1, 1],
[0, 0, 1, 1, 1, 1],
[1, 0, 1, 1, 1, 1],
[0, 1, 1, 1, 1, 1],
[1, 1, 1, 1, 1, 1],
],dtype=int)

smap7 = np.array([
[0, 0, 0, 0, 0, 0, 0],
[1, 0, 0, 0, 0, 0, 0],
[0, 1, 0, 0, 0, 0, 0],
[1, 1, 0, 0, 0, 0, 0],
[0, 0, 1, 0, 0, 0, 0],
[1, 0, 1, 0, 0, 0, 0],
[0, 1, 1, 0, 0, 0, 0],
[1, 1, 1, 0, 0, 0, 0],
[0, 0, 0, 1, 0, 0, 0],
[1, 0, 0, 1, 0, 0, 0],
[0, 1, 0, 1, 0, 0, 0],
[1, 1, 0, 1, 0, 0, 0],
[0, 0, 1, 1, 0, 0, 0],
[1, 0, 1, 1, 0, 0, 0],
[0, 1, 1, 1, 0, 0, 0],
[1, 1, 1, 1, 0, 0, 0],
[0, 0, 0, 0, 1, 0, 0],
[1, 0, 0, 0, 1, 0, 0],
[0, 1, 0, 0, 1, 0, 0],
[1, 1, 0, 0, 1, 0, 0],
[0, 0, 1, 0, 1, 0, 0],
[1, 0, 1, 0, 1, 0, 0],
[0, 1, 1, 0, 1, 0, 0],
[1, 1, 1, 0, 1, 0, 0],
[0, 0, 0, 1, 1, 0, 0],
[1, 0, 0, 1, 1, 0, 0],
[0, 1, 0, 1, 1, 0, 0],
[1, 1, 0, 1, 1, 0, 0],
[0, 0, 1, 1, 1, 0, 0],
[1, 0, 1, 1, 1, 0, 0],
[0, 1, 1, 1, 1, 0, 0],
[1, 1, 1, 1, 1, 0, 0],
[0, 0, 0, 0, 0, 1, 0],
[1, 0, 0, 0, 0, 1, 0],
[0, 1, 0, 0, 0, 1, 0],
[1, 1, 0, 0, 0, 1, 0],
[0, 0, 1, 0, 0, 1, 0],
[1, 0, 1, 0, 0, 1, 0],
[0, 1, 1, 0, 0, 1, 0],
[1, 1, 1, 0, 0, 1, 0],
[0, 0, 0, 1, 0, 1, 0],
[1, 0, 0, 1, 0, 1, 0],
[0, 1, 0, 1, 0, 1, 0],
[1, 1, 0, 1, 0, 1, 0],
[0, 0, 1, 1, 0, 1, 0],
[1, 0, 1, 1, 0, 1, 0],
[0, 1, 1, 1, 0, 1, 0],
[1, 1, 1, 1, 0, 1, 0],
[0, 0, 0, 0, 1, 1, 0],
[1, 0, 0, 0, 1, 1, 0],
[0, 1, 0, 0, 1, 1, 0],
[1, 1, 0, 0, 1, 1, 0],
[0, 0, 1, 0, 1, 1, 0],
[1, 0, 1, 0, 1, 1, 0],
[0, 1, 1, 0, 1, 1, 0],
[1, 1, 1, 0, 1, 1, 0],
[0, 0, 0, 1, 1, 1, 0],
[1, 0, 0, 1, 1, 1, 0],
[0, 1, 0, 1, 1, 1, 0],
[1, 1, 0, 1, 1, 1, 0],
[0, 0, 1, 1, 1, 1, 0],
[1, 0, 1, 1, 1, 1, 0],
[0, 1, 1, 1, 1, 1, 0],
[1, 1, 1, 1, 1, 1, 0],
[0, 0, 0, 0, 0, 0, 1],
[1, 0, 0, 0, 0, 0, 1],
[0, 1, 0, 0, 0, 0, 1],
[1, 1, 0, 0, 0, 0, 1],
[0, 0, 1, 0, 0, 0, 1],
[1, 0, 1, 0, 0, 0, 1],
[0, 1, 1, 0, 0, 0, 1],
[1, 1, 1, 0, 0, 0, 1],
[0, 0, 0, 1, 0, 0, 1],
[1, 0, 0, 1, 0, 0, 1],
[0, 1, 0, 1, 0, 0, 1],
[1, 1, 0, 1, 0, 0, 1],
[0, 0, 1, 1, 0, 0, 1],
[1, 0, 1, 1, 0, 0, 1],
[0, 1, 1, 1, 0, 0, 1],
[1, 1, 1, 1, 0, 0, 1],
[0, 0, 0, 0, 1, 0, 1],
[1, 0, 0, 0, 1, 0, 1],
[0, 1, 0, 0, 1, 0, 1],
[1, 1, 0, 0, 1, 0, 1],
[0, 0, 1, 0, 1, 0, 1],
[1, 0, 1, 0, 1, 0, 1],
[0, 1, 1, 0, 1, 0, 1],
[1, 1, 1, 0, 1, 0, 1],
[0, 0, 0, 1, 1, 0, 1],
[1, 0, 0, 1, 1, 0, 1],
[0, 1, 0, 1, 1, 0, 1],
[1, 1, 0, 1, 1, 0, 1],
[0, 0, 1, 1, 1, 0, 1],
[1, 0, 1, 1, 1, 0, 1],
[0, 1, 1, 1, 1, 0, 1],
[1, 1, 1, 1, 1, 0, 1],
[0, 0, 0, 0, 0, 1, 1],
[1, 0, 0, 0, 0, 1, 1],
[0, 1, 0, 0, 0, 1, 1],
[1, 1, 0, 0, 0, 1, 1],
[0, 0, 1, 0, 0, 1, 1],
[1, 0, 1, 0, 0, 1, 1],
[0, 1, 1, 0, 0, 1, 1],
[1, 1, 1, 0, 0, 1, 1],
[0, 0, 0, 1, 0, 1, 1],
[1, 0, 0, 1, 0, 1, 1],
[0, 1, 0, 1, 0, 1, 1],
[1, 1, 0, 1, 0, 1, 1],
[0, 0, 1, 1, 0, 1, 1],
[1, 0, 1, 1, 0, 1, 1],
[0, 1, 1, 1, 0, 1, 1],
[1, 1, 1, 1, 0, 1, 1],
[0, 0, 0, 0, 1, 1, 1],
[1, 0, 0, 0, 1, 1, 1],
[0, 1, 0, 0, 1, 1, 1],
[1, 1, 0, 0, 1, 1, 1],
[0, 0, 1, 0, 1, 1, 1],
[1, 0, 1, 0, 1, 1, 1],
[0, 1, 1, 0, 1, 1, 1],
[1, 1, 1, 0, 1, 1, 1],
[0, 0, 0, 1, 1, 1, 1],
[1, 0, 0, 1, 1, 1, 1],
[0, 1, 0, 1, 1, 1, 1],
[1, 1, 0, 1, 1, 1, 1],
[0, 0, 1, 1, 1, 1, 1],
[1, 0, 1, 1, 1, 1, 1],
[0, 1, 1, 1, 1, 1, 1],
[1, 1, 1, 1, 1, 1, 1],
],dtype=int)

def get_smap(nstates):
    if nstates == 2:
        return smap2
    if nstates == 3:
        return smap3
    if nstates == 4:
        return smap4
    if nstates == 5:
        return smap5
    if nstates == 6:
        return smap6
    if nstates == 7:
        return smap7
