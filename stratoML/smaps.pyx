import numpy as np


cdef long[:,:] smap2,smap3,smap4

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

def get_smap(int nstates):
    if nstates == 2:
        return smap2
    if nstates == 3:
        return smap3
    if nstates == 4:
        return smap4
