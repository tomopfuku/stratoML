import numpy as np
cimport numpy as np
np.import_array()
cimport cython
import sys
from scipy.linalg import expm

#@cython.wraparound(False)
#@cython.boundscheck(False)

cdef long[:,:] twostate,threestate,fourstate
twostate = np.array([[0,0],[0,1],[1,0],[1,1]],dtype=int)
threestate = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,0,1],[1,1,1]],dtype=int)
fourstate = np.array([
    [0,0,0,0],
    [1,0,0,0],
    [0,1,0,0],
    [0,0,1,0],
    [0,0,0,1],
    [1,1,0,0],
    [1,0,1,0],
    [1,0,0,1],
    [0,1,1,0],
    [0,0,1,1],
    [0,1,0,1],
    [1,1,1,0],
    [0,1,1,1],
    [1,0,1,1],
    [1,1,0,1],
    [1,1,1,1]
],dtype=int)

def mat_mult(double[:,:] mat, double m):
    cdef int i,j,nrow
    nrow = len(mat)
    cdef double[:,:] newmat = np.empty((nrow,nrow),dtype = np.double)
    for i in range(len(mat)):
        for j in range(len(mat)):
            newmat[i][j] = mat[i][j] * m
    return newmat

def calc_p_matrix(double[:,:] qmat, double t):
    cdef double[:,:] pmat
    #for i in mat_mult(qmat,t):
    #    print(list(i))
    pmat = expm(mat_mult(qmat,t))
    return pmat

def create_q_matrix(int nstates,double m, double l):
    cdef double[:,:] qmat
    cdef int matsize
    matsize = 2**nstates
    qmat = np.zeros((matsize,matsize),dtype=np.double)
    smap = statemap(nstates)
    update_mut_loss(qmat,smap,m,l)
    #for row in qmat:
    #    print(list(row))                    
    
    return qmat, smap

def statemap(int nstates):
    cdef long[:,:] smap
    if nstates == 2:
        smap = twostate
    elif nstates == 3:
        smap = threestate
    elif nstates == 4:
        smap = fourstate
    return smap



def update_mut_loss(double[:,:] qmat, long[:,:] smap, double m, double l):
    cdef int i, j, ii, ndiff, nrow, nsti, nstj, nstates 
    nstates = len(smap[0])
    nrow = len(qmat)
    for i in range(nrow):
        if i == 0:
            continue
        for j in range(nrow):
            if i != j:
                ndiff = 0
                for ii in range(nstates):
                    if smap[i][ii] == smap[j][ii]:
                        continue
                    else:
                        ndiff += 1
                nsti = np.add.reduce(smap[i])
                nstj = np.add.reduce(smap[j])
                if nstj - nsti < 0:
                    if ndiff == 1:
                        qmat[i][j] = l
                elif nstj - nsti > 0:
                    if ndiff == 1:
                        qmat[i][j] = m

    for i in range(nrow):
        qmat[i][i] = -np.add.reduce(qmat[i])


