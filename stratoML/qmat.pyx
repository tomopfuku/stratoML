import sys
import numpy as np
cimport numpy as np
np.import_array()
cimport cython
import smaps
import mfc
#@cython.wraparound(False)
#@cython.boundscheck(False)

cdef class Qmat:
    cdef public double[:,:] twostate,threestate,fourstate,fivestate,sixstate,sevenstate

    def __init__(self, double start_m, double start_l):
        #self.data = {}
        self.twostate = np.zeros((2**2,2**2),dtype=np.double)
        self.threestate = np.zeros((2**3,2**3),dtype=np.double)
        self.fourstate = np.zeros((2**4,2**4),dtype=np.double)
        self.fivestate = np.zeros((2**5,2**5),dtype=np.double)
        self.sixstate = np.zeros((2**6,2**6),dtype=np.double)
        self.sevenstate = np.zeros((2**7,2**7),dtype=np.double)
        self.init_q_matrices(start_m,start_l)

    def init_q_matrices(self, double m, double l):
        cdef int i
        #cdef double[:,:] qmat
        #cdef long[:,:] smap
        for i in range(2,5): ## TODO: need to increase stop to 8 after implementing full 7 state smaps
            #smap = smaps.get_smap(i)
            #qmat = self.get_qmat(i)
            self.update_mut_loss(i,m,l)

    def calc_p_mats(self, double t, int max_states = 6):
        cdef double[:,:,:] mats = np.zeros((max_states-1,2**max_states,2**max_states),dtype=np.double)
        cdef double[:,:] curp 
        cdef int i,j,k
        for i in range(2,max_states+1):
            curp = mfc.calc_p_matrix(self.get_qmat(i),t)
            for j in range(len(curp)):
                for k in range(len(curp)):
                    mats[i-2][j][k] = curp[j][k]
        return mats

    def update_all_qmats(self, double m, double l):
        cdef int i
        for i in range(2,5):
            self.update_mut_loss(i ,m ,l)

    def update_mut_loss(self, long nstates, double m, double l):
        cdef int i, j, ii, ndiff, nrow, nsti, nstj 
        cdef double[:,:] qmat = self.get_qmat(nstates)
        cdef long[:,:] smap = smaps.get_smap(nstates)
        cdef double scalar = 2.0 / float(nstates)
        
        nrow = len(qmat)
        for i in range(nrow):
            qmat[i][i] = 0.0
            if i == 0:
                continue
            for j in range(nrow):
                if i != j:
                    ndiff = 0
                    for ii in range(nstates):
                        if smap[i][ii] != smap[j][ii]:
                            ndiff += 1
                    nsti = np.add.reduce(smap[i])
                    nstj = np.add.reduce(smap[j])
                    if nstj - nsti < 0:
                        if ndiff == 1:
                            qmat[i][j] = l * scalar
                    elif nstj - nsti > 0:
                        if ndiff == 1:
                            qmat[i][j] = m * scalar

        for i in range(nrow):
            qmat[i][i] = -np.add.reduce(qmat[i])

    def get_qmat(self, long nstate):
        if nstate == 2:
            return self.twostate
        elif nstate == 3:
            return self.threestate
        elif nstate == 4:
            return self.fourstate
        elif nstate == 5:
            return self.fivestate
        elif nstate == 6:
            return self.sixstate
        elif nstate == 7:
            return self.sevenstate
        else:
            print("currently supports only up to 7-state characters")
            sys.exit()