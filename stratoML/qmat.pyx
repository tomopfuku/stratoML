import sys
import numpy as np
cimport numpy as np
np.import_array()
cimport cython
import smaps
import mfc, bd
import lam_mat
cimport lam_mat
from scipy.linalg import expm
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
        for i in range(2,8): ## TODO: need to increase stop to 8 after implementing full 7 state smaps
            #smap = smaps.get_smap(i)
            #qmat = self.get_qmat(i)
            self.update_mut_loss(i,m,l)

    def calc_single_p_mat(self, double t, int ss):
        cdef double[:,:] p
        p = mfc.calc_p_matrix(self.get_qmat(ss),t)
        return p

    def calc_single_m_mat(self, double dt, lam_mat.lam_mat lam_mats, double[:] bds_rates, int ss):
        cdef double[:,:] curq = self.get_qmat(ss)
        cdef double[:,:] curlam = lam_mats.get_ratemat(ss)
        cdef np.ndarray[np.float64_t, ndim=2] curm = np.asarray(curq).copy()
        cdef np.ndarray[np.float64_t, ndim=2] curM
        cdef double E_t = bd.calc_extinction_prob_eq(bds_rates[0], bds_rates[1], bds_rates[2])
        cdef double diag
        cdef int j, k, nstates = 1 << ss

        for j in range(nstates):
            for k in range(nstates):
                curm[j, k] += curlam[j, k] * E_t

            if j == 0:
                diag = 0.0
            else:
                diag = curq[j, j] - (bds_rates[0] + bds_rates[1] + bds_rates[2]) + (bds_rates[0] * E_t)
            curm[j, j] = diag

        curM = expm(curm * dt)
        return curM


    def calc_m_mats(self, double dt, lam_mat.lam_mat lam_mats, double[:] bds_rates, int max_states = 7):
        cdef double[:,:,:] mats = np.zeros((int(max_states)-1,int(2**max_states),int(2**max_states)),dtype=np.double)
        cdef double[:,:] curM
        cdef int i,j,k

        for i in range(2,max_states+1):
            curM = self.calc_single_m_mat(dt, lam_mats, bds_rates, i)
            for j in range(len(curM)):
                for k in range(len(curM[j])):
                    mats[i-2][j][k] = curM[j, k]
        return mats


    def calc_m_mats_loop(self, double dt, lam_mat.lam_mat lam_mats, double[:] bds_rates, int max_states = 7):
        cdef double[:,:,:] mats = np.zeros((int(max_states)-1,int(2**max_states),int(2**max_states)),dtype=np.double)
        cdef double[:,:] curm, curq, curlam
        cdef double diag_penalty, diag_reward
        cdef int i,j,k, N = int(2**max_states)
        cdef double[:] lam_tot = np.zeros(N, dtype=np.float64)
        E_t = bd.calc_extinction_prob_eq(bds_rates[0], bds_rates[1], bds_rates[2])

        for i in range(2,max_states+1):
            curm = np.zeros((N, N), dtype=np.float64)
            curq = self.get_qmat(i)
            curlam = lam_mats.get_ratemat(i)
            lam_tot = np.zeros(N, dtype=np.float64)
            for j in range(N):
                for k in range(N):
                    lam_tot[j] += curlam[i,j]
            #curm = self.get_qmat(i) + ( lam_mats.get_ratemat(i).T * E_t )
            for j in range(N):
                diag_penalty = -(lam_tot[i] + bds_rates[1] + bds_rates[2])
                diag_reward = lam_tot[i] * E_t
                for k in range(N):
                    if k > 2**i or j > 2**i:
                        continue
                    curm[j, k] = curq[j, k] + ( curlam[j, k] * E_t )
                    if j == k:
                        curm[j, k] += diag_penalty + diag_reward
        return mats


    def calc_p_mats(self, double t, int max_states = 7):
        cdef double[:,:,:] mats = np.zeros((int(max_states)-1,int(2**max_states),int(2**max_states)),dtype=np.double)
        cdef double[:,:] curp
        cdef int i,j,k
        for i in range(2,max_states+1):
            curp = self.calc_single_p_mat(t, i)
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
        cdef double scalar = 1.0 #2.0 / float(nstates)

        nrow = len(qmat)
        for i in range(nrow):
            qmat[i][i] = 0.0
            if i == 0:
                continue
            for j in range(nrow):
                if i != j:
                    if j == 0:  ## remove null state
                        continue
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
