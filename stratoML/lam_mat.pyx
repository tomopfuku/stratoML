import sys
import numpy as np
cimport numpy as np
np.import_array()
cimport cython
import smaps
import mfc
#@cython.wraparound(False)
#@cython.boundscheck(False)

def check_all_anc_in_desc(long[:] anc_map, long[:] desc_map):
    cdef int i#, j
    for i in range(len(desc_map)):
        if anc_map[i] == 0 and desc_map[i] == 1:
            return False
    return True        

def nchoosek(long n, long k):
    if k > n:
        return 0
    if k == 0 or k == n:
        return 1
        
    if k > n - k:
        k = n - k
        
    cdef long result = 1
    cdef long i
    
    for i in range(1, k + 1):
        result = result * (n - k + i) // i
        
    return result

cdef class lam_mat:
    cdef public double[:,:] twostate,threestate,fourstate,fivestate,sixstate,sevenstate

    def __init__(self, double start_lam_glob, double start_lam_sub, double start_lam_exp = 0, double start_lam_jump = 0):
        self.twostate = np.zeros((2**2,2**2),dtype=np.double)
        self.threestate = np.zeros((2**3,2**3),dtype=np.double)
        self.fourstate = np.zeros((2**4,2**4),dtype=np.double)
        self.fivestate = np.zeros((2**5,2**5),dtype=np.double)
        self.sixstate = np.zeros((2**6,2**6),dtype=np.double)
        self.sevenstate = np.zeros((2**7,2**7),dtype=np.double)
        """
        max_states = 7
        self.twostate = np.zeros((2**max_states,2**max_states),dtype=np.double)
        self.threestate = np.zeros((2**max_states,2**max_states),dtype=np.double)
        self.fourstate = np.zeros((2**max_states,2**max_states),dtype=np.double)
        self.fivestate = np.zeros((2**max_states,2**max_states),dtype=np.double)
        self.sixstate = np.zeros((2**max_states,2**max_states),dtype=np.double)
        self.sevenstate = np.zeros((2**max_states,2**max_states),dtype=np.double)
        """
        self.init_lam_matrices(start_lam_glob, start_lam_sub, start_lam_exp, start_lam_jump)

    def init_lam_matrices(self, double lglob, double lsub, double lexp, double ljump):
        cdef int i
        for i in range(2,8): ## TODO: need to increase stop to 8 after implementing full 7 state smaps
            self.update_rates(i ,lglob, lsub, lexp, ljump)

    """
    def calc_single_p_mat(self, double t, int ss):
        cdef double[:,:] p
        p = mfc.calc_p_matrix(self.get_qmat(ss),t)
        return p

    def calc_p_mats(self, double t, int max_states = 7):
        cdef double[:,:,:] mats = np.zeros((int(max_states)-1,int(2**max_states),int(2**max_states)),dtype=np.double)
        cdef double[:,:] curp
        cdef int i,j,k
        for i in range(2,max_states+1):
            curp = mfc.calc_p_matrix(self.get_qmat(i),t)
            for j in range(len(curp)):
                for k in range(len(curp)):
                    mats[i-2][j][k] = curp[j][k]
        return mats
    """

    def update_all_mats(self, double lglob, double lsub, double lexp=0.0, double ljump=0.0):
        cdef int i
        for i in range(2, 8):
            self.update_rates(i ,lglob, lsub, lexp, ljump)

    def update_rates(self, long nstates, double lglob, double lsub, double ljump = 0.0, double lexp = 0.0):
        cdef int i, j, ii, ndiff, nrow, nsti, nstj, nsub, n_tier, subs 
        cdef double[:,:] ratemat = self.get_ratemat(nstates)
        cdef long[:,:] smap = smaps.get_smap(nstates)
        cdef double scalar = 1.0 / float(nstates)

        nrow = len(ratemat)
        #nrow = 2**nstates
        for i in range(nrow):
            ratemat[i][i] = 0.0
            if i == 0:
                continue

            nsti = np.add.reduce(smap[i])
            n_tier = nsti - 1 
            #nsub = 0
            for j in range(nrow):
                if i != j:
                    if j == 0:  ## remove null state
                        continue
                    ndiff = 0
                    for ii in range(nstates):
                        if smap[i][ii] != smap[j][ii]:
                            ndiff += 1
                    nstj = np.add.reduce(smap[j])
                    if nstj - nsti < 0 and check_all_anc_in_desc(smap[i], smap[j]):
                        #if ndiff == 1:
                        #nsub += 1
                        subs = nchoosek(nsti, nstj)
                        weight = 1.0 / ( float(n_tier) * float(subs) )
                        #if nsti == 4:# and nstj == 2:
                        #    print(nsti, nstj, subs, weight)
                        #    exit()
                        ratemat[i][j] = lsub * weight 
                    elif nstj - nsti > 0:
                        #if ndiff == 1:
                        ratemat[i][j] = lexp * scalar
                        if lexp > 0:
                            print("ERROR: CLADOGENETIC EXPANSION NOT IMPLEMENTED YET")
                            exit()

                    elif nsti == nstj == 1 and ljump > 0.0:
                        ratemat[i][j] = ljump * float(nstates-1)
            """
            for j in range(nrow):
                if i != j:
                    if j == 0:  ## remove null state
                        continue

                    nstj = np.add.reduce(smap[j])
                    if nstj - nsti < 0:
                        ratemat[i][j] *= 1.0 / float(nsub)
            """
            if nsti > 1:
                ratemat[i][i] = lglob - lsub 
            else:
                ratemat[i][i] = lglob - ljump # lsym + lsub ## this needs to be modified for lexp and ljump


    def get_ratemat(self, long nstate):
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
