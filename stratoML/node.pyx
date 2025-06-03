import numpy as np
cimport numpy as np
np.import_array()
cimport cython
import qmat
cimport qmat
import sys
#import mfc


cdef:
    unsigned int PREORDER = 0
    unsigned int POSTORDER = 1

#@cython.wraparound(False)
#@cython.boundscheck(False)

cdef class Node:
    #cdef public dict data
    cdef public bint istip
    cdef public double height,length,upper,lower,midpoint
    cdef public unsigned int num_occurrences
    cdef public str label
    cdef public object parent
    cdef public list children
    #cdef public long[:] disc_starts
    cdef public double[:,:] disc_traits
    cdef public double[:,:,:] budd_marginals  # for simulation
    cdef public double[:,:,:] timeslice_lv
    cdef public double[:,:] scaling_factors # first dim is children, second is traits
    cdef public double[:] strat
    cdef public double[:,:,:] pmats
    #cdef public double[:] cont_traits
    cdef public int index, subtree, index_from_parent, parent_lv_index, midpoint_lv_index

    def __init__(self):
        #self.data = {}
        self.istip = False
        self.label = ""
        self.length = 0.
        self.parent = None
        self.children = []
        #self.height = 0.0
        self.upper = 0.0
        self.lower = 0.0
        self.midpoint = 0.0
        #self.num_occurrences = 0
        self.strat = np.array([0.0,0.0],dtype=np.double)
        self.disc_traits = np.array([[]],dtype=np.double)
        self.scaling_factors = np.array([[]],dtype=np.double)
        self.budd_marginals = np.array([[[]]],dtype=np.double) # in same order as self.children, these give marginal probs of being in each state at each budding point along branch
        self.timeslice_lv = np.array([[[]]],dtype=np.double) # ordered from tip toward root, with budding descendants and midpoint along the way 
        #self.cont_traits = np.array([],dtype=np.double)
        self.pmats = np.array([[[]]],dtype=np.double) # these are Pmats from either the start or midpoint of the branch
        #self.budd_pmats = np.array([[[]]],dtype=np.double) # these are Pmats from  
        self.index = 0
        self.index_from_parent = 0
        self.parent_lv_index = 0
        self.midpoint_lv_index = 0

    def update_pmat(self, qmat.Qmat ratemats, int maxstates, str mode="end"):
        cdef double t
        cdef unsigned int i
        #for i in range(len(self.children)):

        if mode != "end" and mode != "mid":
            print("transition probability matrices must be calculated assuming traits are observed at either the endpoint or midpoint.")
            print("valid options are \"end\" and \"mid\".")
            sys.exit()

        if mode == "end" or self.istip == False:
            t = self.length
        elif mode == "mid":
            t = self.length / 2.0
        self.pmats = ratemats.calc_p_mats(t, maxstates)


    def add_disc_traits(self, list traitls, long[:] ss):
        cdef double trait_freq
        cdef int i, cur_trait, j, nstate
        cdef double[:,:] trait_probs = np.zeros((len(traitls),128),dtype=np.double)
        for i in range(len(traitls)):
            cur_trait = traitls[i]
            if cur_trait != -9:
                trait_probs[i][cur_trait] = 1.0
            elif cur_trait == -9: # plug in flat priors for missing traits
                nstate = 2 ** int(ss[i])
                trait_freq = 1.0 #/ float(nstate)
                for j in range(len(trait_probs[i])):
                    if j == nstate:
                        break
                    trait_probs[i][j] = trait_freq
        self.disc_traits = trait_probs


    def get_newick_repr(self, bint showbl=False):
        cdef unsigned int i
        cdef str ret
        ret = ""
        for i in range(len(self.children)):
            if i == 0:
                ret += "("
            ret += self.children[i].get_newick_repr(showbl)
            if i == len(self.children)-1:
                ret += ")"
            else:
                ret += ","
        if self.label != None:
            ret += self.label
        if showbl == True:
            ret += ":" + str(self.length)
        return ret

    def add_child(self, object child):
        assert child not in self.children
        self.children.append(child)
        child.parent = self
        #self.nchildren += 1

    def remove_child(self, object child):
        assert child in self.children
        self.children.remove(child)
        child.parent = None
        self.nchildren -= 1


    def prune_from_node(self):
        for i in self.descendants("POSTORDER"):
            if len(self.children) == 0:
                self.prune()

    def leaves(self):
        return [ n for n in self.iternodes() if n.istip ]

    
    def get_sib(self):
        cdef int i
        if self.parent == None:
            return "cannot return sibling for root node"
        if len(self.parent.children) > 2:
            return "cannot return single sibling for multifurcating node"
        if len(self.parent.children) == 1:
            return None
        for i in range(2):
            if self.parent.children[i] != self:
                return self.parent.children[i]

    def inorder(self):
        if len(self.children)>0:
            for i,child in enumerate(self.children):
                if i == 0:
                    for d in child.inorder():
                        yield d
                    yield self
                else:
                    for d in child.inorder():
                        yield d
        else:
            yield self


    def iternodes(self, unsigned int order=PREORDER):#, v=None):
        cdef:
            unsigned int i
            Node d

        if order == PREORDER:
            yield self
        #print [i.label for i in self.children]
        for i in range(len(self.children)):
            for d in self.children[i].iternodes(order):
                yield d #self.children[i][d]
        if order == POSTORDER:
            yield self

    def prune(self):
        p = self.parent
        if p:
            p.remove_child(self)
        return p

        

