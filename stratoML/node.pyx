#from libcpp cimport bool
#from copy import deepcopy
import numpy as np
cimport numpy as np
np.import_array()
cimport cython
#from cpython cimport array
#import array


cdef:
    unsigned int PREORDER = 0
    unsigned int POSTORDER = 1

@cython.wraparound(False)
@cython.boundscheck(False)

cdef class Node:
    #cdef public dict data
    cdef public bint istip
    cdef public double height,length,upper,lower
    cdef public unsigned int num_occurrences
    cdef public str label
    cdef public object parent
    cdef public list children
    #cdef public long[:] disc_starts
    cdef public long[:] disc_traits
    cdef public double[:] strat
    #cdef public double[:] cont_traits

    def __init__(self):
        #self.data = {}
        self.istip = False
        self.label = ""
        self.length = 0.
        self.parent = None
        self.children = []
        self.height = 0.0
        self.upper = 0.0
        self.lower = 0.0
        self.num_occurrences = 0
        self.strat = np.array([0.0,0.0],dtype=np.double)
        #self.disc_starts = np.array([],dtype=int)
        self.disc_traits = np.array([],dtype=int)
        #self.cont_traits = np.array([],dtype=np.double)

    def add_disc_traits(self, list traitls):
        self.disc_traits = np.array(traitls,dtype=int)
        print(traitls)
        print(self.disc_traits)

    """def add_disc_traits(self, list traitls):
        N = sum(map(len, traitls))
        #starts = np.empty(N, dtype=int) 
        starts = np.empty(len(traitls)+1, dtype=int) 
        traits = np.empty(N, dtype=int)

        starts[0], cnt = 0, 0
        for i,states in enumerate(traitls):
            for el in states:
                traits[cnt] = el
                cnt += 1       # update index in the flattened array for the next element
            starts[i+1] = cnt  # remember the start of the next list
        
        print(starts)
        print(type(starts))
        print(traits)
        print(np.array([starts,traits],dtype=int))
        #self.disc_traits = np.array([starts,traits],dtype=int)
        self.disc_starts = starts
        self.disc_states = traits"""

    def get_newick_repr(self,showbl=False,show_rate=False):
        ret = ""
        for i in range(len(self.children)):
            if i == 0:
                ret += "("
            ret += self.children[i].get_newick_repr(showbl,show_rate)
            if i == len(self.children)-1:
                ret += ")"
            else:
                ret += ","
        if self.label != None:
            ret += self.label
        if showbl == True:
            ret += ":" + str(self.length)
        if show_rate ==True:
            ret += ":" + str(self.sigsq)
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

    def iternodes(self, unsigned int order=PREORDER):#, v=None):
        cdef:
            object i, d

        if order == PREORDER:
            yield self
        #print [i.label for i in self.children]
        for i in range(len(self.children)):
            for d in self.children[i].iternodes(order):
                yield d
        if order == POSTORDER:
            yield self

    def prune(self):
        p = self.parent
        if p:
            p.remove_child(self)
        return p


