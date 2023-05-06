#from libcpp cimport bool
#from copy import deepcopy
import numpy as np

cdef:
    unsigned int PREORDER = 0
    unsigned int POSTORDER = 1

cdef class Node:
    #cdef public dict data
    cdef public bint istip
    cdef public double height,length,upper,lower
    cdef public unsigned int num_occurrences
    cdef public str label
    cdef public object parent
    cdef public list children
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
        #self.cont_traits = np.array([],dtype=np.double)

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
        self.nchildren += 1

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
            object child, d

        if order == PREORDER:
            yield self
        #print [i.label for i in self.children]
        for child in self.children:
            for d in child.iternodes(order):
                yield d
        if order == POSTORDER:
            yield self


    def descendants(self, order=PREORDER, v=None):
        if v is None:
            v = []
        #assert order in ("PREORDER", "POSTORDER")
        for child in self.children:
            if order == PREORDER:
                v.append(child)
            else:
                v.insert(0, child)
            if child.children:
                child.descendants(order, v)
        return v

    def prune(self):
        p = self.parent
        if p:
            p.remove_child(self)
        return p


