cimport cython
#from libcpp.list cimport list as cpplist
from libcpp.vector cimport vector

@cython.wraparound(False)
@cython.boundscheck(False)

cdef class morph:
    #cdef public dict data
    cdef public str label
    cdef public cpplist[str] seq

    def __init__(self):
        self.label = ""
        self.seq = seq