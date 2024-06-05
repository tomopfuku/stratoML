cdef class Node:
    cdef public bint istip
    cdef public double height,length,upper,lower
    cdef public unsigned int num_occurrences
    cdef public str label
    cdef public object parent
    cdef public list children
    cdef public double[:,:] disc_traits
    cdef public double[:] strat
    cdef public double[:,:,:] pmats
    cdef public int index, subtree
