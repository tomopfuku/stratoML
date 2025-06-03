cdef class Node:
    cdef public bint istip
    cdef public double height,length,upper,lower,midpoint
    cdef public unsigned int num_occurrences
    cdef public str label
    cdef public object parent
    cdef public list children
    cdef public double[:,:] disc_traits
    cdef public double[:,:,:] budd_marginals
    cdef public double[:,:,:] timeslice_lv 
    cdef public double[:,:] scaling_factors # first dim is children, second is traits
    cdef public double[:] strat
    cdef public double[:,:,:] pmats
    cdef public int index, subtree, index_from_parent, parent_lv_index, midpoint_lv_index
