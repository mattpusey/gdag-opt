"""Conditional independeces and semigraphoid closure

Based on Milian Studeny, "Complexity of Structural Models"
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.47.4379
"""
# Copyright (C) 2014 Matthew F. Pusey -- see LICENSE.txt

import cython
from bitfuncs import val2nice

cdef class Triplet:
    """Represents a conditional indepdence relation"""
    cdef int A
    cdef int B
    cdef int C
    cdef int myhash

    def __init__(self, int A, int B, int C):
        """Return < A , B | C > for bitstring A, B, C"""
        self.A = A
        self.B = B
        self.C = C
        self.myhash = hash((A, B, C))

    cdef Triplet sym(Triplet self):
        """Swap A and B"""
        return self.__class__(self.B, self.A, self.C);

    def __repr__(self):
        return "< %s , %s | %s >" % (val2nice(self.A),
            val2nice(self.B), val2nice(self.C))

    def __mul__(Triplet u, Triplet v):
        """* operator according to Definition 8"""
        cdef int A = u.A, B = u.B, C = u.C
        cdef int I = v.A, J = v.B, K = v.C
        cdef int ABC = A|B|C, IJK = I|J|K

        if (C & ~IJK) or (K & ~(ABC)):
            return None
        if ((A & I) == 0) or (((J & ~C) | (B & IJK)) == 0):
            return None

        return u.__class__(A & I, (J & ~C) | (B & IJK), C | (A & K))

    def __hash__(self):
        return self.myhash

    def __richcmp__(Triplet u, Triplet v, int op):
        cdef int X = u.A, Y = u.B, Z = u.C
        cdef int A = v.A, B = v.B, C = v.C
        # Dominating relation according to definition 7
        if op == 0:
            return ((X & ~A) == 0) and ((Y & ~B) == 0) and \
                ((C & ~Z) == 0) and ((Z & ~(A|B|C)) == 0)
        # Equality
        if op == 2:
            return X == A and Y == B and Z == C
        # Inequality
        if op == 3:
            return X != A or Y != B or Z != C

def closure(input_trips):
    """Dominating set for semigraphoid closure of Triples input_trips,
    according to Theorem 1"""
    trips = set(input_trips)
    old = set()
    cdef Triplet u, v, ustv
    while old != trips:
        old = trips.copy()

        for u in old:
            trips.add(u.sym())

        cpy = trips.copy()
        for u in cpy:
            for v in cpy:
                ustv = u * v;
                if ustv is not None:
                    trips.add(ustv)

        cpy = trips.copy()
        for u in cpy:
            for v in cpy:
                if u < v and u != v and v in trips:
                    trips.discard(u)

    return trips
