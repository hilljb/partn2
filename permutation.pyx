#!/usr/bin/python

# Partn2 - A permutation group and partition backtrack Python module
# permutation.pyx - permutation and permutation group base classes


###################################################################################################
# Copyright and license                                                                           #
###################################################################################################

# Copyright (C) 2013-2014 Jason B. Hill <jason@jasonbhill.com> <www.github.com/hilljb/partn2>

# This program is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details. See: <http://www.gnu.org/licenses/>.


###################################################################################################
# Imports                                                                                         #
###################################################################################################

import cython

from libc.stdlib cimport malloc, free
from libc.stdint cimport uint32_t


###################################################################################################
# Imports                                                                                         #
###################################################################################################

cdef class Permutation:
    """
    A permutation represented as an action on a set of natural numbers.
    """
    # A permutation $\sigma:i\mapsto\sigma(i)$ is recorded as a list [$\sigma(i)$]. Hence, the
    # permutation (1,2,5,6) acting on the set {1,2,3,4,5,6} is recorded as [2,5,3,4,6,1].
    cdef uint32_t              *permutation

    # The domain of a permutation is the set on which the permutation acts. This is consistent with
    # existing literature and programs such as GAP, but not entirely consistent with permutations
    # and permutation groups in Sage. This will be discussed more in the documentation.
    cdef public                 domain # implemented as a Python set of Python ints and longs

    # The degree of a permutation is the size of the domain on which the permutation acts.
    cdef public uint32_t               degree

    # self.maxSupport is the maximum point in self.domain under the usual ordering.
    cdef uint32_t               maxSupport

    # self.order is the order of the permutation
    cdef public uint32_t        order

    def __cinit__(self, perm, verbose=False):
        cdef uint32_t   i_32

        # initialize a set for the domain
        self.domain = set([])

        # perm is a list
        # case 1: perm is a list of disjoint tuples (cycles) of positive integers
        # example: Permutation([(1,2),(4,5)])
        if isinstance(perm, list) and all(isinstance(t, tuple) for t in perm):
            # verify that tuples have the correct structure
            for t in perm:
                for p in t:
                    # verify every point is an integer
                    if not isinstance(p, int) and not isinstance(p, long):
                        raise TypeError("\'%s\' not an int or long in tuple %s" % (p, t))
                    # verify every point is greater than 0
                    if p < 1:
                        raise ValueError("\'%s\' not greater than zero in tuple %s" % (p, t))
                    # verify every point is unique (integers cannot appear twice in cycle notation)
                    if p in self.domain:
                        raise ValueError("cycles not disjoint: %s appears at least twice")
                    else:
                        self.domain.add(p)

            # get the maximum point in the support of the permutation
            self.maxSupport = max(self.domain)

            # get the degree of the permutation
            self.degree = len(self.domain)

            # form the permutation map
            permutation = <uint32_t *>malloc((self.maxSupport + 1) * sizeof(uint32_t))
            # set to the identity
            for i_32 in range(self.maxSupport + 1):
                permutation[i_32] = i_32
            # now set any non-identity action
            for t in perm:
                for i in range(len(t)-1):
                    permutation[t[i]] = <uint32_t>t[i+1]
                permutation[t[len(t)-1]] = <uint32_t>t[0]

            s = '['
            for i in range(1, self.maxSupport):
                s += '%s,' % permutation[i]
            s += '%s]' % permutation[self.maxSupport]
            print s











