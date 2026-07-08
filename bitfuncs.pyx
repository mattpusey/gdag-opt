"""Helper functions for dealing with bitstrings"""

# Copyright (C) 2014 Matthew F. Pusey
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.

import string

def getindex(int bit):
    """Get the position of the lowest 1 in bit"""
    cdef int i
    for i in range(26):
        if bit & 1<<i:
            return i

def val2nice(int val):
    """Output a bitstring val as letters, e.g. 101 becomes 'AC'"""
    return "".join([string.ascii_uppercase[i] for i in bitsof(val)])

def nice2val(nice):
    """Turn letters nice into a bitstring, e.g. 'AC' becomes 101"""
    cdef int val = 0
    for i in nice:
        val |= 1<<string.ascii_uppercase.find(i)
    return val

def bitsof(int left):
    """Generator of the locations of 1s in bitstring left"""
    cdef int x
    while left:
        x = left & ~(left - 1)
        left &= ~x
        yield getindex(x)

def nonemptysubsetsof(int space):
    """Generator of the non-empty subsets of bitsring space"""
    cdef int x = space
    while x:
        yield x
        x = (x-1) & space

def subsetsof(int space):
    """Generator of the subsets of bitsring space"""
    cdef int x
    for x in nonemptysubsetsof(space):
        yield x
    yield 0

def propersubsetsof(int space):
    """Generator of the proper subsets of bitsring space"""
    cdef int x = space
    while x:
        x = (x-1) & space
        yield x
