# Topology and graph calculations
#
#  Created on Apr 23, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    May 31, 2020:
#    - Fixed a bug in function xyzNonH() for choosing non-hydrogen atoms
#

import numpy as np

# Only choose non-hydrogen atoms from a set of xyz coordinates:
def xyzNonH( elem, xyz ):
    at = []
    iat = -1
    for e in elem:
        iat += 1
        ecap = e.upper()
        if ecap == '1' or ecap == 'H':
            continue
        at.append( iat )
    xyz1 = xyz[ at, : ]
    return xyz1, at
# enddef xyzNonH()


# Adjacency matrix from a set of xyz coordinates:
#   xyz: a ndarray storing a set of Cartesian coordinates
def adjMat( xyz ):

    MAX_BONDLEN = 1.85;

    N = xyz.shape[0] # Number of atoms
    A = np.zeros( ( N, N ), dtype='i' )

    for i in range( N ):
        for j in range( i+1, N ):
            r = np.linalg.norm( xyz[i,:] - xyz[j,:] )
            if r <= MAX_BONDLEN:
                A[i,j] = 1
                A[j,i] = 1

    return A
# enddef adjMat()


# Get bond list from an adjacent matrix
#   A: Adjacent matrix
# Return:
#   List, each entry indicating indices (starting from 0) of the two atoms
def bonds( A ):

    N = len( A )
    B = []
    for i in range( N ):
        for j in range( i+1, N ):
            if A[i,j] == 1:
                B.append( np.array([ i, j ]) )

    return B
# enddef bonds()


# Remove all bonds involving atoms in bnd[] from bond list:
def removeBonds( B, bnd ):

    # Number of bonds in the list:
    NB = len( B )

    Bnew = []
    for i in range(NB):
        if B[i][0] == bnd[0] or B[i][1] == bnd[0] or \
           B[i][0] == bnd[1] or B[i][1] == bnd[1]:
            continue
        Bnew.append( B[i] )

    return Bnew
# enddef removeBonds()
