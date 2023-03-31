# Topology and graph calculations
#
#  Created on Apr 23, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    Feb 4, 2021:
#    - Added function readRings()
#
#    Dec 27, 2020:
#    - Added functions to get connectivity and neighbors of rings
#
#    Dec 26, 2020:
#    - Added functions to determine the rings
#
#    May 31, 2020:
#    - Fixed a bug in function xyzNonH() for choosing non-hydrogen atoms
#

import numpy as np

# Only choose non-hydrogen atoms from a set of xyz coordinates:
#
#  Return:
#           xyz1: Cartesian coordinates of non-hydrogen atoms
#             at: List of indices of atoms ( starting from 0 )
#
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
#
#    A: Adjacent matrix
#   at: List of indices of atoms (starting from 0)
#
# Return:
#   List, each entry indicating indices (starting from 0) of the two atoms
def bonds( A, at=[] ):

    N = len( A )
    B = []
    for i in range( N ):
        for j in range( i+1, N ):
            if A[i,j] == 1:
                B.append( np.array([ i, j ]) )

    # Assign back the original atomic labels in B[]:
    if len( at ) == 0:
        return B

    B0 = []
    for bnd in B:
        B0.append( np.array([ at[bnd[0]], at[bnd[1]] ]) )
    return B0

# enddef bonds()


# Get neighbors of each atom
#
#    A: Adjacent matrix
#   at: List of indices of atoms (starting from 0)
#
# Return:
#   List, each entry indicating indices (starting from 0) of all neighbors of 
#   each atom
def neighbors( A, at=[] ):
    N = len( A )
    nblist = []
    for i in range( N ):
        ix = [ j for j, a in enumerate(A[i,:]) if a > 0 ] 
        nblist.append( ix )

    # Assign back the original atomic labels in B[]:
    if len( at ) == 0:
        return nblist

    nblist0 = []
    for nb in nblist:
        x = [ at[a] for a in nb ]
        nblist0.append( x )
    return nblist0
# enddef neighbors()


# Remove EXACTLY the bonds specified by list bnd[] from the whole bond list B[]:
# NOTE: The two atomic indices in bnd[] are in ascending order, and so are the
# two atomic indices in each element of list B[]
def removeExactBonds( B, bnd ):

    Bnew = []
    for bd in B:
        if np.array_equal( bd, bnd ):
            continue
        Bnew.append( bd )

    return Bnew
# enddef removeExactBonds()


# Remove all bonds INVOLVING atoms in bnd[] from bond list:
# NOTE: Note the difference between removeBonds() and removeExactBonds()
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


# Remove all bonds INVOLVING given atoms:
def removeBonds_by_Atoms( B, Atoms ):
 
    # Number of bonds in the list:
    NB = len( B )
    # Number of given atoms:
    if not isinstance( Atoms, list ):
        Atoms = [ Atoms ]
    NA = len( Atoms )
 
    Bnew = []
    for i in range( NB ):
        flag = False
        for j in range( NA ):
            # Found the bond involving atom j:
            if B[i][0] == Atoms[j] or B[i][1] == Atoms[j]:
                flag = True
                break
        if flag:
            continue
        Bnew.append( B[i] )
 
    return Bnew
# enddef removeBonds_by_Atoms()


# Remove all involved atoms in bnd[] from atom list:
def removeAtoms_from_Bond( A, bnd ):
 
    # Number of atoms in the list:
    NAt = len( A )
 
    Anew = []
    for i in range(NAt):
        if A[i] == bnd[0] or A[i] == bnd[1]:
            continue
        Anew.append( A[i] )
 
    return Anew
# enddef removeAtoms_from_Bonds()


# Check if a ring is monocyclic
#
#  rg: a ring specified by a list of indices of atoms (starting from 0)
#  rglist: List of rings to store the results of determination; For initial
#          call, set it an empty list
def isMonocyclic( rg, nblist ):
    nrg = len( rg )

    for j in range( nrg ):
        iat = rg[j]
        if j == 0:
            iat_1 = rg[nrg-1]
            iat1 = rg[1]
        elif j == nrg-1:
            iat_1 = rg[j-1]
            iat1 = rg[0]
        else:
            iat_1 = rg[j-1]
            iat1 = rg[j+1]
        nb = nblist[ iat ]
    
        for a in nb:
            if a == iat_1 or a == iat1:
                continue
            # The rest of neighbors of iat:
            if a in rg:
                return False

    return True
# enddef isMonocyclic()


# Canonical numbering of a ring. The rule is as follows:
#   - 1. Keep the first two numbers as small as possible.
#   - 2. Being clockwise or counterclockwise does not matter.
def canonicalRing( rg ):
    nrg = len( rg )
    ix = rg.index( min(rg) )
    if ix == 0:
        at_1 = rg[ nrg-1 ]
    else:
        at_1 = rg[ ix-1 ]

    if ix == nrg-1:
        at1 = rg[0]
    else:
        at1 = rg[ ix+1 ]

    num = list( range( nrg ) )
    loop = num[ ix+1: ] + num[ 0 : ix ]
    if at_1 < at1: # Reverse
        loop.reverse()

    loop = [ ix ] + loop

    return [ rg[i] for i in loop ]
# enddef canonicalRing()


# Get rings from a given atom
#
#  rglist: List of rings to store the results of determination; For initial
#          call, set it an empty list
#  rg: Current ring; For initial call, rg = [iat], where iat is the index
#      of the given atom
#  iat: Index of the given atom, which must be from 0, 1, 2, ..., NAt
#       NOTE: It may not the actual atomic index in the whole molecule that
#       contains hydrogens; it is the appearance index among non-H atoms
#  nblist: List, each entry indicating indices (starting from 0) of all 
#          neighbors of each atom
#
# Return:
#   List, each entry indicating indices (starting from 0) of each ring
def rings_from_an_atom( rglist, rg, nblist, iat ):

    # Number of atoms in current rg[]:
    nrg = len( rg )
    nb = nblist[ iat ]

    for at_next in nb:
        # Check if at_next is already in rg[]:
        #ix = [ i for i, a in enumerate(rg) if a == at_next ]
        #ix = rg.index( at_next )
        rg_curr = rg.copy()

        # Avoid large rings:
        if len( rg_curr ) > 7:
            continue

        if at_next in rg:
            #print( at_next, rg )
            #input()
            ix = rg.index( at_next )
            # Backward atom:
            if ix == nrg-2:
                continue
            # Ring complete:
            elif rg.index( at_next ) == 0:
                if isMonocyclic( rg_curr, nblist ):
                    rg_curr = canonicalRing( rg_curr )
                    # Check if rg is a new ring:
                    ifNew = True
                    for rg_k in rglist:
                        if rg_curr == rg_k:
                            ifNew = False
                            break
                    if ifNew:
                        rglist.append( rg_curr )
                        # print( rglist )

            continue
        else:
            # Ring growth:
            rg_curr.append( at_next )

        rings_from_an_atom( rglist, rg_curr, nblist, at_next )

    return rglist
# enddef rings_from_an_atom()


# Get ring list of a molecular graph
#
#  nblist: List, each entry indicating indices (starting from 0) of all 
#          neighbors of each atom
#  at: List of indices of atoms (starting from 0)
#
# Return:
#   List, each entry indicating indices (starting from 0) of atoms in each ring
def rings( nblist, at=[] ):
    NAt = len( nblist )

    rglist = []
    for iat in range( NAt ):
        rglist = rings_from_an_atom( rglist, [ iat ], nblist, iat )

    # Assign back the original atomic labels in B[]:
    if len( at ) == 0:
        return rglist

    rglist0 = []
    for rg in rglist:
        x = [ at[a] for a in rg ]
        rglist0.append( x )
    return rglist0
# enddef rings()


# Determine if two rings are adjacent
#
def ifAdjacetRings( rg1, rg2 ):
    nrg1 = len( rg1 )
    nrg2 = len( rg2 )
    for i in range( nrg1 ):
        if i < nrg1-1:
            bnd1 = [ rg1[i], rg1[i+1] ]
        else:
            bnd1 = [ rg1[i], rg1[0] ]
        # Reorder the two atom labels in bnd1:
        if bnd1[1] < bnd1[0]:
            bnd1 = [ bnd1[1], bnd1[0] ]

        for j in range( nrg2 ):
            if j < nrg2-1:
                bnd2 = [ rg2[j], rg2[j+1] ]
            else:
                bnd2 = [ rg2[j], rg2[0] ]
            # Reorder the two atom labels in bnd2:
            if bnd2[1] < bnd2[0]:
                bnd2 = [ bnd2[1], bnd2[0] ]

            if bnd1 == bnd2:
                return True

    return False
# enddef ifAdjacetRings()


# Get neighbors of each ring
#
#   rglist: List of indices of atoms (starting from 0) of each ring
#
# Return:
#   List, each entry indicating indices (starting from 0) of all neighbors of 
#   each ring
def ringNeighbors( rglist ):
    N = len( rglist )

    # Get connectivity between rings:
    A = np.zeros( ( N, N ), dtype='i' )
    for i in range( N ):
        for j in range( i+1, N ):
            A[ i, j ] = ifAdjacetRings( rglist[i], rglist[j] )
            A[ j, i ] = A[ i, j ]
    #print( A )

    nblist = []
    for i in range( N ):
        ix = [ j for j, a in enumerate(A[i,:]) if a > 0 ] 
        nblist.append( ix )
    return nblist
# enddef RingNeighbors()


# Read rings from *.rings file and return them to RG[]
#   ringsFileName: Name of *.rings file where atomic indices of rings are stored
def readRings( ringsFileName ):
    with open( ringsFileName ) as f:
        RG = [ [int(x) for x in line.split()] for line in f ]

    return RG
# enddef readRings()
