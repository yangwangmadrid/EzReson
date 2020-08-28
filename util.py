# General utilities
#
#  Created on Feb, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#

import numpy as np

#==============================================================================
# Conversion of units
def au2kcal( E ):
    return np.array( E ) * 627.50947427719404

#==============================================================================


#==============================================================================
# Write a string
#
#   fileName: name of the file to write. Tips: it can be "/dev/stdout"
#   M: the matrix
#   strFmt: string format
#   mode: 'w' or 'a'
def writeStr( string, fileName='/dev/stdout', end='', strFmt='%s', mode='w' ):
    with open( fileName, mode ) as writer:    
        writer.write( strFmt % string )
        writer.write( end )
#==============================================================================
# enddef writeStr()

#==============================================================================
# Write a matrix
#
#   fileName: name of the file to write. Tips: it can be "/dev/stdout"
#   M: the matrix
#   valFmt: string format for the written values
#   mode: 'w' or 'a'
def writeMat( M, fileName='/dev/stdout', valFmt=' %5.2f', mode='w' ):

    with open( fileName, mode ) as writer:    
        for row in M:
            if np.isscalar( row ):
                row = np.array( [ row ] )
            for x in row:
                writer.write( valFmt % x )
            writer.write( '\n' )


#==============================================================================
# enddef writeMat()


# Split a given array into degenerate groups, each of which contains
# identical (within a given tolerence) values
#
# Return:  A list containg arrays of indices corresponding to degenerate groups 
#
def degenerateGroups( arr, tol ):
    N = len( arr )
    iG = 0 # Current index of group
    G = [ [0] ] # A list containing groups of indices
    a_last = arr[0] # Last element in arr[]

    for j in range( 1, N ):
        if abs( arr[j] - a_last ) < tol:
            G[iG].append( j )
        else:
            G.append( [j] )
            iG += 1
            a_last = arr[j]

    return G

# enddef degenerateGroups()

# Get the maximum eigenvalue and eigenvector of a given matrix:
# NOTE: ONLY returns ONE eigenvector even for degenerate cases
# NOTE: This function is only valid for REAL eigenvalues/eigenvectors
#       So, make sure that the input A is a REAL diagonalizable matrix
def maxEig( A ):
    if A.shape[0] != A.shape[1]:
        raise ValueError( 'The input matrix must be a square matrix' )
    #U, E, _ = np.linalg.svd( A )
    E, U = np.linalg.eigh( A )
    E = E.real # Only take the real part
    U = U.real # Only take the real part
    # Sort in descending order:
    ix = np.argsort( -E )
    E = E[ ix ]
    U = U[ :, ix ]
    return ( U[ :, [0] ], E[0] )
# enddef maxEig()

# Get all eigenvalues (in descending order) and eigenvectors of a given matrix:
# NOTE: This function is only valid for REAL eigenvalues/eigenvectors
#       So, make sure that the input A is a REAL diagonalizable matrix
def eigen( A, order='descend' ):
    if A.shape[0] != A.shape[1]:
        raise ValueError( 'The input matrix must be a square matrix' )
    E, U = np.linalg.eigh( A )
    E = E.real # Only take the real part
    U = U.real # Only take the real part
    if order == 'ascend':
        ix = np.argsort( E )
        E = E[ ix ]
        U = U[ :, ix ]
    elif order == 'descend':
        ix = np.argsort( -E )
        E = E[ ix ]
        U = U[ :, ix ]
    else:
        raise ValueError( 'Invalid value (%s) for argument order' % order )
    return ( U, E )
# enddef eigen()

# Get all nonzero (within a rounding error, tol) eigenvalues and eigenvectors 
# of a given matrix A:
def nonzeroEig( A, tol ):
    E, U = np.linalg.eigh( A )
    E = E.real # Only take the real part
    U = U.real # Only take the real part
    ix = np.nonzero( abs(E) > tol )[0]
    return ( U[:,ix], E[ix] )
# enddef nonzeroEig()

# Determine if matrix is diagonal or not:
#   Return True if the input is a diagonal matrix
def isdiagonal( A, tol ):
    if A.shape[0] != A.shape[1]:
        raise ValueError( 'The input matrix must be a square matrix' )
    dev = np.linalg.norm( A - np.diag( np.diag(A) ) ) / len(A) / len(A)
    if dev < tol:
        result = True
    else:
        result = False
    return ( result, dev )
# enddef isdiagonal()

