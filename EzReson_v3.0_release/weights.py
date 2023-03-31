# Weights of constituting wave functions
#
#  Created on Apr 9, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    Feb 18, 2022:
#    - Added P2WSO scheme for the WFRT weight is approx. proportional to the 
#    square of WF-projection, not the projection.
#    More importantly, for open-shell systems, the projection of resonance 
#    structure onto the actual structure can be negative. Thus, the squared
#    form seems more reasonable and natural.
#
#    Jul 23, 2020:
#    - Improved the PWSO scheme to avoid the problem that a pure covalent
#    actual wave function for a 2c-2e system has a covalent Lewis structure
#    contribution less than 100% (which is also a problem for Lowdin weights)
#

import numpy as np
from util import *

# Tolerance for validations:
TOL = 1E-5


# Mulliken weights:
#   C: Coefficient vector
#   S: Overlap matrix
def weightMulliken( C, S ):
    W = C * S @ C
    if abs( W.sum() - 1. ) > TOL:
        raise ValueError( 'Sum of weights (%E) greater than tolerence (%E) in '
                'weightMulliken' % ( abs( W.sum() - 1. ), TOL ) )
    return W
# enddef weightMulliken()


# Lowdin weights:
#   C: Coefficient vector
#   S: Overlap matrix
def weightLowdin( C, S ):
    Sd, U = np.linalg.eigh( S )
    Q = U @ np.diag( 1/np.sqrt(Sd) ) @ U.T
    C_Lowd = np.linalg.pinv( Q ) @ C
    W = C_Lowd * C_Lowd
    if abs( W.sum() - 1. ) > TOL:
        raise ValueError( 'Sum of weights (%E) greater than tolerence (%E) in '
                'weightLowdin' % ( abs( W.sum() - 1. ), TOL ) )
    return W
# enddef weightLowdin()


# Bickelhaupt weights:
#   C: Coefficient vector
#   S: Overlap matrix
def weightBickelhaupt( C, S ):
    C2 = C * C # Squares of coeff.
    N = len( C )
    W = np.zeros( ( N, 1 ) )
    for i in range( N ):
        if C[i] == 0: # To avoid singularity
            continue
        for j in range( N ):
            if C[j] == 0: # To avoid singularity
                continue
            W[i] += 2*C2[i]/(C2[i]+C2[j]) * C[i]*C[j]*S[i,j]
    if abs( W.sum() - 1. ) > TOL:
        raise ValueError( 'Sum of weights (%E) greater than tolerence (%E) in '
                'weightBickelhaupt' % ( abs( W.sum() - 1. ), TOL ) )
    return W
# enddef weightBickelhaupt()


# Ros-Schuit weights:
#   C: Coefficient vector
def weightRos_Schuit( C ):
    C2 = C * C # Squares of coeff.
    W = C2 / C2.sum()
    if abs( W.sum() - 1. ) > TOL:
        raise ValueError( 'Sum of weights (%E) greater than tolerence (%E) in '
                'weightRos_Schuit' % ( abs( W.sum() - 1. ), TOL ) )
    return W
# enddef weightRos_Schuit()



# THIS IS AN OLD VERSION:
#     The projection is directly used for the weighting in SO.
#     The problem is that, e.g. a 2c-2e system, if the acutal wave function is 
#   pure covalent, the covalent resonance structure is not 100% (as in
#   Mulliken weights), but less than 100% with small contribution from both 
#   ionic resonance structures.
#
# Projection-Weighted Symmetric Orthogonalization (PWSO) weights:
#   C: Coefficient vector
#   S: Overlap matrix
#   P: Projection vector
def weightPWSO_old( C, S, P ):
    a = 1E-4 # To avoid component with zero projection leading to singularity
    Pd = np.diagflat( (P+a) ) # Diagonal matrix with projections
    #Pd = np.eye( len(P) ) # Diagonal matrix with projections
    E, U = np.linalg.eigh( Pd @ S @ Pd )
    eps = 1E-15
    Q = Pd @ U @ np.diag( 1/np.sqrt(E+eps) ) @ np.linalg.pinv( U )
    #Snew = Q.T @ S @ Q
    #writeMat( Snew )
    C_orth = np.linalg.pinv( Q ) @ C
    W = C_orth * C_orth
    if abs( W.sum() - 1. ) > TOL:
        raise ValueError( 'Sum of weights (%E) greater than tolerence (%E) in '
                'weightPWSO' % ( abs( W.sum() - 1. ), TOL ) )
    return W
# enddef weightLowdin()


# Projection-Weighted Symmetric Orthogonalization (PWSO) weights:
#   C: Coefficient vector
#   S: Overlap matrix
#   P: Projection vector
def weightPWSO( C, S, P ):
    a = 1E-8 # To avoid component with zero projection leading to singularity

    # Regulate projections so that when P_i == 1, all P_j (j!=i) are zero:
    P = 1/(1-P-a)

    Pd = np.diagflat( (P+a) ) # Diagonal matrix with projections
    #Pd = np.eye( len(P) ) # Diagonal matrix with projections
    E, U = np.linalg.eigh( Pd @ S @ Pd )
    eps = 1E-15
    Q = Pd @ U @ np.diag( 1/np.sqrt(E+eps) ) @ np.linalg.pinv( U )
    #Snew = Q.T @ S @ Q
    #writeMat( Snew )
    C_orth = np.linalg.pinv( Q ) @ C
    W = C_orth * C_orth
    if abs( W.sum() - 1. ) > TOL:
        raise ValueError( 'Sum of weights (%E) greater than tolerence (%E) in '
                'weightPWSO' % ( abs( W.sum() - 1. ), TOL ) )
    return W
# enddef weightPWSO()


# Squared-Projection-Weighted Symmetric Orthogonalization (P2WSO) weights:
#   C: Coefficient vector
#   S: Overlap matrix
#   P: Projection vector
def weightP2WSO( C, S, P ):
    a = 1E-8 # To avoid component with zero projection leading to singularity

    # Regulate projections so that when P_i == 1, all P_j (j!=i) are zero:
    P = 1 / ( 1 - (P**2) - a )

    Pd = np.diagflat( (P+a) ) # Diagonal matrix with projections
    #Pd = np.eye( len(P) ) # Diagonal matrix with projections
    E, U = np.linalg.eigh( Pd @ S @ Pd )
    eps = 1E-15
    Q = Pd @ U @ np.diag( 1/np.sqrt(E+eps) ) @ np.linalg.pinv( U )
    #Snew = Q.T @ S @ Q
    #writeMat( Snew )
    C_orth = np.linalg.pinv( Q ) @ C
    W = C_orth * C_orth
    if abs( W.sum() - 1. ) > TOL:
        raise ValueError( 'Sum of weights (%E) greater than tolerence (%E) in '
                'weightP2WSO' % ( abs( W.sum() - 1. ), TOL ) )
    return W
# enddef weightP2WSO()
