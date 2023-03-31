# Matrix and vector transformations
#
#  Created on Feb, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    Jun 12, 2020:
#    - Allowed automatically loading or dumping LMOs for spin-unrestricted cases
#    in function genLMOs()
#

import numpy as np
import os.path
from writeFchkOrb import *
from scipy.linalg import *
from math import *
from util import *

# Density matrix calculated from CMO coefficient matrix:
# Note:
# For spin-restricted:
#   C is a single coefficient matrix and NE is the total number of electrons
# For spin-unrestricted:
#   C is a tuple: C = ( Ca, Cb )
#   NE is a tuple: NE = ( NEa, NEb )
#
def densityMatrix( C, NE ):
    # Determine whether this is spin-unrestricted or not:
    if type(C) is tuple:
        spin = 1 # spin-unrestricted
        Ca, Cb = C
        NEa, NEb = NE
        NBas, NMO = Ca.shape
        Da = np.zeros( (NBas, NBas) )
        Db = np.zeros( (NBas, NBas) )
    else:
        spin = 0 # spin-restricted
        Ca = C
        if NE % 2 == 1:
            raise ValueError( 'Odd number of electrons (%i) is found '
                    'for a spin-restricted calculation' % NE )
        NEa = NE//2
        ( NBas, NMO ) = Ca.shape
        Da = np.zeros( (NBas, NBas) )

    # A faster to compute DM:
    for i in range( NEa ):
        # Note: Numpy 1D array has no transpose. Use outer() to get tensor
        Da += np.outer( Ca[:,i], Ca[:,i] ) #Note: Numpy 1Darray has no transpose
    # For spin-unrestricted, get beta DM as well:
    if spin == 1:
        for i in range( NEb ):
            # Note: Numpy 1D array has no transpose. Use outer() to get tensor
            Db += np.outer( Cb[:,i], Cb[:,i] )

    # Return:
    if spin == 0:
        return 2*Da
    else:
        return ( Da, Db )
# enddef densityMatrix()


# Overlap matrix in AO basis calculated from CMO coefficient matrix C:
# NOTE: Only if no linear dependence is present in C (i.e. C is invertible)
def overlapMatrix( C ):
    # Determine whether this is spin-unrestricted or not:
    if type(C) is tuple:
        spin = 1 # spin-unrestricted
        Ca, Cb = C
        NBas, NMO = Ca.shape
    else:
        spin = 0 # spin-restricted
        Ca = C
        NBas, NMO = Ca.shape

    # Check linear dependence in the basis:
    if NBas != NMO:
        raise ValueError( 'Linear dependence in the basis is found. '
                'The overlap matrix must be explicitly calculated.' )

    # For alpha electrons:
    invC = inv( Ca )
    Sa = invC.T @ invC
    # writeMat( Sa, '/dev/stdout', ' %12.6f' )

    # For spin-unrestricted, get Sb for beta electrons as well:
    if spin == 1:
        invC = inv( Cb )
        Sb = invC.T @ invC

    # Return:
    if spin == 0:
        return Sa
    else:
        return ( Sa + Sb )/2
# enddef overlapMatrix()


# Fock matrix in AO basis calculated from CMO coefficient matrix C and CMO
# energies:
# NOTE: Only if no linear dependence is present in C (i.e. C is invertible)
def FockMatrix( C, E ):
    # Determine whether this is spin-unrestricted or not:
    if type(C) is tuple:
        spin = 1 # spin-unrestricted
        Ca, Cb = C
        NBas, NMO = Ca.shape
        Ea, Eb = E
    else:
        spin = 0 # spin-restricted
        Ca = C
        NBas, NMO = Ca.shape
        Ea = E

    # Check linear dependence in the basis:
    ifBasDep = False
    if NBas != NMO:
        #raise ValueError( 'Linear dependence in the basis is found. '
        #        'The Fock matrix must be explicitly calculated.' )
        print( 'WARNING: Linear dependence in the basis is found. '
               'The Fock matrix ought to be explicitly calculated.' )
        ifBasDep = True

    # For alpha electrons:
    if not ifBasDep:
        invC = inv( Ca )
    else:
        invC = pinv( Ca )
#        invC = linalg.pinv( Ca )
    Fa = invC.T @ np.diag(Ea) @ invC
    # writeMat( FAOa, '/dev/stdout', ' %12.6f' )

    # For spin-unrestricted, get Sb for beta electrons as well:
    if spin == 1:
        invC = inv( Cb )
        Fb = invC.T @ np.diag(Eb) @ invC

    # Return:
    if spin == 0:
        return Fa
    else:
        return ( Fa, Fb )
# enddef FockMatrix()

# Atomic population calculated from a denisty matrix in ORTHOGONAL basis
def atomicPopulation( D, orbIx ):
    # For spin-unrestricted:
    if type(D) is tuple:
        Dd = np.diag(D[0]) + np.diag(D[1]) # Merge alpha + beta densities
    # For spin-restricted:
    else:
        Dd = np.diag( D )
    q = np.zeros( len(orbIx) )

    iat = 0
    for ix in orbIx:
        q[iat] = sum( Dd[ix-1] )
        iat += 1

    return q
# enddef atomicPopulation()


# Atomic population calculated from a denisty matrix in NON-ORTHOGONAL basis,
# based on Mulliken's scheme
def atomicPopulation_nonortho( D, S, orbIx ):
    # For spin-unrestricted:
    if type(D) is tuple:
        D = D[0] + D[1] # Merge alpha + beta densities
    return atomicPopulation( D@S, orbIx )
# enddef atomicPopulation_nonortho()


# Atomic charges calculated from a denisty matrix in ORTHOGONAL basis
def atomicCharge( D, orbIx, Zeff ):
    return Zeff - atomicPopulation( D, orbIx )
# enddef atomicCharge()


# Atomic charges calculated from a denisty matrix in NON-ORTHOGONAL basis,
# based on Mulliken's scheme
def atomicCharge_nonortho( D, S, orbIx, Zeff ):
    return Zeff - atomicPopulation_nonortho( D, S, orbIx )
# enddef atomicCharge_nonortho()


# Lowdin's symmetric orthogonalization
def LowdinOrtho( S ):
    NBas = len( S )
    # Check linear dependence of the original basis set:
    if np.linalg.matrix_rank(S) < NBas:
        raise NotImplementedError( 'Linear dependence is found in the basis. '
                'Lowdin orthogonalization aborted.' )

    # NOTE: Using svd to get eigenvalues is MUCH MORE accurate than using eig()
    # function: the accuracy (see the err variable below) by the svd method is 
    # MANY orders of magnitude higher than that by the eig method. And the 
    # former is also slightly faster than the latter. 
    # So, we use svd() instead of using, e.g., Sd, U = eig( S )
    U, Sd, _ = svd( S )
    #print( U @ U.T )
    Q = U @ np.diag( 1/np.sqrt(Sd) ) @ U.T

    # Verify the new basis is really an orthonormal set
    err = norm( Q.T @ S @ Q - np.eye(NBas) )/ NBas
    if err > 1e-10:
        raise AssertionError( 'Accuray of Lowdin orthogonalization is not ' 
        'satisfactory: (%E) \nAborted' % err )

    return Q
# enddef LowdinOrtho()


# Fast version of PipekMezeyLMO_1spin()
#
# Pipek-Mezey localization of molecular orbitals for a single spin (alpha or
# beta)
#   C: MO coefficient matrix in orthonormal AO basis (e.g., NAOs, Lowdin AOs)
#   F: Fock matrix in orthonormal AO basis (e.g., NAOs, Lowdin AOs)
#   aoIx: AO indices for each atom
#   NE: Number of electrons
#   maxInter: Maximum iterations to be done
def PipekMezeyLMO_1spin_fast( C, F, aoIx, NE, maxIter=9000 ):
    # Tolerance for convergence:
    TOL = 1E-10
    lowTOL = 1E-4 # Low accuracy for trial runs

    NBas, NMO = C.shape
    if NE > NMO:
        raise ValueError( 'Number of electrons (%i) must not be greater than '
                'the number of molecular orbitals (%i)' % (NE, NMO) )

    # Initiate variables:
    C_LMO = np.copy( C[:,0:NE] )

    # Array of contributions of MOs to the objective functional:
    f_arr = np.zeros( NE )
    for i in range( NE ):
        for ix in aoIx:
            # Atomic population of each atom
            # qA = sum( C_LMO[ix-1,i] * C_LMO[ix-1,i]  )
            # The following vectorization is preferred:
            qA = C_LMO[ix-1,i].T @ C_LMO[ix-1,i]
            f_arr[i] += qA * qA
    # Sort by f_arr[], in descending order (i.e. from more localized to less)
    ix_sort_f = np.argsort( -f_arr )

    # Generate initial pair list for indices (s,t) of MOs:
    npair = NE * (NE-1) // 2
    ix_pair = np.zeros( (npair, 2 ), dtype='i' )
    iPair = 0
    for ix_s in np.arange( NE ):
        s = ix_sort_f[ ix_s ]
        for ix_t in np.arange( ix_s+1, NE ):
            t = ix_sort_f[ ix_t ]
            ix_pair[ iPair, : ] = [ s, t ]
            iPair += 1
    ix_pair0 = np.copy( ix_pair )

    # Objective functional, D^{-1}, to be maximized:
    f0 = 0
    for iIter in range( maxIter ):
        print( 'Pipek-Mezey localization iteration %i ...' % (iIter+1) )

        # print( ix_pair.shape )
        iPair = 0
        for (s,t) in ix_pair: # a pair of MOs: s and t
            A = 0
            B = 0
            for ix in aoIx:
                # Integrals < s| P_iat | t > for atom iat:
                # For orthogonal basis, we have:
                # Pst = sum_{u in A}( c_u,s*c_u,t )
                # or in matrix form as:
                # Pst = C_s[rows in A]' @ C_t[:]
                # Note that vectorization can speed up significantly,
                # compared to element-wise summation.
                Pss = C_LMO[ix-1,s].T @ C_LMO[ix-1,s]
                Ptt = C_LMO[ix-1,t].T @ C_LMO[ix-1,t]
                Pst = C_LMO[ix-1,s].T @ C_LMO[ix-1,t]
                dPst = Pss - Ptt
                A += Pst * Pst  -  dPst * dPst /4 
                B += Pst * dPst

            # Determine the rotation angle:
            # Trivial case: A = B = 0, no rotation required
            if A == 0 and B == 0:
                continue
            # General case for maximization: (gamma == alpha => a)
            R = sqrt( A*A + B*B )
            a = atan2( B/R, -A/R ) / 4
            sina = sin(a)
            cosa = cos(a)
            # Remove already converged MO pairs:
            if abs(sina) < lowTOL*100:
                ix_pair = np.delete( ix_pair, iPair, 0 )
                continue
            Q = np.array([ [ cosa, -sina], [sina, cosa] ])
            # Do the rotation and update coefficient matrix C_LMO:
            C_LMO[:,(s,t)] = C_LMO[:,(s,t)] @ Q

            iPair += 1
        # endfor (s,t)

        # Calculate the objective functional:
        f = 0
        for i in range( NE ):
            for ix in aoIx:
                # Atomic population of each atom
                # qA = sum( C_LMO[ix-1,i] * C_LMO[ix-1,i]  )
                # The following vectorization is preferred:
                qA = C_LMO[ix-1,i].T @ C_LMO[ix-1,i]
                f += qA * qA

        # Check convergence:
        dev = ( f - f0 ) / NE
        if abs(dev) < TOL:
            print( 'Pipek-Mezey localization converged with an accuracy of '
                   '%E after %i iterations' % ( abs(dev), iIter+1 ) )
            break
        elif dev < 0:
            raise RuntimeError( 'Failure in Pipek-Mezey localization: '
                    'Unexpected decrease of the objective functional: '
                    '%E < %E' %  ( f, f0 ) )
        elif iIter == maxIter-1:
            raise RuntimeError( 'Failure in Pipek-Mezey localization: '
                    'maximum number of iteractions exceeded (%i) and the '
                    'convergence criterion is not met (%E > %E)' % 
                    (maxIter, abs(dev), TOL) )
        # Update f0:
        f0 = f
        print( '    Current objective function = %.12f' % f )

        # Refinement when the low accuracy is fulfilled:
        if abs(dev) < lowTOL:
            ix_pair = np.copy( ix_pair0 )
            lowTOL = TOL
    # endfor iIter


    # CMO --> LMO transformation matrix:
    Q = lstsq( C[:,0:NE], C_LMO )[0]
    # Check if the left inverse is correctly computed by lstsq():
    assert( norm( C[:,0:NE] @ Q - C_LMO ) / NE < 1E-14 )
    # Fock matrix in LMO basis:
    # NOTE: The Fock matrix of CMOs in AO basis is diagonal, which is equal to 
    #       np.diag( E[0:NE] ). However, in general cases, where either the CMO 
    #       coefficients are expressed in NAO/Lowdin basis, or the input MOs
    #       are not CMOs, the Fock matrix is NOT diagonal.
    F_MO = C[:,0:NE].T @ F @ C[:,0:NE]  # AO --> MO
    F_LMO = Q.T @ F_MO @ Q
    # Get LMO energies:
    E_LMO = np.diag( F_LMO )

    # Sort LMOs by energies (from lower to higher):
    ix = np.argsort( E_LMO )
    C_LMO = C_LMO[:,ix]
    E_LMO = E_LMO[ix]

    return ( C_LMO, E_LMO )
# enddef PipekMezeyLMO_1spin_fast()


# Pipek-Mezey localization of molecular orbitals for a single spin (alpha or
# beta)
#   C: MO coefficient matrix in orthonormal AO basis (e.g., NAOs, Lowdin AOs)
#   F: Fock matrix in orthonormal AO basis (e.g., NAOs, Lowdin AOs)
#   aoIx: AO indices for each atom
#   NE: Number of electrons
#   maxInter: Maximum iterations to be done
def PipekMezeyLMO_1spin( C, F, aoIx, NE, maxIter=9000 ):
    # Tolerance for convergence:
    TOL = 1E-10

    NBas, NMO = C.shape
    if NE > NMO:
        raise ValueError( 'Number of electrons (%i) must not be greater than '
                'the number of molecular orbitals (%i)' % (NE, NMO) )

    # Initiate variables:
    C_LMO = np.copy( C[:,0:NE] )

    # Array of contributions of MOs to the objective functional:
    f_arr = np.zeros( NE )
    for i in range( NE ):
        for ix in aoIx:
            # Atomic population of each atom
            # qA = sum( C_LMO[ix-1,i] * C_LMO[ix-1,i]  )
            # The following vectorization is preferred:
            qA = C_LMO[ix-1,i].T @ C_LMO[ix-1,i]
            f_arr[i] += qA * qA
    # Sort by f_arr[], in descending order (i.e. from more localized to less)
    ix_sort_f = np.argsort( -f_arr )

    # Objective functional, D^{-1}, to be maximized:
    f0 = 0
    for iIter in range( maxIter ):
        print( 'Pipek-Mezey localization iteration %i ...' % (iIter+1) )

        for ix_s in np.arange( NE ): # MO s
            s = ix_sort_f[ ix_s ]
            for ix_t in np.arange( ix_s+1, NE ): # MO t
                t = ix_sort_f[ ix_t ]

                A = 0
                B = 0
                for ix in aoIx:
                    # Integrals Pst = < s| P_iat | t > for atom iat:
                    # For orthogonal basis, we have:
                    # Pst = sum_{u in A}( c_u,s*c_u,t )
                    # or in matrix form as:
                    # Pst = C_s[rows in A]' @ C_t[:]
                    # Note that vectorization can speed up significantly,
                    # compared to element-wise summation.
                    Pss = C_LMO[ix-1,s].T @ C_LMO[ix-1,s]
                    Ptt = C_LMO[ix-1,t].T @ C_LMO[ix-1,t]
                    Pst = C_LMO[ix-1,s].T @ C_LMO[ix-1,t]
                    dPst = Pss - Ptt
                    A += Pst * Pst  -  dPst * dPst /4 
                    B += Pst * dPst
 
                # Determine the rotation angle:
                # Trivial case: A = B = 0, no rotation required
                if A == 0 and B == 0:
                    continue
                # General case for maximization: (gamma == alpha => a)
                R = sqrt( A*A + B*B )
                a = atan2( B/R, -A/R ) / 4
                sina = sin(a)
                cosa = cos(a)
                Q = np.array([ [ cosa, -sina], [sina, cosa] ])
                # Do the rotation and update coefficient matrix C_LMO:
                C_LMO[:,(s,t)] = C_LMO[:,(s,t)] @ Q
            # endfor t
        # endfor s

        # Calculate the objective functional:
        f = 0
        for i in range( NE ):
            for ix in aoIx:
                # Atomic population of each atom
                # qA = sum( C_LMO[ix-1,i] * C_LMO[ix-1,i]  )
                # The following vectorization is preferred:
                qA = C_LMO[ix-1,i].T @ C_LMO[ix-1,i]
                f += qA * qA

        # Check convergence:
        dev = ( f - f0 ) / NE
        if abs(dev) < TOL:
            print( 'Pipek-Mezey localization converged with an accuracy of '
                   '%E after %i iterations' % ( abs(dev), iIter+1 ) )
            break
        elif dev < 0:
            raise RuntimeError( 'Failure in Pipek-Mezey localization: '
                    'Unexpected decrease of the objective functional: '
                    '%E < %E' %  ( f, f0 ) )
        elif iIter == maxIter-1:
            raise RuntimeError( 'Failure in Pipek-Mezey localization: '
                    'maximum number of iteractions exceeded (%i) and the '
                    'convergence criterion is not met (%E > %E)' % 
                    (maxIter, abs(dev), TOL) )
        # Update f0:
        f0 = f
        print( '    Current objective function = %.12f' % f )
    # endfor iIter


    # CMO --> LMO transformation matrix:
    Q = lstsq( C[:,0:NE], C_LMO )[0]
    # Check if the left inverse is correctly computed by lstsq():
    assert( norm( C[:,0:NE] @ Q - C_LMO ) / NE < 1E-14 )
    # Fock matrix in LMO basis:
    # NOTE: The Fock matrix of CMOs in AO basis is diagonal, which is equal to 
    #       np.diag( E[0:NE] ). However, in general cases, where either the CMO 
    #       coefficients are expressed in NAO/Lowdin basis, or the input MOs
    #       are not CMOs, the Fock matrix is NOT diagonal.
    F_MO = C[:,0:NE].T @ F @ C[:,0:NE]  # AO --> MO
    F_LMO = Q.T @ F_MO @ Q
    # Get LMO energies:
    E_LMO = np.diag( F_LMO )

    # Sort LMOs by energies (from lower to higher):
    ix = np.argsort( E_LMO )
    C_LMO = C_LMO[:,ix]
    E_LMO = E_LMO[ix]

    return ( C_LMO, E_LMO )
# enddef PipekMezeyLMO_1spin()


# Fast version of PipekMezeyLMO_1spin_nonortho()
#
# Pipek-Mezey localization of molecular orbitals for a single spin (alpha or
# beta)
#   C: MO coefficient matrix in AO basis
#   S: Overlap matrix between AO basis
#   F: Fock matrix in AO basis
#   aoIx: AO indices for each atom
#   NE: Number of electrons
#   maxInter: Maximum iterations to be done
def PipekMezeyLMO_1spin_nonortho_fast( C, S, F, aoIx, NE, maxIter=9000 ):
    # Tolerance for convergence:
    TOL = 1E-10
    lowTOL = 1E-4 # Low accuracy for trial runs

    NBas, NMO = C.shape
    if NE > NMO:
        raise ValueError( 'Number of electrons (%i) must not be greater than '
                'the number of molecular orbitals (%i)' % (NE, NMO) )

    # Initiate variables:
    C_LMO = np.copy( C[:,0:NE] )

    # Array of contributions of MOs to the objective functional:
    f_arr = np.zeros( NE )
    for i in range( NE ):
        for ix in aoIx:
            # Atomic population of each atom
            # The following vectorization is preferred:
            qA = C_LMO[ix-1,i].T @ S[ix-1,:] @ C_LMO[:,i]
            f_arr[i] += qA * qA
    # Sort by f_arr[], in descending order (i.e. from more localized to less)
    ix_sort_f = np.argsort( -f_arr )

    # Generate initial pair list for indices (s,t) of MOs:
    npair = NE * (NE-1) // 2
    ix_pair = np.zeros( (npair, 2 ), dtype='i' )
    iPair = 0
    for ix_s in np.arange( NE ):
        s = ix_sort_f[ ix_s ]
        for ix_t in np.arange( ix_s+1, NE ):
            t = ix_sort_f[ ix_t ]
            ix_pair[ iPair, : ] = [ s, t ]
            iPair += 1
    ix_pair0 = np.copy( ix_pair )

    # Objective functional, D^{-1}, to be maximized:
    f0 = 0
    for iIter in range( maxIter ):
        print( 'Pipek-Mezey localization iteration %i ...' % (iIter+1) )

        # print( ix_pair.shape )
        iPair = 0
        for (s,t) in ix_pair: # a pair of MOs: s and t
            A = 0
            B = 0
            for ix in aoIx:
                # Integrals < s| P_iat | t > for atom iat:
                # For nonorthogonal basis, we have Pipek-Mezey formula in
                # matrix form. 
                # Note that vectorization can speed up significantly.
                # Pst = ( C_t[rows in A]' @ S[rows in A] @ C_s[:] + 
                #         C_s[rows in A]' @ S[rows in A] @ C_t[:] ) / 2
                # Pss = C_s[rows in A]' @ S[rows in A] @ C_s[:]
                Pss = C_LMO[ix-1,s].T @ S[ix-1,:] @ C_LMO[:,s]
                Ptt = C_LMO[ix-1,t].T @ S[ix-1,:] @ C_LMO[:,t]
                Pst = ( C_LMO[ix-1,t].T @ S[ix-1,:] @ C_LMO[:,s] + \
                        C_LMO[ix-1,s].T @ S[ix-1,:] @ C_LMO[:,t] ) / 2
                dPst = Pss - Ptt
                A += Pst * Pst  -  dPst * dPst /4 
                B += Pst * dPst

            # Determine the rotation angle:
            # Trivial case: A = B = 0, no rotation required
            if A == 0 and B == 0:
                continue
            # General case for maximization: (gamma == alpha => a)
            R = sqrt( A*A + B*B )
            a = atan2( B/R, -A/R ) / 4
            sina = sin(a)
            cosa = cos(a)
            # Remove already converged MO pairs:
            if abs(sina) < lowTOL*100:
                ix_pair = np.delete( ix_pair, iPair, 0 )
                continue
            Q = np.array([ [ cosa, -sina], [sina, cosa] ])
            # Do the rotation and update coefficient matrix C_LMO:
            C_LMO[:,(s,t)] = C_LMO[:,(s,t)] @ Q

            iPair += 1
        # endfor (s,t)

        # Calculate the objective functional:
        f = 0
        for i in range( NE ):
            for ix in aoIx:
                # Atomic population of each atom
                # The following vectorization is preferred:
                qA = C_LMO[ix-1,i].T @ S[ix-1,:] @ C_LMO[:,i]
                f += qA * qA

        # Check convergence:
        dev = ( f - f0 ) / NE
        if abs(dev) < TOL:
            print( 'Pipek-Mezey localization converged with an accuracy of '
                   '%E after %i iterations' % ( abs(dev), iIter+1 ) )
            break
        elif dev < 0:
            raise RuntimeError( 'Failure in Pipek-Mezey localization: '
                    'Unexpected decrease of the objective functional: '
                    '%E < %E' %  ( f, f0 ) )
        elif iIter == maxIter-1:
            raise RuntimeError( 'Failure in Pipek-Mezey localization: '
                    'maximum number of iteractions exceeded (%i) and the '
                    'convergence criterion is not met (%E > %E)' % 
                    (maxIter, abs(dev), TOL) )
        # Update f0:
        f0 = f
        print( '    Current objective function = %.12f' % f )

        # Refinement when the low accuracy is fulfilled:
        if abs(dev) < lowTOL:
            ix_pair = np.copy( ix_pair0 )
            lowTOL = TOL
    # endfor iIter


    # CMO --> LMO transformation matrix:
    Q = lstsq( C[:,0:NE], C_LMO )[0]
    # Check if the left inverse is correctly computed by lstsq():
    assert( norm( C[:,0:NE] @ Q - C_LMO ) / NE < 1E-14 )
    # Fock matrix in LMO basis:
    # NOTE: The Fock matrix of CMOs in AO basis is diagonal, which is equal to 
    #       np.diag( E[0:NE] ). However, in general cases, where either the CMO 
    #       coefficients are expressed in NAO/Lowdin basis, or the input MOs
    #       are not CMOs, the Fock matrix is NOT diagonal.
    F_MO = C[:,0:NE].T @ F @ C[:,0:NE]  # AO --> MO
    F_LMO = Q.T @ F_MO @ Q
    # Get LMO energies:
    E_LMO = np.diag( F_LMO )

    # Sort LMOs by energies (from lower to higher):
    ix = np.argsort( E_LMO )
    C_LMO = C_LMO[:,ix]
    E_LMO = E_LMO[ix]

    return ( C_LMO, E_LMO )
# enddef PipekMezeyLMO_1spin_nonortho_fast()


# Extended version of PipekMezeyLMO_1spin() for MOs in nonorthogonal basis:
#   C: MO coefficient matrix in AO basis
#   S: Overlap matrix between AO basis
#   F: Fock matrix in AO basis
#   aoIx: AO indices for each atom
#   NE: Number of electrons
#   maxInter: Maximum iterations to be done
def PipekMezeyLMO_1spin_nonortho( C, S, F, aoIx, NE, maxIter=9000 ):
    # Tolerance for convergence:
    TOL = 1E-10

    NBas, NMO = C.shape
    if NE > NMO:
        raise ValueError( 'Number of electrons (%i) must not be greater than '
                'the number of molecular orbitals (%i)' % (NE, NMO) )

    # Initiate variables:
    C_LMO = np.copy( C[:,0:NE] )

    # Array of contributions of MOs to the objective functional:
    f_arr = np.zeros( NE )
    for i in range( NE ):
        for ix in aoIx:
            # Atomic population of each atom
            # The following vectorization is preferred:
            qA = C_LMO[ix-1,i].T @ S[ix-1,:] @ C_LMO[:,i]
            f_arr[i] += qA * qA
    # Sort by f_arr[], in descending order (i.e. from more localized to less)
    ix_sort_f = np.argsort( -f_arr )

    # Objective functional, D^{-1}, to be maximized:
    f0 = 0
    for iIter in range( maxIter ):
        print( 'Pipek-Mezey localization iteration %i ...' % (iIter+1) )

        for ix_s in np.arange( NE ): # MO s
            s = ix_sort_f[ ix_s ]
            for ix_t in np.arange( ix_s+1, NE ): # MO t
                t = ix_sort_f[ ix_t ]

                A = 0
                B = 0
                for ix in aoIx:
                    # Integrals < s| P_iat | t > for atom iat:
                    # For nonorthogonal basis, we have Pipek-Mezey formula in
                    # matrix form. 
                    # Note that vectorization can speed up significantly.
                    # Pst = ( C_t[rows in A]' @ S[rows in A] @ C_s[:] + 
                    #         C_s[rows in A]' @ S[rows in A] @ C_t[:] ) / 2
                    # Pss = C_s[rows in A]' @ S[rows in A] @ C_s[:]
                    Pss = C_LMO[ix-1,s].T @ S[ix-1,:] @ C_LMO[:,s]
                    Ptt = C_LMO[ix-1,t].T @ S[ix-1,:] @ C_LMO[:,t]
                    Pst = ( C_LMO[ix-1,t].T @ S[ix-1,:] @ C_LMO[:,s] + \
                            C_LMO[ix-1,s].T @ S[ix-1,:] @ C_LMO[:,t] ) / 2
                    dPst = Pss - Ptt
                    A += Pst * Pst  -  dPst * dPst /4 
                    B += Pst * dPst
 
                # Determine the rotation angle:
                # Trivial case: A = B = 0, no rotation required
                if A == 0 and B == 0:
                    continue
                # General case for maximization: (gamma == alpha => a)
                R = sqrt( A*A + B*B )
                a = atan2( B/R, -A/R ) / 4
                sina = sin(a)
                cosa = cos(a)
                Q = np.array([ [ cosa, -sina], [sina, cosa] ])
                # Do the rotation and update coefficient matrix C_LMO:
                C_LMO[:,(s,t)] = C_LMO[:,(s,t)] @ Q
            # endfor t
        # endfor s

        # Calculate the objective functional:
        f = 0
        for i in range( NE ):
            for ix in aoIx:
                # Atomic population of each atom
                # The following vectorization is preferred:
                qA = C_LMO[ix-1,i].T @ S[ix-1,:] @ C_LMO[:,i]
                f += qA * qA

        # Check convergence:
        dev = ( f - f0 ) / NE
        if abs(dev) < TOL:
            print( 'Pipek-Mezey localization converged with an accuracy of '
                   '%E after %i iterations' % ( abs(dev), iIter+1 ) )
            break
        elif dev < 0:
            raise RuntimeError( 'Failure in Pipek-Mezey localization: '
                    'Unexpected decrease of the objective functional: '
                    '%E < %E' %  ( f, f0 ) )
        elif iIter == maxIter-1:
            raise RuntimeError( 'Failure in Pipek-Mezey localization: '
                    'maximum number of iteractions exceeded (%i) and the '
                    'convergence criterion is not met (%E > %E)' % 
                    (maxIter, abs(dev), TOL) )
        # Update f0:
        f0 = f
        print( '    Current objective function = %.12f' % f )
    # endfor iIter


    # CMO --> LMO transformation matrix:
    Q = lstsq( C[:,0:NE], C_LMO )[0]
    # Check if the left inverse is correctly computed by lstsq():
    assert( norm( C[:,0:NE] @ Q - C_LMO ) / NE < 1E-14 )
    # Fock matrix in LMO basis:
    # NOTE: The Fock matrix of CMOs in AO basis is diagonal, which is equal to 
    #       np.diag( E[0:NE] ). However, in general cases, where either the CMO 
    #       coefficients are expressed in NAO/Lowdin basis, or the input MOs
    #       are not CMOs, the Fock matrix is NOT diagonal.
    F_MO = C[:,0:NE].T @ F @ C[:,0:NE]  # AO --> MO
    F_LMO = Q.T @ F_MO @ Q
    # Get LMO energies:
    E_LMO = np.diag( F_LMO )

    # Sort LMOs by energies (from lower to higher):
    ix = np.argsort( E_LMO )
    C_LMO = C_LMO[:,ix]
    E_LMO = E_LMO[ix]

    return ( C_LMO, E_LMO )
# enddef PipekMezeyLMO_1spin()


# Pipek-Mezey localization of molecular orbitals
#   C: MO coefficient matrix in orthonormal AO basis (e.g., NAOs, Lowdin AOs)
#      C is a single matrix in spin-restricted case
#      C is a tuple, C = (Ca,Cb) in spin-unrestricted case
#   F: Fock matrix in orthonormal AO basis (e.g., NAOs, Lowdin AOs)
#      F is a single matrix in spin-restricted case
#      F is a tuple, F = (Fa,Fb) in spin-unrestricted case
#   aoIx: AO indices for each atom
#   NE: Number of electrons
#      NE is a single number in spin-restricted case
#      NE is a tuple, NE = (NEa,NEb) in spin-unrestricted case
def PipekMezeyLMO( C, F, aoIx, NE, algo='normal', maxIter=9000 ):
    # Algorithm:
    if algo == 'fast':
        ifFast = True
    elif algo == 'normal':
        ifFast = False
    else:
        raise ValueError( 'Unrecognized value for argument algo '
                          'in function PipekMezeyLMO()' )
    # Determine whether this is spin-unrestricted or not:
    if type(C) is tuple:
        spin = 1 # spin-unrestricted
        print( 'Pipek-Mezey localization for spin-unrestricted calculation' )
        Ca, Cb = C
        NBas, NMO = Ca.shape
        Fa, Fb = F
        NEa, NEb = NE
    else:
        print( 'Pipek-Mezey localization for spin-restricted calculation' )
        spin = 0 # spin-restricted
        Ca = np.copy( C ) # pass by value
        NBas, NMO = Ca.shape
        Fa = np.copy( F ) # pass by value
        NEa = NE // 2

    ## Check linear dependence in the basis:
    #if NBas != NMO:
    #    raise ValueError( 'Linear dependence in the basis is found. '
    #            'The overlap matrix must be explicitly calculated.' )

    # For alpha electrons:
    if ifFast:
        C_LMO_a, E_LMO_a = PipekMezeyLMO_1spin_fast( Ca, Fa, aoIx, NEa, \
                                                     maxIter )
    else:
        C_LMO_a, E_LMO_a = PipekMezeyLMO_1spin( Ca, Fa, aoIx, NEa, maxIter )
    Ca = np.hstack( (C_LMO_a, Ca[:,NEa:NMO]) ) # Merge occ. & virt. MOs
    Ea_virt = np.diag( Ca[:,NEa:NMO].T @ Fa @ Ca[:,NEa:NMO] )
    Ea = np.hstack( (E_LMO_a, Ea_virt ) ) # Merge occ. & virt. MOs
    # For spin-unrestricted, run for beta as well:
    if spin == 1:
        if ifFast:
            C_LMO_b, E_LMO_b = PipekMezeyLMO_1spin_fast( Cb, Fb, aoIx, NEb, \
                                                         maxIter )
        else:
            C_LMO_b, E_LMO_b = PipekMezeyLMO_1spin( Cb, Fb, aoIx, NEb, \
                                                         maxIter )
        Cb = np.hstack( (C_LMO_b, Cb[:,NEb:NMO]) ) # Merge occ. & virt. MOs
        Eb_virt = np.diag( Cb[:,NEb:NMO].T @ Fb @ Cb[:,NEb:NMO] )
        Eb = np.hstack( (E_LMO_b, Eb_virt ) ) # Merge occ. & virt. MOs

    # Return:
    if spin == 0:
        return ( Ca, Ea )
    else:
        return ( (Ca, Cb), (Ea, Eb) )
# enddef PipekMezeyLMO()


# Extended version of PipekMezeyLMO() for MOs in nonorthogonal basis:
# Pipek-Mezey localization of molecular orbitals
#   C: MO coefficient matrix in AO basis
#      C is a single matrix in spin-restricted case
#      C is a tuple, C = (Ca,Cb) in spin-unrestricted case
#   S: Overlap matrix between AO basis
#   F: Fock matrix in AO basis
#      F is a single matrix in spin-restricted case
#      F is a tuple, F = (Fa,Fb) in spin-unrestricted case
#   aoIx: AO indices for each atom
#   NE: Number of electrons
def PipekMezeyLMO_nonortho( C, S, F, aoIx, NE, algo='normal', maxIter=9000 ):
    # Algorithm:
    if algo == 'fast':
        ifFast = True
    elif algo == 'normal':
        ifFast = False
    else:
        raise ValueError( 'Unrecognized value for argument algo '
                          'in function PipekMezeyLMO_nonortho()' )
    # Determine whether this is spin-unrestricted or not:
    if type(C) is tuple:
        spin = 1 # spin-unrestricted
        Ca, Cb = C
        NBas, NMO = Ca.shape
        Fa, Fb = F
        NEa, NEb = NE
    else:
        spin = 0 # spin-restricted
        Ca = np.copy( C ) # pass by value
        NBas, NMO = Ca.shape
        Fa = np.copy( F ) # pass by value
        NEa = NE // 2

    # Check linear dependence in the basis:
    if NBas != NMO:
        raise ValueError( 'Linear dependence in the basis is found. '
                'The overlap matrix must be explicitly calculated.' )

    # For alpha electrons:
    if ifFast:
        C_LMO_a, E_LMO_a = PipekMezeyLMO_1spin_nonortho_fast( Ca, S, Fa, \
                                                          aoIx, NEa, maxIter )
    else:
        C_LMO_a, E_LMO_a = PipekMezeyLMO_1spin_nonortho( Ca, S, Fa, \
                                                          aoIx, NEa, maxIter )
    Ca = np.hstack( (C_LMO_a, Ca[:,NEa:NMO]) ) # Merge occ. & virt. MOs
    Ea_virt = np.diag( Ca[:,NEa:NMO].T @ Fa @ Ca[:,NEa:NMO] )
    Ea = np.hstack( (E_LMO_a, Ea_virt ) ) # Merge occ. & virt. MOs
    # For spin-unrestricted, run for beta as well:
    if spin == 1:
        if ifFast:
            C_LMO_b, E_LMO_b = PipekMezeyLMO_1spin_nonortho_fast( Cb, S, Fb, \
                                                          aoIx, NEb, maxIter )
        else:
            C_LMO_b, E_LMO_b = PipekMezeyLMO_1spin_nonortho( Cb, S, Fb, \
                                                          aoIx, NEb, maxIter )
        Cb = np.hstack( (C_LMO_b, Cb[:,NEb:NMO]) ) # Merge occ. & virt. MOs
        Eb_virt = np.diag( Cb[:,NEb:NMO].T @ Fb @ Cb[:,NEb:NMO] )
        Eb = np.hstack( (E_LMO_b, Eb_virt ) ) # Merge occ. & virt. MOs

    # Return:
    if spin == 0:
        return ( Ca, Ea )
    else:
        return ( (Ca, Cb), (Ea, Eb) )
# enddef PipekMezeyLMO_nonortho()

# Generate Pipek-Mezey LMOs in NAO basis:
def genLMOs( FchkInfo, naoInfo, CAONAO, inpMainName, algo='normal', \
             ifWriteLMOs=True ):
    spin = FchkInfo.spin
    # If LMOs have been obtained previously, load the results from saved files:
    if spin == 0: # spin-restricted
        file_lmo = inpMainName + '_CNAOLMO.dat'
        file_elmo = inpMainName + '_ELMO.dat'
        if os.path.exists( file_lmo ) and os.path.exists ( file_elmo ):
            CNAOLMO = np.load( file_lmo, allow_pickle=True )
            ELMO = np.load( file_elmo, allow_pickle=True )
            return ( CNAOLMO, ELMO )
    else:
        file_lmoa = inpMainName + '_CNAOLMOa.dat'
        file_lmob = inpMainName + '_CNAOLMOb.dat'
        file_elmoa = inpMainName + '_ELMOa.dat'
        file_elmob = inpMainName + '_ELMOb.dat'
        if os.path.exists( file_lmoa ) and os.path.exists ( file_elmoa ) \
           and os.path.exists( file_lmob ) and os.path.exists ( file_elmob ):
            CNAOLMOa = np.load( file_lmoa, allow_pickle=True )
            CNAOLMOb = np.load( file_lmob, allow_pickle=True )
            ELMOa = np.load( file_elmoa, allow_pickle=True )
            ELMOb = np.load( file_elmob, allow_pickle=True )
            return (CNAOLMOa, CNAOLMOb), (ELMOa, ELMOb) 
    
    # LMO in NAO basis + LMO energies:
    if spin == 0: # spin-restricted
        CNAO = np.linalg.pinv( CAONAO ) @ FchkInfo.C
        # Fock matrix in AO bais:
        F = FockMatrix( FchkInfo.C, FchkInfo.E )
        # Fock matrix in NAO bais:
        FNAO = CAONAO.T @ F @ CAONAO
        CNAOLMO, ELMO = PipekMezeyLMO( CNAO, FNAO, naoInfo.naoIx, \
                                       naoInfo.NE, algo )
        # Save the results to external files:
        CNAOLMO.dump( file_lmo )
        ELMO.dump( file_elmo )
        # For writing LMOs to a fchk file:
        CLMO = CAONAO @ CNAOLMO # LMO in AO basis
    else: # spin-unrestricted
        CNAOa = np.linalg.pinv( CAONAO ) @ FchkInfo.C[0]
        CNAOb = np.linalg.pinv( CAONAO ) @ FchkInfo.C[1]
        # Fock matrices in AO bais:
        Fa, Fb = FockMatrix( FchkInfo.C, FchkInfo.E )
        # Fock matriices in NAO bais:
        FNAOa = CAONAO.T @ Fa @ CAONAO
        FNAOb = CAONAO.T @ Fb @ CAONAO
        CNAOLMO, ELMO = PipekMezeyLMO( (CNAOa,CNAOb), (FNAOa,FNAOb), \
                              naoInfo.naoIx, (naoInfo.NEa,naoInfo.NEb), algo )
        # Save the results to external files:
        CNAOLMOa, CNAOLMOb = CNAOLMO
        ELMOa, ELMOb = ELMO
        CNAOLMOa.dump( file_lmoa )
        CNAOLMOb.dump( file_lmob )
        ELMOa.dump( file_elmoa )
        ELMOb.dump( file_elmob )
        # For writing LMOs to a fchk file:
        CLMOa = CAONAO @ CNAOLMOa # alpha LMO in AO basis
        CLMOb = CAONAO @ CNAOLMOb # beta LMO in AO basis
        CLMO = ( CLMOa, CLMOb )

    # Write the LMOs in a fchk file:
    if ifWriteLMOs:
        outp = inpMainName + '_LMO.fchk'
        writeFchkOrb( inpMainName +'.fchk', outp, CLMO, ELMO )
        print( 'File %s written' % outp  )

    return ( CNAOLMO, ELMO )
# enddef genLMOs()


# Generate CMOs in NAO basis:
def genCMOs( FchkInfo, CAONAO ):
    # CMO in NAO basis + CMO energies:
    CNAO = np.linalg.pinv( CAONAO ) @ FchkInfo.C

    return ( CNAO, FchkInfo.E )
# enddef genCMOs()


# Density matrix built from localized molecular orbitals (LMOs):
# Note: the basis should be orthonormal set, such as NAOs, Lowdin AOs.
def densityMatrixLMO( C, spin=0 ):
    if type(C) is tuple: # spin-unrestricted
        Ca, Cb = C
        NEa = Ca.shape[1]
        NEb = Cb.shape[1]
        NE = ( NEa, NEb )
    else: # spin-restricted
        NE = C.shape[1] * 2
    return densityMatrix( C, NE )
# enddef densityMatrixLMO()

# Total energy as sum of energies of given localized molecular orbitals (LMOs):
#   ELMO:  Energies of LMOs
#   lmoIx: Indices of the given set of LMOs (e.g., pi-LMOs, or LMOs of interest)
#          NOTE: lmoIx starts from 0 (not 1)
#   spin:  =0 --> spin-unrestricted;  =1 --> spin-restricted
def LMOEnergy( ELMO, lmoIx, spin=0 ):
    # Doubly occupied LMOs (for closed-shell system):
    if spin == 0:
        return ELMO[ lmoIx ].sum() * 2
    else:
        return ELMO[ lmoIx ].sum()
# enddef LMOEnergy()

# Make alpha and beta LMOs match each other
def matchAlphaBetaLMOs( C_LMOa, C_LMOb ):
    U, _, V = svd( C_LMOa.T @ C_LMOb )
    # Set ones and zeros in  U and V:
    U0 = np.zeros( U.shape )
    i = 0
    for U_col in U.T:
        maxix = abs( U_col ).argmax()
        U0[ maxix, i ] = np.sign( U_col[maxix] )
        i += 1
    # Assert that U0 is unitary:
    assert( np.array_equal( U0.T @ U0, np.eye( U0.shape[0], U0.shape[1] )  ) )
    #print( U )
    print( U0 )
    V0 = np.zeros( V.shape )
    i = 0
    for V_col in V.T:
        maxix = abs( V_col ).argmax()
        V0[ maxix, i ] = np.sign( V_col[maxix] )
        i += 1
    # Assert that V0 is unitary:
    assert( np.array_equal( V0.T @ V0, np.eye( V0.shape[0], V0.shape[1] )  ) )
    #print( V )
    print( V0 )

    return C_LMOa @ U0, C_LMOb @ V0

# enddef matchAlphaBetaLMOs()
