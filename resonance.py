# Resonance theory based on decomposition of density matrix
#
#  Created on Mar 5, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    Jul 23, 2020:
#    - Add function rao_uni_manual() allowing for setting RAO phases manually
#
#    Jun 13, 2020:
#    - Found a remaining bug (still!) in function rao_uni(), which 
#    wrongly determines the phase match between two RAOs, in cases like ozone
#    molecule, where the two terminal O atoms have positive pi bond orders but
#    the associated interatomic block of DM is dominated by negative values
#    (This is because in the pi HOMO the atomic orbitals of the two atoms are
#    off-phase)
#    This issue is yet to be solved. Note that O3 is really a very special
#    case, which merits a detailed future study!!!
#
#    Jun 3, 2020:
#    - Found and fixed a minor bug remaining in function rao_uni(), which 
#    wrongly determines the phase match between two RAOs, due to the reversed 
#    ordering of the RAO-iat and RAO-iat1 in denisty matrix block D12
#
#    May 31, 2020:
#    - Fixed the still remaining bug in function rao_uni() for determination
#    of RAO phases, by applying a new, simpler algorithm
#    In more detail, the in-phase or off-phase characteristics between two
#    RAOs is just the 1st element of the transformed off-diagonal block of the
#    diatomic block of DM, which is associated with the two RAOs and can be
#    directly obtain as U1.T @ D12_off @ U2 (where D12_off is the off-diagonal
#    block of the diatomic block of DM; U1 and U2 are the transform vector of
#    RAO-1 and RAO-2, respectively
#
#    - Fixed the remaining bug in function rao_uni(), for determination of RAO
#    phases was erroneous for non-planar molecules
#    This is due to the fact that the actual DM off-diagonal element is not
#    always positive for an LMO, e.g., the pi-LMO of Al4(2+)
#
#    Apr 30, 2020:
#    - Fixed the bug in function rao_uni(), that is the determination of RAO
#    phases was erroneous for non-planar molecules
#
#    From Apr 28, 2020:
#    - Added implementation using wave functions obtained from simple Hueckel  
#      Molecular Orbital theory
#
#    Apr 24, 2020:
#      - Allow for straightforward specification of Lewis structures parsed
#        from strings, in function dmrt_spec()
#
#    From Mar 27, 2020:
#      - Added functions rao_uni_general() allowing RAOs to be represented by 
#        multi-AO-basis functions for each Lewis strcuture
#

import numpy as cp
import copy as cp
from transform import *
from comb import *
from util import *
from writeFchkOrb import *


# Gross atomic population and Wiberg bond orders, obtained from density matrix:
# Note: the basis should be orthonormal set, such as NAOs, Lowdin AOs.
#   D:     Density matrix 
#   aoIx:  AO (forming the orthonormal basis) indices for each atom
#   q:     Gross atomic population
#   b:     BO-Valence matrix, corresponding to off-diagonal and diagonal entries
def popDM( D, aoIx ):
    NAt = len( aoIx )
    q = np.zeros( (NAt,1) )
    b = np.zeros( (NAt,NAt) )
    for i in range( NAt ):
        ix = aoIx[i]-1 # AO indices for atom i
        q[i] = np.trace( D[np.ix_( ix,ix )] )
        # Wiberg bond orders:
        for j in range( i+1, NAt ):
            jx = aoIx[j]-1 # AO indices for atom j
            b[i,j] = np.trace( D[np.ix_( ix,jx )] @ D[np.ix_( jx,ix )] )
            b[j,i] = b[i,j]

    # Atomic valences as diagonal entries of b[]:
    for i in range( NAt ):
        b[i,i] = b[i,:].sum()

    return ( q, b )
# enddef popLMO


# Resonance molecular orbitals (RMOs) corresponding to a given Lewis strcuture:
#   D: Density matrix of the whole system
#   aoIx: AO (forming the orthonormal basis) indices for each atom
#   LP: Array of atomic indices of lone pairs
#   BD: Array of atomic indices of bonds
def rmo( D, aoIx, LP, BD ):
    nLP = len( LP )
    nBD = BD.shape[0]
    NBas = len( D )

    # RMOs associated with LPs:
    rmoLP = np.zeros( ( NBas, nLP ) )
    k = 0
    #for k in range( nLP ):
    for iat in LP-1:
        ix = aoIx[iat]-1 # AO indices for atom iat
        # Diagonalization of subblock of each atom to get LPs:
        U, _ = maxEig( D[np.ix_( ix, ix )] )
        rmoLP[ ix, k ] = U
        # print( rmoLP[:,k] )
        k += 1

    # RMOs associated with BDs:
    rmoBD = np.zeros( ( NBas, nBD ) )
    k = 0
    for iat1, iat2 in BD-1:
        # AO indices for atoms iat1 and iat2
        ix = concatenate( [ aoIx[iat1]-1, aoIx[iat2]-1 ] )
        # Diagonalization of subblock of each atom pair to get BDs:
        U, _ = maxEig( D[np.ix_( ix, ix )] )
        rmoBD[ ix, k ] = U
        # print( rmoBD[:,k] )
        k += 1

    # Coefficient matrix of RMOs in AO basis:
    return hstack( ( rmoLP, rmoBD ) )
# enddef rmo()


# Resonance atomic orbitals (RAOs) corresponding to a given Lewis strcuture:
# NOTE: This set of latter RAOs depends on each individual Lewis strcuture
#
#   D: Density matrix of the whole system in the AO basis (NAO, Lowdin's, etc.)
#   aoIx: AO (forming the orthonormal basis) indices for each atom
#   LP: Array of atomic indices of lone pairs
#   BD: Array of atomic indices of bonds
def rao( D, aoIx, LP, BD ):
    nLP = len( LP )
    if len( BD ) == 0:
        nBD = 0
    else:
        nBD = BD.shape[0]
    NBas = len( D )
    Nc = nLP+nBD*2  # Number of centers

    # Initialize RAO coefficient matrix C, for the Nc centers:
    C = np.zeros( ( NBas, Nc ) )
    # Initialize density matrix in RAO basis:
    Drao = np.eye( Nc )

    # RAOs associated with LPs:
    k = 0 # To record index of RAO
    for iat in LP-1:
        ix = aoIx[iat]-1 # AO indices for atom iat
        # Diagonalization of subblock of each atom to get LPs:
        C[np.ix_( ix, [k] )], _ = maxEig( D[np.ix_( ix, ix )] )
        #print( C[:,k] )
        # Density matrix in RAO basis (for closed-shell systems):
        Drao[k,k] = 2 # Exact double occupancy (a lone pair)
        k += 1

    # RAOs associated with BDs:
    ic = nLP # To record the appearance index of centers
    for ib in range( nBD ):
        iat1, iat2 = BD[ib,:]-1
        # AO indices for atoms iat1 and iat2
        ix1 = aoIx[iat1]-1
        ix2 = aoIx[iat2]-1
        # Interatomic blocks between iat1 and iat2:
        D12 = D[np.ix_( ix1, ix2 )]
        D21 = D[np.ix_( ix2, ix1 )]
        # Diagonalization of bond order matrices to get RAOs for iat1 and iat2:
        C1, E1 = maxEig( D12 @ D21 )
        C2, E2 = maxEig( D21 @ D12 )
        # Assure that the two RAOs are in phase:
        D12_new = C1.T @ D12 @ C2
        #D21_new = C2.T @ D21 @ C1  # which should be equal to D12_new
        #print( BD[ib,:], D12_new )
        # If off phase, then reverse RAO-2:
        if D12_new < 0:
            C2 = -C2

        # Update C[] vectors:
        C[np.ix_( ix1, [k] )] = C1
        C[np.ix_( ix2, [k+1] )] = C2
        #writeMat( C[:,0], '/dev/stdout', ' %8.4f' )
        #input()
        k += 2

        # Density matrix in RAO basis (for closed-shell systems):
        Drao[ic,ic+1] = 1 # Exact covalent bond (equal sharing of electrons)
        Drao[ic+1,ic] = 1 # Exact covalent bond (equal sharing of electrons)
        ic += 2


    # Obtain density matrix in AO basis (NAOs, Lowdin AOs, etc.):
#    writeMat( C, '/dev/stdout', ' %8.4f' )
    D_RAO = C @ Drao @ C.T
    #print( np.trace(D_RAO) )

    return ( C, D_RAO )
# enddef rao()


# Universal resonance atomic orbitals (RAOs) for a given molecule:
# NOTE: Mind the difference between rao_uni() and rao(), that is, the former
# RAOs are the same for all Lewis strcutures whereas the latter RAOs depend on
# each individual Lewis strcuture
#   D: Density matrix of the whole system in the AO basis (NAO, Lowdin's, etc.)
#   aoIx: AO (forming the orthonormal basis) indices for each atom
#   flipAtoms: Indices of atoms whose RAO phase is flipped
#              -  Default: None --> The phases are determined automatically
#              -  if flipAtoms == [ 0 ], then no flip is done, nor automatic
#                 detection is performed
def rao_uni( D, aoIx, flipAtoms=None ):

    if flipAtoms is not None:
        return rao_uni_manual( D, aoIx, flipAtoms )

    Nc = len( aoIx ) # Number of centers = number of atoms
    NBas = len( D )

    # Bond order matrices:
    BO = np.zeros( ( Nc, Nc ) )
    for i in range( Nc ):
        ix = aoIx[i]-1 # AO indices for atom i
        # Negative sign represents atomic population:
        BO[i,i] = -( D[np.ix_( ix, ix )] @ D[np.ix_( ix, ix )] ).trace()
        for j in range( i+1, Nc ):
            jx = aoIx[j]-1 # AO indices for atom j
            BO[i,j] = ( D[np.ix_( ix, jx )] @ D[np.ix_( jx, ix )] ).trace()
            BO[j,i] = BO[i,j]
    #writeMat( BO, '/dev/stdout', ' %12.6f' )
    # Build up connectivity chain via the nearest possbile neighbors:
    # A chain connecting all atoms through the nearest distances between them
    chain =  [ 0 ] # Atom-1 is the first atom on the chain
    pool = list( range( 1, Nc ) ) # Atoms not been picked yet
    while True:
        # Determine the nearest atom pair between pool[] and chain[]:
        maxBO = 0.
        nearestIx = -1
        for i in pool:
            for j in chain:
                if BO[i,j] > maxBO:
                    maxBO = BO[i,j]
                    nearestIx = i
        # Update chain[] and pool[]:
        chain.append( nearestIx )
        pool.remove( nearestIx )
        # Completed:
        if len( pool ) == 0:
            break
    assert( len(chain) == Nc ) # All atoms must be in the chain
    #print( np.array(chain)+1 )
    #input()
    # Get the list of reference atoms for determining RAO phases:
    refList = [ -1 ] # List of reference atoms; 1st atom has no ref
    for i in range( 1, Nc ):
        ix = chain.index( i )
        BOi = BO[ i, chain[0:ix] ]
        refList.append( chain[  np.nonzero( BOi == max(BOi) )[0][0]  ] )
    #print( np.arange(Nc)+1 )
    #print( np.array(refList)+1 )
    #input()


    # Initialize RAO coefficient matrix C, for the Nc centers:
    C = np.zeros( ( NBas, Nc ) )

    # RAOs associated with each atom:
    Uref = [] # RAO vectors for reference atoms
    for iat in chain:
        ix = aoIx[iat]-1 # AO indices for atom iat
        # Diagonalization of subblock of each atom to get the associated RAO:
        U, _ = maxEig( D[np.ix_( ix, ix )] )

        # For the first atom:
        if iat == 0:
            # Reverse the phase of RAO-1 to make maximum component positive
            if U[ np.argmax(abs(U)) ] < 0:
                U = -U
                print( 'RAO #%-3i reversed' %( iat+1 ) )

        # For all other atoms, assure that the two RAOs are in phase:
        else:
            iat1 =  refList[iat]
            print( 'Determining phase of RAO #%i using reference RAO #%i ...'
                    % ( iat+1, iat1+1 ) )
            U1 = Uref[ chain.index(iat1) ]
            #writeMat( U, '/dev/stdout', ' %12.6f' )
            #print( '-'*60 )
            #writeMat( U1, '/dev/stdout', ' %12.6f' )
            len1 = len( U1 )
            len2 = len( ix )
            len12 = len1 + len2
            ix12 = np.r_[ aoIx[iat1]-1, ix ]
            D12 = D[np.ix_( ix12, ix12 )]
            #writeMat( D12, '/dev/stdout', ' %8.3f' )
            # Off-diagonal block of D12:
            D12_off = D12[ 0:len1, len1:len12 ] 
            # The transformed element of D12_off associated with the bonding
            # between RAO-iat1 and RAO-iat:
            phase = ( U1.T @ D12_off @ U )[0][0]
            #print( phase )
            # If off phase, then reverse the RAO of iat:
            if phase < 0:
                U = -U
                print( 'RAO #%-3i reversed' %( iat+1 ) )
                #writeMat( U, '/dev/stdout', ' %8.4f' )
                #input()

        C[np.ix_( ix, [iat] )] = U
        Uref.append( U )

        #print( C[:,iat] )
    return C
# enddef rao_uni()


# Manual mode of rao_uni():
#
# Universal resonance atomic orbitals (RAOs) for a given molecule:
# NOTE: Mind the difference between rao_uni() and rao(), that is, the former
# RAOs are the same for all Lewis strcutures whereas the latter RAOs depend on
# each individual Lewis strcuture
#   D: Density matrix of the whole system in the AO basis (NAO, Lowdin's, etc.)
#   aoIx: AO (forming the orthonormal basis) indices for each atom
#   flipAtoms: Indices of atoms whose RAO phase is flipped
#              if flipAtoms == [ 0 ], then no flip is done, nor automatic
#              detection is performed
def rao_uni_manual( D, aoIx, flipAtoms ):

    print( "Setting the RAO phases manually ... " )

    Nc = len( aoIx ) # Number of centers = number of atoms
    NBas = len( D )

    # Initialize RAO coefficient matrix C, for the Nc centers:
    C = np.zeros( ( NBas, Nc ) )

    # RAOs associated with each atom:
    for iat in range( Nc ):
        ix = aoIx[iat]-1 # AO indices for atom iat
        # Diagonalization of subblock of each atom to get the associated RAO:
        U, _ = maxEig( D[np.ix_( ix, ix )] )
       
        # Flip the phase manually:
        if flipAtoms[0] > 0:
            if (iat+1) in flipAtoms:
                U = -U
                print( 'RAO #%-3i reversed' %( iat+1 ) )
       
        C[np.ix_( ix, [iat] )] = U

    return C
# enddef rao_uni_manual()

# Density matrix associtaed with a given Lewis structure in AO basis 
# (NAOs, Lowdin AOs, etc.), based on universal resonance atomic orbitals (RAOs) 
def DM_RAO_uni( C_RAO, LP, BD ):
    nLP = len( LP )
    if len( BD ) == 0:
        nBD = 0
    else:
        nBD = BD.shape[0]

    # Initialize density matrix in RAO basis:
    Drao = np.eye( nLP+nBD*2 )

    for k in range( nLP ):
        # Density matrix in RAO basis (for closed-shell systems):
        Drao[k,k] = 2 # Exact double occupancy (a lone pair)

    atIx = list( LP-1 ) # A list to record index of atoms
    ic = nLP # A single variable to record the appearance index of centers
    for ib in range( nBD ):
        atIx.append( BD[ib,0]-1 )
        atIx.append( BD[ib,1]-1 )
        # Density matrix in RAO basis (for closed-shell systems):
        Drao[ic,ic+1] = 1 # Exact covalent bond (equal sharing of electrons)
        Drao[ic+1,ic] = 1 # Exact covalent bond (equal sharing of electrons)
        ic += 2
    #print( Drao )

    # Obtain density matrix in AO basis (NAOs, Lowdin AOs, etc.):
    C = C_RAO[ :, atIx ]
    #writeMat( C, '/dev/stdout', ' %8.4f' )
    D_RAO = C @ Drao @ C.T
    #print( np.trace(D_RAO) )

    return D_RAO
# enddef DM_RAO_uni()



# General and universal resonance atomic orbitals (RAOs) for a given molecule:
#   This is a general version of rao_uni(), allowing RAOs to be represented by 
#   multi-AO-basis functions for each Lewis strcuture.
#   Hopefully, this general set of RAO basis can COMPLETELY reproduce the
#   actual density matrix of the system in question
def rao_uni_general( D, aoIx ):
    TOL = 1E-10 # Tolerance to determine if an AO should be considered or not
    Nc = len( aoIx ) # Number of centers = number of atoms
    NBas = len( D )

    # Initialize RAO coefficient matrix C, for the Nc centers:
    C = np.zeros( ( NBas, 0 ) )
    iRAO = 0 # Index of RAOs
    raoIx = [] # Indices of RAOs for each atom

    # RAOs associated with each atom:
    for iat in range( Nc ):
        ix = aoIx[iat]-1 # AO indices for atom iat
        # Diagonalization of subblock of each atom to get the associated RAO:
        U, _ = nonzeroEig( D[np.ix_( ix, ix )], TOL )
        #U, _ = maxEig( D[np.ix_( ix, ix )] )
        nRAOi = U.shape[1] # Number of RAOs for atom iat
        print( '%2i RAOs found for atom #%i' % ( nRAOi, iat+1 ) )

        # Assure that all RAOs are in phase:
        if iat > 0:
            len2 = len(ix) # Number of basis function for atom iat
            len12 = len1 + len2
            Q = np.zeros( ( len12, 2 ) )
            Q[ len2:len12, [1] ] = U1
            for k in range( nRAOi ):
                Q[ 0:len2, [0] ] = U[:,[k]]
                Dtest = Q @ np.array([[1,1],[1,1]]) @ Q.T
                phase = 0
                for i in range( len2 ):
                    for j in range( len2, len12 ):
                        phase += Dtest[i,i] * Dtest[j,j] * Dtest[i,j]
                if phase < 0:
                    U[:,[k]] = -U[:,[k]]
                    print( '     RAO #%-2i of atom #%-3i reversed' 
                            % (k+1, iat+1 ) )
        else:
            # Reverse the phase of RAO-1 to make maximum component positive
            U1 = U[:,[0]] # 1st RAO of 1st atom
            if U1[ np.argmax(abs(U1)) ] < 0:
                U[:,[0]] = -U1 # NOTE: By this operation, U1 is also reversed
                print( '     RAO #1  of atom #1   reversed' )
            ix1 = np.copy( ix )
            U1 = np.copy( U[:,[0]] )
            len1 = len(ix1) # Number of basis function for 1st atom

            # Assure that rest of 1st atom's RAOs are in phase with 1st atom's
            # 1st RAO, U1:
            len12 = len1 * 2
            Q = np.zeros( ( len12, 2 ) )
            Q[ len1:len12, [1] ] = U1
            for k in range( 1, nRAOi ):
                Q[ 0:len1, [0] ] = U[:,[k]]
                Dtest = Q @ np.array([[1,1],[1,1]]) @ Q.T
                phase = 0
                for i in range( len1 ):
                    for j in range( len1, len12 ):
                        phase += Dtest[i,i] * Dtest[j,j] * Dtest[i,j]
                if phase < 0:
                    U[:,[k]] = -U[:,[k]]
                    print( '     RAO #%-2i of atom #%-3i reversed' 
                            % (k+1, iat+1 ) )

        # Append all RAOs of atom iat to C matrix:
        C = np.c_[ C, np.zeros( (NBas,nRAOi) ) ] # Add nRAOi columns to C matrix
        C[np.ix_( ix, np.arange(iRAO,iRAO+nRAOi) )] = U
        raoIx.append( np.arange( iRAO, iRAO+nRAOi ) )
        iRAO += nRAOi
        #writeMat( C, '/dev/stdout', ' %6.3f' )
        #print( raoIx )
        #input()


    return ( C, raoIx )
#enddef rao_uni_general()


# Density matrix associtaed with a given Lewis structure in AO basis 
# (NAOs, Lowdin AOs, etc.), based on general universal resonance atomic 
# orbitals (RAOs), allowing multi-RAO basis for each atom
#   This is a general version of function DM_RAO_uni()
def DM_RAO_uni_general( C_RAO, raoIx_all, LP, BD ):
    nLP = len( LP )
    if len( BD ) == 0:
        nBD = 0
    else:
        nBD = BD.shape[0]

    # Initialize density matrix in RAO basis:
    Drao = np.eye( nLP+nBD*2 )

    for k in range( nLP ):
        # Density matrix in RAO basis (for closed-shell systems):
        Drao[k,k] = 2 # Exact double occupancy (a lone pair)

    atIx = list( LP-1 ) # A list to record index of atoms
    ic = nLP # A single variable to record the appearance index of centers
    for ib in range( nBD ):
        atIx.append( BD[ib,0]-1 )
        atIx.append( BD[ib,1]-1 )
        # Density matrix in RAO basis (for closed-shell systems):
        Drao[ic,ic+1] = 1 # Exact covalent bond (equal sharing of electrons)
        Drao[ic+1,ic] = 1 # Exact covalent bond (equal sharing of electrons)
        ic += 2
    #print( Drao )

    # Number of RAOs for each atom:
    raoIx = []
    raoLen = [] # Number of RAOs for each atom
    for iat in atIx:
        raoIx.append( raoIx_all[iat] )
        raoLen.append( len(raoIx_all[iat]) )

    # Enumerate all combinations of RAOs of all atoms:
    #print( '==============', atIx )
    #print( '              ', raoIx, '      ', raoLen )
    NAt = len( atIx )
    ixRAO_list = [] # A list to record all combinations
    ixRAO = np.zeros( NAt, dtype='i' ) # Record current combination
    combRAO( NAt, raoIx, ixRAO, ixRAO_list )
    ixRAO_list = np.array( ixRAO_list ) # Convert to a numpy 2D array
    # Remove duplicated combinations:
    ixRAO_list = np.unique( ixRAO_list, axis=0 )
    #print( ixRAO_list )
    # Check if the total number of combinations is consistent:
    nComb = 1
    for n in raoLen:
        nComb *= n
    if nComb != ixRAO_list.shape[0]:
        raise AssertionError( 'Number of enumerated combinations by combRAO() '
                'is not correct: should be %i instead of %i' % 
                ( nComb, ixRAO_list.shape[0] ) )

    # Obtain density matrix in AO basis (NAOs, Lowdin AOs, etc.):
    D_RAO_list = []

    for ixRAO in ixRAO_list:
        C = C_RAO[ :, ixRAO ]
        #print( ixRAO )
        #writeMat( C, '/dev/stdout', ' %8.4f' )
        D_RAO = C @ Drao @ C.T
        #print( np.trace(D_RAO) )
        D_RAO_list.append( D_RAO )

    return D_RAO_list
# enddef DM_RAO_uni_general()


# Enumerate all combinations of RAOs of all atoms:
def combRAO( NAt, raoIx, ixRAO, ixRAO_list ):
    if NAt == 0:
        # print( 'ixRAO =', ixRAO )
        ixRAO_list.append( cp.deepcopy(ixRAO) )
        return

    for iat in range( NAt-1, -1, -1 ):
        for irao in raoIx[iat]:
            ixRAO[iat] = irao
            combRAO( NAt-1, raoIx, ixRAO, ixRAO_list )

    return
# enddef


# Fock matrix in RAO basis:
#   C_RAO: <NBas*Nc>   Coefficient matrix of RAOs in NAO basis
#   F_NAO: <NBas*NBas> Fock matrix in NAO basis
# Return:
#   <Nc*Nc> Fock matrix, whose diagonal entries give RAO energies
def FockRAO( C_RAO, F_NAO ):
    return C_RAO.T @ F_NAO @ C_RAO
# enddef FockRAO()


# Energy of a given Lewis structure:
#   F_RAO: <Nc*Nc> Fock matrix in RAO basis
#   LP:    <1*nLP> 1d array for the LPs of the Lewis structure
#   BD:    <1*nBD> 1d array for the BDs of the Lewis structure
def LewisEnergy( F_RAO, LP, BD, raoType='uni' ):
    nLP = len( LP )
    nBD = len( BD )
    E = 0

    if raoType == 'uni':
        # Lone pairs:
        for k in range( nLP ):
            iat = LP[k]-1 # Index of the atom for the LP
            # Doubly occupied RAO for the kth-LP (for closed-shell systems):
            E += F_RAO[iat,iat] * 2
        # Bonds:
        for ib in range( nBD ):
            iat1 = BD[ib,0]-1
            iat2 = BD[ib,1]-1
            E += F_RAO[iat1,iat1] + F_RAO[iat2,iat2] + F_RAO[iat1,iat2] * 2

    elif raoType == 'var':
        # Lone pairs:
        for k in range( nLP ):
            # Doubly occupied RAO for the kth-LP (for closed-shell systems):
            E += F_RAO[k,k] * 2
        # Bonds:
        ic = nLP # To record the appearance index of centers
        for ib in range( nBD ):
            ic2 = ic + 1
            E += F_RAO[ic,ic] + F_RAO[ic2,ic2] + F_RAO[ic,ic2] * 2
            ic += 2

    else:
        raise ValueError( 'Unrecognized value of argument raoType for '
                          'function LewisEnergy()' )

    return E
# enddef LewisEnergy()

# Inner (dot) product between two density matrices:
# This is a much more efficient realization, utilizing the fact that
# we need only the trace and both D1 and D2 are symmetric matrices
def innerDM( D1, D2 ):
    return ( D1*D2 ).sum()
# enddef innerDM()

# Inner (dot) product between two density matrices:
# This is a standard way for more general D1 and D2 matrices, but inefficient
def innerDM_standard( D1, D2 ):
    return np.trace( D1 @ D2 )
# enddef innerDM_standard()


# Outer product (overlap matrix) for an array of density matrices:
def outerDM( D_arr ):
    n = len( D_arr )
    S = np.zeros( (n,n) )

    for i in range(n):
        S[i,i] = innerDM( D_arr[i], D_arr[i] )
        for j in range( i+1, n ):
            S[i,j] = innerDM( D_arr[i], D_arr[j] )
            S[j,i] = S[i,j]

    return S
# enddef outerDM()

# Outer product (overlap matrix) for a set of Lewis resonace structures:
# NOTE: This is a specialized case of outerDM(), where all DM in D_arr[]
#       have the same trace, which is equal to 2*NE
#       Therefore, we use a more efficient way to calculate the overlap
def overlapDM_Lewis( D_arr ):
    NE2 = round( np.trace(D_arr[0]) )*2 # Two times of number of electrons
    n = len( D_arr )
    S = np.zeros( (n,n) )

    for i in range(n):
        S[i,i] = NE2
        for j in range( i+1, n ):
            S[i,j] = innerDM( D_arr[i], D_arr[j] )
            S[j,i] = S[i,j]

    # Determine the rank of S:
    _, E, _ = svd( S )
    rank = len( np.where( E > 1E-6 )[0] ) 

    return ( S, rank )
# enddef outerDM()


# Projection of a set of density matrices (basis) onto a given density matrix:
# D_arr:  List of component density matrices (i.e., the basis)
# D0: The density matrix to be projected onto
def projectDM( D_arr, D0 ):
    n = len( D_arr )
    P = np.zeros( (n,1) )

    for i in range(n):
        P[i] = innerDM( D_arr[i], D0 )

    return P
# enddef projectDM()


# Determine if a group of Lewis structures is essential (linearly independent 
# of all others) or not:
#  S: Overlap matrix between all considered Lewis structures
#  ix_g: list containg arrays of indices for degenerate Lewis structures
#  ix: <1d array> Indices of the given group of Lewis structures
def isEssentialLewis( S, ix_g, ix ):
    n_g = len( ix_g ) # Number of considered degenerate groups
    N = len( S ) # Number of all considered Lewis structures

    # Merge symmetrically equivalent columns of S:
    A0 = np.zeros( ( N, n_g ) )
    for j in range( n_g ):
        g = ix_g[j]
        for x in g:
            A0[:,j] += S[:,x]
    # Determine the rank:
    _, E, _ = svd( A0 )
    rank0 = len( np.where( E > 1E-6 )[0] ) 
    print( 'rank0 = ', rank0 )


    # Merge symmetrically equivalent columns of S, but excluding ix[]:
    A = np.zeros( ( N, n_g-1 ) )
    k = 0
    for j in range( n_g ):
        g = ix_g[j]
        if np.intersect1d(g,ix).size > 0:
            print( 'XXXXXX', g )
            continue
        for x in g:
            A[:,k] += S[:,x]
        k += 1
    # Determine the rank:
    _, E, _ = svd( A )
    rank = len( np.where( E > 1E-6 )[0] ) 
    print( 'size =', A.size, 'rank = ', rank )

    return
# enddef isEssentialLewis()



# Expansion of a given density matrix into sum of component density matrices:
# D0: The density matrix to be expanded
# D:  List of component density matrices (i.e., the basis)
# S:  Overlap (covariance) matrix between D[]
# P:  Projection vector of D[] into D0
# isNormalized: If impose the constraint that the sum of weights equals one
def expandDM( D0, S, P, isNormalized ):

    if isNormalized:
        N = len( S ) # Number of Lewis structures
        S1 = np.ones( ( N, N ) )
        P1 = np.ones( ( N, 1 ) )
        for i in range( 1, N ):
            S1[i,:] = S[i,:] - S[0,:]
            P1[i] = P[i] - P[0]

    #W = solve( S, P ) # Solved weights
    # Ordinary linear regression:
    #W = np.linalg.pinv(S@S) @ S @ P
    # Ridge linear regression:
    a = 1E-16
    if isNormalized:
        W = np.linalg.pinv( S1@S1 + a*np.eye(len(P1)) ) @ S1 @ P1
    else:
        W = np.linalg.pinv( S@S + a*np.eye(len(P)) ) @ S @ P
    
    R = P.T @ W / ( 2*np.trace(D0) )

    return ( W, R )
# enddef expandDM()


# Expansion of a given density matrix into sum of component density matrices
# with symmetry constraint, that is symmetrically equivalent Lewis structures
# must have the same weight
#
# D0:           The density matrix to be expanded
# D:            List of component density matrices (i.e., the basis)
# S:            Overlap (covariance) matrix between D[]
# P:            Projection vector of D[] into D0
# ix_g:         List containg arrays of indices for degenerate Lewis structures
# isNormalized: If impose the constraint that the sum of weights equals one
def expandDM_symm( D0, S, P, ix_g, isNormalized ):
    n_g = len( ix_g ) # Number of considered degenerate groups
    N = len( S ) # Number of Lewis structures

    if isNormalized:
        S1 = np.ones( ( N, N ) )
        P1 = np.ones( ( N, 1 ) )
        for i in range( 1, N ):
            S1[i,:] = S[i,:] - S[0,:]
            P1[i] = P[i] - P[0]

    # Merge symmetrically equivalent columns of S:
    A = np.zeros( ( N, n_g ) )
    for j in range( n_g ):
        ix = ix_g[j]
        for x in ix:
            if isNormalized:
                A[:,j] += S1[:,x]
            else:
                A[:,j] += S[:,x]

    # Ordinary linear regression:
    #Wsymm = np.linalg.pinv(A.T@A) @ A.T @ P
    # Ridge linear regression:
    a = 1E-6 #1E-10
    if isNormalized:
        Wsymm = np.linalg.pinv( A.T@A + a*np.eye(n_g) ) @ A.T @ P1
    else:
        Wsymm = np.linalg.pinv( A.T@A + a*np.eye(n_g) ) @ A.T @ P

    # Assign back Ws to W using symmetry:
    W = np.zeros( (N, 1) )
    for j in range( n_g ):
        ix = ix_g[j]
        for x in ix:
            W[x] = Wsymm[j]

    # Reproducibility:
    R = P.T @ W / ( 2*np.trace(D0) )

    # Rank deficiency:
    _, E, _ = svd( A )
    rank = len( np.where( E > 1E-6 )[0] ) 
    defi = n_g - rank
    print( 'Rank = %i    Deficiency = %i' % ( rank, defi ) )
    
    return ( W, R, defi )
# enddef expandDM_symm()


# Automated expansion of a given density matrix into sum of component density
# matrices
#
# D0:   The density matrix to be expanded
# D:    List of component density matrices (i.e., the basis)
# P:    Projection vector of D[] into D0
#
# NOTE: D must be pre-arranged in descending order of P
#
def expandDM_auto( D0, D, P, maxLSPercent, fixedLS, excludedLS, \
                   prec, degCri, isNormalized ):

    fixedLS = np.array( fixedLS ) - 1
    excludedLS = np.array( excludedLS ) - 1

    ix_gr = degenerateGroups( P, degCri*np.trace(D0)*2 )
    #print( ix_gr )
    n_gr = len( ix_gr )
    print( '%i degenerate groups of Lewis structures detected, ' 
            % n_gr, end='' )
    n_gr = round( n_gr * maxLSPercent );
    ix_gr = [ ix_gr[i] for i in range( n_gr ) ]
    print( 'only %i to be considered' % n_gr )

    # Reduce the Lewis strcutures set:
    ix_LS = []
    for ix in ix_gr:
        for j in ix:
            ix_LS.append( j )
    Dr = [ D[i] for i in ix_LS ] # D is list
    Pr = P[ ix_LS ] # P is numpy.ndarray

    # Calculate overlap matrix:
    print( 'Calculating overlap matrix between the Lewis structure DMs ...' )
    Sr, _ = overlapDM_Lewis( Dr )


    # Determine which groups of Lewis structures are essential (linearly
    # independent of all others):
    #for ix in ix_gr:
    #    print( ix )
    #    isEssentialLewis( Sr, ix_gr, ix )
    #    input()

    # Exclude excludedLS[]:
    ix_gr_trial = list( range( n_gr ) )
    j_removed = []
    for j in range( n_gr-1, -1, -1 ):
        # Remove the whole group if it contains any element of excludedLS[]:
        if np.intersect1d(ix_gr[j],excludedLS).size > 0:
            ix_gr_trial.remove( j )
            j_removed.append( j )
            print( 'Excluded Lewis structures: ', end='' )
            for k in ix_gr[j]:
                print( ' %i' % (k+1), end='' )
            print()

    # Check if the LS set is complete:
    W, R, defi  = expandDM_symm( D0, Sr, Pr, ix_gr, isNormalized )
    print( 'Reproducibility = %.3f%% for the initial set of Lewis structures'
            % (R*100) )
    if R < prec:
        print( 'Reproducibility is not sufficient ( < %.3f%%)' % (prec*100) )
        print( 'The Lewis structure (LS) set is incomplete.' )
        print( 'Try to increase the number of LSs '
                '(via parameter maxLSPercent or MAX_LP).' )
        exit( 1 )

    if defi == 0:
        print( 'Converged' )
        ix_gr_fin = [ ix_gr[i] for i in ix_gr_trial ]
        return ( W, R, ix_gr_fin )


    # Progressive elemination of dependent LSs:
    for j in range( n_gr-1, -1, -1 ):
        # Check if group j contains any atom in fixedLS:
        if np.intersect1d(ix_gr[j],fixedLS).size > 0:
            continue

        if j in j_removed:
            continue

        ix_gr_tmp = cp.deepcopy( ix_gr_trial )
        ix_gr_tmp.remove( j )
        ix_gr_curr = [ ix_gr[i] for i in ix_gr_tmp ]
        W, R, defi = expandDM_symm( D0, Sr, Pr, ix_gr_curr, isNormalized )
        print( 'Reproducibility = %.3f%% for the trial set of Lewis structures'
                % (R*100) )

        # Not sufficient ==> Keep the removed group of LSs:
        if R < prec:
            continue

        ix_gr_trial = ix_gr_tmp # Update ix_gr_trial

        # Converged:
        if defi == 0:
            print( 'Converged' )
            break
    
    ix_gr_fin = [ ix_gr[i] for i in ix_gr_trial ]
    return ( W, R, ix_gr_fin )
# enddef expandDM_auto()


# Density matrix based resonance theory (DMRT) analysis:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   atIx:         Indices of all atoms involved in the resonance analysis
#   Etot_LMO:     Total energy of LMOs == sum of energies of LMOs
#   MAX_LP:       Maximum number of lone pairs to be considered
#   maxLSPercent: Maximum percentage of Lewis structures to be considered
#   fixedLS:      Indices of Lewis structures that must not be eliminated
#   excludedLS:   Indices of Lewis structures that must be excluded
#   prec:         Minimum reproducibility that must be achieved
#   degCri:       Criterium to determine degenerate Lewis structures
#   raoType:      'uni'--> universal RAOs;  'var'--> varying RAOs
#   isNormalized: If impose the constraint that the sum of weights equals one
def dmrt( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
          MAX_LP=inf, maxLSPercent=1.0, fixedLS=[], excludedLS=[], \
          prec=0.99, degCri=1E-3, raoType='uni', isNormalized=True ):
    # Vectorize atIx and make moIx strat from ZERO (not 1):
    atIx = np.array( atIx1 ) # Starting from 1
    moIx = np.array( moIx1 ) - 1 # Starting from 0

    # Number of atoms:
    NAt = len( atIx )
    # Number of electrons:
    NE = len( moIx ) * 2  # Closed-shell system

    # Choose all LMOs involved in the resoncance subsystem:
    C_LMO = CLMO[:,moIx]

    if type( naoInfo ) is list: # Interface for dmrt_hmo()
        ifHMO = True
        aoIx = naoInfo
    else:                       # Normal usage
        ifHMO = False
        aoIx = naoInfo.naoIx

    # Density matrix from pi-LMOs in NAO basis:
    D0 = densityMatrixLMO( C_LMO )
    # Total energy of the LMOs for the resonance subsystem:
    Etot_LMO = LMOEnergy( ELMO, moIx )
    # Fock matrix in the original AO basis:
    F = FockMatrix( FchkInfo.C, FchkInfo.E )
    # Fock matrix in NAO basis:
    FNAO = CAONAO.T @ F @ CAONAO

    # ---- Gross atomic population and bond orders: ----#
    q, BO = popDM( D0, aoIx )
    # Print out atomic population:
    NAtAll = len( q )
    print( '='*80 )
    print( 'Gross atomic population from the LMOs:' )
    print()
    for i in range( 0, NAtAll, 10 ):
        kmax = min( i+10, NAtAll )
        for k in range( i, kmax ):
            print( ' %7s' % (k+1), end='' )
        print()
        for k in range( i, kmax ):
            print( ' %7s' % ('-'*7), end='' )
        print()
        for k in range( i, kmax ):
            print( ' %7.3f' % q[k], end='' )
        print()
        print()
    print( '='*80 )
    # Print out bond orders and valences:
    print( '='*80 )
    print( 'Wiberg bond orders and valences from the LMOs:' )
    print()
    for j in range( NAtAll ):
        for i in range( 0, NAtAll, 10 ):
            kmax = min( i+10, NAtAll )
            print( '%4s ' % ' ', end='' )
            for k in range( i, kmax ):
                print( ' %6s' % (k+1), end='' )
            print()
            print( '%4s ' % ' ', end='' )
            for k in range( i, kmax ):
                print( ' %6s' % ('-'*6), end='' )
            print()
            print( '%-4i ' % (j+1), end='' )
            for k in range( i, kmax ):
                print( ' %6.3f' % BO[j,k], end='' )
            print()
        print()
    print( '='*80 )


    # Get all resonance structures:
    LP0, BD0 = lewis_read( NAt, NE, MAX_LP )

    # Set actual indices of atoms:
    atIx = np.array( atIx )
    LP = list( map(lambda x: atIx[x-1], LP0 ) )
    BD = []
    for arr in BD0:
        BD.append(  np.array( list( map(lambda x: atIx[x-1], arr) ) )  )

    NLS = len( LP )  # Total number of Lewis structures


    # Resonance atomic orbitals (RAOs) in AO basis (NAOs, Lowdin's AOs etc.):
    D = [] # Density matrix of RAOs in NAO (or Lowdin etc.) basis
    E = [] # Energies of Lewis structures
    if raoType == 'uni': # Universal RAOs that are independent of LSs
        CRAO = rao_uni( D0, aoIx )
        for k in range( NLS ):
            D.append( DM_RAO_uni( CRAO, LP[k], BD[k] ) )

        #Alt: Calculate energies of Lewis structures:
        #Alt: CAORAO = CAONAO @ CRAO # RAO coefficient matrix in original AO basis
        #Alt: F_RAO = CAORAO.T @ F @ CAORAO # Fock matrix in RAO basis
        #Alt:for k in range( NLS ):
            #Alt: E.append( LewisEnergy( F_RAO, LP[k], BD[k], 'uni' ) )
    elif raoType == 'var': # Varying RAOs that depend on specific LSs:
        for k in range( NLS ):
            CRAOk, Dk = rao( D0, aoIx, LP[k], BD[k] )
            D.append( Dk )
            #Alt: Calculate energy of this Lewis structure:
            #Alt: CAORAOk = CAONAO @ CRAOk # RAO coeff. matrix in original AO basis
            #Alt: F_RAOk = CAORAOk.T @ F @ CAORAOk # Fock matrix in RAO basis
            #Alt: E.append( LewisEnergy( F_RAOk, LP[k], BD[k], 'var' ) )
    elif raoType == 'uni-gen': # General universal RAOs
        CRAO, raoIx = rao_uni_general( D0, aoIx )
        DM_Ix = [] # A list to record indices of DM for each Lewis structure
        iDM = 0
        for k in range( NLS ):
            Dk_list = DM_RAO_uni_general( CRAO, raoIx, LP[k], BD[k] )
            for Dk in Dk_list:
                D.append( Dk )
            DM_Ix.append( range( iDM, iDM+len(Dk_list) ) )
            iDM += len(Dk_list)
    else:
        raise ValueError( 
                'Unrecognized value for argument raoType in function dmrt()' )

    # Total number of density matrix:
    NDM = len( D )  # Note: NDM == NLS for mono-RAO-basis cases

    # Energies of Lewis structures:
    for k in range( NDM ):
        E.append( np.trace( D[k] @ FNAO ) )
    E = np.c_[ E ]

    # Projections onto the actual DM:
    P = projectDM( D, D0 )

    if raoType == 'uni-gen':
        print( 'Expanding the actual density matrix '
               'using multi-basis RAOs (%i in total)...' % len(D) )
        # Calculate overlap matrix:
        print( 'Calculating overlap matrix between Lewis structure DMs ...' )
        S, _ = overlapDM_Lewis( D )
        # Expand D0 into sum of D[]:
        W, R = expandDM( D0, S, P, True )
        print( 'The final reproducibility is %.6f%%' % (R*100) )
        print( 'Sum of weights is %.6f' % W.sum() )

        # Sum of Lewis structure energies, compared with the LMO total energy:
        Etot_LS = ( W * E ).sum()
        print( 'Total HF energy of reference structure is %.8E a.u.' 
                % Etot_LMO )
        print( 'Sum of energies of all %i Lewis structures is %.8E a.u.' 
                % ( len(D), Etot_LS ) )

        for i in range( NLS ):
            str_LS = lewis_str( LP[i], BD[i] ) 
            for ix in DM_Ix[i]:
                print( '%-5i %15.10f    %s' % (i+1, W[ix], str_LS) )
            print( '-'*80 )


        # Energies of Lewis structures:
        # NOTE: This part is theoretically not reasonable since the expansion
        # is not EXACT. As a result, the energy as well as deviation of DM is
        # not even higher (larger) when using multi-RAO-basis than that when
        # using a single principal RAO basis.
        #E_LS = np.empty( NLS )
        #R_LS = np.empty( NLS )
        #iLS = 0
        #for ix in DM_Ix:
        #    W_LS, R_LS[iLS] = expandDM( D0, S[np.ix_( ix,ix )], 
        #                                P[ix], True )
        #    #print( W_LS.sum() )
        #    E_LS[iLS] = ( W_LS * E[ix] ).sum()
        #    print( E_LS[iLS], R_LS[iLS]  )
        #    iLS += 1
        #Erel_LS = au2kcal( E_LS ) - Etot_LMO
        #print( Erel_LS )
        #print( 'Pauling--Wheland resonance energy is %.2f kcal/mol' % Ereson )

        return

    # Sort by P:
    ix_sort = np.argsort( -P[:,0] )
    P = P[ ix_sort ]

    # Update LP[], BD[], D and E:
    LP = [ LP[i] for i in ix_sort ]
    BD = [ BD[i] for i in ix_sort ]
    D = [ D[i] for i in ix_sort ]
    E = [ E[i] for i in ix_sort ]

    # Lewis structure energies relative to the LMO total energy:
    if ifHMO:
        Erel = E - Etot_LMO
    else:
        Erel = au2kcal( E - Etot_LMO )
    # The lowest Erel is the Pauling--Wheland resonance energy
    Ereson = min( Erel )
    ix_min_LS = []
    for k in range( NLS ):
        if abs( Erel[k] - Ereson ) < 0.001:
            ix_min_LS.append( k+1 )

    if ifHMO:
        print( 'Total HF energy of reference structure is %.4f |beta|' 
                % Etot_LMO )
        print( 'Pauling--Wheland resonance energy is %.4f |beta|' % Ereson )
    else:
        print( 'Total HF energy of reference structure is %.8E a.u.' 
                % Etot_LMO )
        print( 'Pauling--Wheland resonance energy is %.2f kcal/mol' % Ereson )
    if len( ix_min_LS ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_LS:
        print( ' #%i' % i, end='' )
    print()

    # Get strings for Lewis structures:
    print( 'The complete set of Lewis structures for a %ic-%ie system:' %
            ( NAt, NE ) )
    NE2 = NE * 2
    STR = []
    for i in range( NLS ):
        STR.append( lewis_str( LP[i], BD[i] ) )
        print( '%-5i %15.10f %8.2f    %s' % (i+1, P[i]/NE2, Erel[i], STR[i]) )
    print()

    # Automated expansion of density matrix:
    W, R, ix_gr = expandDM_auto( D0, D, P, maxLSPercent, fixedLS, excludedLS, \
                                 prec, degCri, isNormalized )
    print( 'The final reproducibility is %.3f%%' % (R*100) )

    Wsum = W[ np.where(W>0) ].sum()

    # Print out the analysis result:
    print( 'Final result:' )
    print( '-'*80 )
    Etot_LS = 0.
    Wsum_all = 0.
    for ix in ix_gr:
        for i in ix:
            if W[i] > 0:
                weight = '%5.2f%%' % ( W[i]*100 / Wsum )
            else:
                weight = ' -----'
            print( '%-5i %15.10f %8.2f ( %9.6f )  [ %s ]    %s' 
                    % ( i+1, P[i]/NE2, Erel[i], W[i], weight, STR[i] ) )
            Etot_LS += W[i]*E[i]
            Wsum_all += W[i]
    print( '-'*80 )
    print( 'Reproducibility = %10.6f%%' % ( R*100 ) )
    Esum = Etot_LS
    print( 'Sum of all weights is %.8f' % Wsum_all )

    if ifHMO:
        print( 'Total HF energy of reference structure is %.4f |beta|' 
                % Etot_LMO )
        print( 'Total energy of expanded wave function is %.4f |beta|' 
                % Esum )
        print( '  Percent error = %.6f%% ' 
                % ( ( Esum/Etot_LMO - 1 ) * 100 * np.sign( Etot_LMO ) ) )
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.4f |beta|' % Ereson )
        print( 'Pauling--Wheland resonance energy is %.4f |beta|' % Ereson )
        Ereson_expan = min(E) - Esum
    else:
        print( 'Total HF energy of reference structure is %.8E a.u.' 
                % Etot_LMO )
        print( 'Total energy of expanded wave function is %.8E a.u.' 
                % Esum )
        print( '  Percent error = %.6f%% ' 
                % ( ( Esum/Etot_LMO - 1 ) * 100 * np.sign( Etot_LMO ) ) )
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.2f kcal/mol' % Ereson )
        print( 'Pauling--Wheland resonance energy is %.2f kcal/mol' % Ereson )
        Ereson_expan = au2kcal( min(E) - Esum )
    print( '  Percent error = %.6f%% ' % ( ( 1 - Ereson_expan/Ereson )*100 ) )
    if len( ix_min_LS ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_LS:
        print( ' #%i' % i, end='' )
    print()

    return
# enddef dmrt()


# DMRT analysis with a specified set of Lewis structures:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   atIx:         Indices of all atoms involved in the resonance analysis
#   Etot_LMO:     Total energy of LMOs == sum of energies of LMOs
#   ixLS:         Lewis structures specified in either of the following ways:
#                 (1) List of indices of the specified Lewis structures
#                 This choice involves enumeration of all possbile LSs. If the
#                 enumeration has been done previously, file LEWIS_Nc_Me.dat
#                 will be read in.
#                 NOTE: 1) Indices start from 1 (not 0)
#                       2) The list is from ordering of projections, determined
#                          by calling the dmrt() function prior to calling this
#                          function
#                 (2) Straightfoward specification given by LP[] and BD[],
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[] and BD[], atomic indices start from 1 
#                       (not 0), and they correspond to ACTUAL indices of 
#                       atoms, which may not be consecutive, not as those in 
#                       LP0[] and BD0[] (also see notes in body of this func)
#   raoType:      'uni'--> universal RAOs;  'var'--> varying RAOs
def dmrt_spec( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, ixLS, \
               raoType='uni' ):

    # Vectorize atIx and make moIx strat from ZERO (not 1):
    atIx = np.array( atIx1 ) # Starting from 1
    moIx = np.array( moIx1 ) - 1 # Starting from 0

    # Number of atoms:
    NAt = len( atIx )
    # Number of electrons:
    NE = len( moIx ) * 2  # Closed-shell system

    # Choose all LMOs involved in the resoncance subsystem:
    C_LMO = CLMO[:,moIx]

    if type( naoInfo ) is list: # Interface for dmrt_hmo()
        ifHMO = True
        aoIx = naoInfo
    else:                       # Normal usage
        ifHMO = False
        aoIx = naoInfo.naoIx

    # Density matrix from pi-LMOs in NAO basis:
    D0 = densityMatrixLMO( C_LMO )
    # Total energy of the LMOs for the resonance subsystem:
    Etot_LMO = LMOEnergy( ELMO, moIx )
    # Fock matrix in the original AO basis:
    F = FockMatrix( FchkInfo.C, FchkInfo.E )
    # Fock matrix in NAO basis:
    FNAO = CAONAO.T @ F @ CAONAO

    # ---- Gross atomic population and bond orders: ----#
    q, BO = popDM( D0, aoIx )
    # Print out atomic population:
    NAtAll = len( q )
    print( '='*80 )
    print( 'Gross atomic population from the LMOs:' )
    print()
    for i in range( 0, NAtAll, 10 ):
        kmax = min( i+10, NAtAll )
        for k in range( i, kmax ):
            print( ' %7s' % (k+1), end='' )
        print()
        for k in range( i, kmax ):
            print( ' %7s' % ('-'*7), end='' )
        print()
        for k in range( i, kmax ):
            print( ' %7.3f' % q[k], end='' )
        print()
        print()
    print( '='*80 )
    # Print out bond orders and valences:
    print( '='*80 )
    print( 'Wiberg bond orders and valences from the LMOs:' )
    print()
    for j in range( NAtAll ):
        for i in range( 0, NAtAll, 10 ):
            kmax = min( i+10, NAtAll )
            print( '%4s ' % ' ', end='' )
            for k in range( i, kmax ):
                print( ' %6s' % (k+1), end='' )
            print()
            print( '%4s ' % ' ', end='' )
            for k in range( i, kmax ):
                print( ' %6s' % ('-'*6), end='' )
            print()
            print( '%-4i ' % (j+1), end='' )
            for k in range( i, kmax ):
                print( ' %6.3f' % BO[j,k], end='' )
            print()
        print()
    print( '='*80 )


    #-------------------- Lewis structures ----------------------------------
    #-- OPTION 1: Automatic enumeration ----------
    if type(ixLS) is not tuple:
        ifEnume = True
        # Get all resonance structures:
        #  LP0:   List of LPs for all considered Lewis structures
        #  BD0:   List of BDs for all considered Lewis structures
        #  Suffix 0 means the indices of atoms stored in LP0[] and BD0[]
        #  follow the order of appearance, i.e., 1, 2, 3, ..., which is NOT 
        #  necessarily the actural indices of atoms appeared in the given 
        #  molecule. The latter indices are stored in LP[] and BD[]
        LP0, BD0 = lewis_read( NAt, NE, inf )

        # Set actual indices of atoms:
        atIx = np.array( atIx )
        LP = list( map(lambda x: atIx[x-1], LP0 ) )
        BD = []
        for arr in BD0:
            BD.append(  np.array( list( map(lambda x: atIx[x-1], arr) ) )  )

    #-- OPTION 2: Straightforward specification ----------
    else:
        ifEnume = False
        LP, BD = ixLS

        # Set sequential indices of appearance of atoms:
        atIx = np.array( atIx, dtype='i' )
        sqIx = -np.ones( max(atIx), dtype='i' )
        for i in range( NAt ):
            sqIx[ atIx[i]-1 ] = i
        LP0 = list( map(lambda x: sqIx[x-1]+1, LP ) )
        BD0 = []
        for arr in BD:
            BD0.append(  np.array( list( map(lambda x: sqIx[x-1]+1, arr) ), \
                                    dtype='i' )  )
        # Make sure that indices of two atoms in each bond are in ascending
        # order:
        BD0tmp = []
        for bds in BD0:
            for i in range( len(bds) ):
                if bds[i,0] > bds[i,1]:
                    bds[i,:] = [ bds[i,1], bds[i,0] ]
            BD0tmp.append( bds )
        BD0 = BD0tmp
    #------------------------------------------------------------------------

    NLS = len( LP )  # Total number of Lewis structures
    print( '%i Lewis structures in total to be considered' % NLS )

    # Resonance atomic orbitals (RAOs) in AO basis (NAOs, Lowdin's AOs etc.):
    D = [] # Density matrix of RAOs in NAO (or Lowdin etc.) basis
    if raoType == 'uni': # Universal RAOs that are independent of LSs
        CRAO = rao_uni( D0, aoIx )
        for k in range( NLS ):
            D.append( DM_RAO_uni( CRAO, LP[k], BD[k] ) )
    elif raoType == 'var': # Varying RAOs that depend on specific LSs:
        for k in range( NLS ):
            CRAOk, Dk = rao( D0, aoIx, LP[k], BD[k] )
            D.append( Dk )
    elif raoType == 'uni-gen': # General universal RAOs
        CRAO, raoIx = rao_uni_general( D0, aoIx )
        DM_Ix = [] # A list to record indices of DM for each Lewis structure
        iDM = 0
        for k in range( NLS ):
            Dk_list = DM_RAO_uni_general( CRAO, raoIx, LP[k], BD[k] )
            for Dk in Dk_list:
                D.append( Dk )
            DM_Ix.append( range( iDM, iDM+len(Dk_list) ) )
            iDM += len(Dk_list)
    else:
        raise ValueError( 
                'Unrecognized value for argument raoType in function dmrt()' )

    # Total number of density matrix:
    NDM = len( D )  # Note: NDM == NLS for mono-RAO-basis cases

    # Energies of Lewis structures:
    E = [] # Energies of Lewis structures
    for k in range( NDM ):
        E.append( np.trace( D[k] @ FNAO ) )
    E = np.array( E )

    # Projections onto the actual DM:
    P = projectDM( D, D0 )

    if raoType == 'uni-gen':
        print( 'Expanding the actual density matrix '
               'using multi-basis RAOs (%i in total)...' % len(D) )
        # Calculate overlap matrix:
        print( 'Calculating overlap matrix between Lewis structure DMs ...' )
        S, _ = overlapDM_Lewis( D )
        # Expand D0 into sum of D[]:
        W, R = expandDM( D0, S, P, True )
        print( 'The final reproducibility is %.6f%%' % (R*100) )
        print( 'Sum of weights is %.6f' % W.sum() )

        # Sum of Lewis structure energies, compared with the LMO total energy:
        Etot_LS = ( W * E ).sum()
        print( 'HF Energy of reference structure is %.3f kcal/mol' % Etot_LMO )
        print( 'Total HF energy of reference structure is %.8E a.u.'
                % Etot_LMO )
        print( 'Sum of energies of all %i Lewis structures is %.8E a.u.' 
                % ( len(D), Etot_LS ) )

        for i in range( NLS ):
            str_LS = lewis_str( LP[i], BD[i] ) 
            for ix in DM_Ix[i]:
                print( '%-5i %15.10f    %s' % (i+1, W[ix], str_LS) )
            print( '-'*80 )

        return

    # Sort by P:
    ix_sort = np.argsort( -P[:,0] )
    P = P[ ix_sort ]

    # Update LP[], BD[], D and E:
    LP = [ LP[i] for i in ix_sort ]
    BD = [ BD[i] for i in ix_sort ]
    D = [ D[i] for i in ix_sort ]
    E = E[ ix_sort ]

    # Lewis structure energies relative to the LMO total energy:
    if ifHMO:
        Erel = E - Etot_LMO
    else:
        Erel = au2kcal( E - Etot_LMO )
    # The lowest Erel is the Pauling--Wheland resonance energy
    Ereson = min( Erel )
    ix_min_LS = []
    for k in range( NLS ):
        if abs( Erel[k] - Ereson ) < 0.001:
            ix_min_LS.append( k+1 )

    if ifHMO:
        print( 'Total HF energy of reference structure is %.4f |beta|' 
                % Etot_LMO )
        print( 'Pauling--Wheland resonance energy is %.4f |beta|' % Ereson )
    else:
        print( 'Total HF energy of reference structure is %.8E a.u.' 
                % Etot_LMO )
        print( 'Pauling--Wheland resonance energy is %.2f kcal/mol' % Ereson )
    if len( ix_min_LS ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_LS:
        print( ' #%i' % i, end='' )
    print()

    # Get strings for Lewis structures:
    print( 'The complete set of Lewis structures for a %ic-%ie system:' %
            ( NAt, NE ) )
    NE2 = NE * 2
    str_LS = []
    for i in range( NLS ):
        str_LS.append( lewis_str( LP[i], BD[i] ) )
        print( '%-5i %15.10f %8.2f    %s' 
                % (i+1, P[i]/NE2, Erel[i], str_LS[i]) )
    print()

    # Only choose the specified Lewis structures:
    if ifEnume:
        ix_reducLS = np.array( ixLS, dtype='i' ) - 1
        print( 'Using the %i specified Lewis structures to expand the '
               'actual wave function' % ( len( ix_reducLS ) ) )
    else:
        ix_reducLS = np.arange( NLS )
        print( 'Using %i straightforwardly specified Lewis structures to '
               'expand the actual wave function' % ( len( ix_reducLS ) ) )

    # Reduce D[] and P[]:
    Dr = [ D[i] for i in ix_reducLS ] # D is list
    Pr = P[ ix_reducLS ] # P is numpy.ndarray

    # Calculate overlap matrix:
    print( 'Calculating overlap matrix between the Lewis structure DMs ...' )
    Sr, rank = overlapDM_Lewis( Dr )
    rank_defi = len( ix_reducLS ) - rank
    if rank_defi > 0:
        print( 'Warning: The input Lewis strucutres are not all linearly'
                ' independent\nRank deficiency is %i' % rank_defi )
    #if rank_defi > 0:
    #    raise AssertionError( 'The input Lewis strucutres are not all linearly'
    #                 ' independent\nRank deficiency is %i' % rank_defi )
        #print( 'WARNING: ', end='' )
        #print( 'The input Lewis strucutres are not all linearly'
        #             ' independent\nRank deficiency is %i' % rank_defi )


    Wr, R  = expandDM( D0, Sr, Pr, True )
    W2r = (Wr*Wr) / (Wr*Wr).sum()

    # Recover the full coefficient vector by setting zeros for unchosen LSs:
    W = np.zeros( NLS ) # Initialize
    W2 = np.zeros( NLS ) # Initialize
    for i in range( NLS ):
        ix_reduc = np.nonzero( ix_reducLS == i )[0]
        if len( ix_reduc ) > 0:
            W[i] = Wr[ ix_reduc ]
            W2[i] = W2r[ ix_reduc ]

    # Energy of the admixture of Lewis structures (resonance hybrid):
    Etot_LS = ( W * E ).sum()


    #------------------------------------------------------------------
    # Print out the final results:
    #------------------------------------------------------------------
    Pn = P / ( 2 * NE ) # Normalized projections
    print( '-'*80 )
    print( '%5s %10s %8s %11s' % 
            ( 'No.', 'Projection', 'Erel', 'Weight' ), end='' )
    print( ' %10s' % 'Weight^2%', end='' )
    print( '  Lewis structure' )
    print( '-'*80 )
    for i in range( NLS ):
        # Only show Lewis structures in the reduced set:
        if np.count_nonzero( ix_reducLS == i ) == 0:
            continue
        print( '%5i %10.6f %8.2f %11.7f' 
                % ( i+1, Pn[i], Erel[i], W[i] ), end='' )
        print( ' %9.2f%%' % ( W2[i]*100 ), end='' )
        print( '  %s' % str_LS[i] )
    print( '-'*80 )
    print( 'Reproducibility = %10.6f%%' % ( R*100 ) )
    Esum = Etot_LS
    if ifHMO:
        print( 'Total HF energy of reference structure is %.4f |beta|' 
                % Etot_LMO )
        print( 'Total energy of expanded wave function is %.4f |beta|' 
                % Esum )
        print( '  Percent error = %.6f%% ' 
                % ( ( Esum/Etot_LMO - 1 ) * 100 * np.sign( Etot_LMO ) ) )
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.4f |beta|' % Ereson )
        print( 'Pauling--Wheland resonance energy is %.4f |beta|' % Ereson )
        Ereson_expan = min(E) - Esum
        print( 'Pauling--Wheland resonance energy given by expanded wave '
               'function is %.4f |beta|' % Ereson_expan )
        print( '  Percent error = %.6f%% ' 
                % ( ( 1 - Ereson_expan/Ereson )*100 ) )
    else:
        print( 'Total HF energy of reference structure is %.8E a.u.' 
                % Etot_LMO )
        print( 'Total energy of expanded wave function is %.8E a.u.' 
                % Esum )
        print( '  Percent error = %.6f%% ' 
                % ( ( Esum/Etot_LMO - 1 ) * 100 * np.sign( Etot_LMO ) ) )
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.2f kcal/mol' % Ereson )
        print( 'Pauling--Wheland resonance energy is %.2f kcal/mol' % Ereson )
        Ereson_expan = au2kcal( min(E) - Esum )
        print( 'Pauling--Wheland resonance energy given by expanded wave '
               'function is %.2f kcal/mol' % Ereson_expan )
        print( '  Percent error = %.6f%% ' 
                % ( ( 1 - Ereson_expan/Ereson )*100 ) )
    if len( ix_min_LS ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_LS:
        print( ' #%i' % i, end='' )
    print()
    #------------------------------------------------------------------

    return
# enddef dmrt_spec()



# DMRT analysis in the simple Hueckel Molecular Orbital (HMO) framework:
#   hmoSol:      Solution of simple HMO theory
#   MAX_LP:       Maximum number of lone pairs to be considered
#   maxLSPercent: Maximum percentage of Lewis structures to be considered
#   fixedLS:      Indices of Lewis structures that must not be eliminated
#   excludedLS:   Indices of Lewis structures that must be excluded
#   prec:         Minimum reproducibility that must be achieved
#   degCri:       Criterium to determine degenerate Lewis structures
#   raoType:      'uni'--> universal RAOs;  'var'--> varying RAOs
#   isNormalized: If impose the constraint that the sum of weights equals one
def dmrt_hmo( hmoSol, MAX_LP=inf, maxLSPercent=1., fixedLS=[], excludedLS=[], \
              prec=0.99, degCri=1E-3, isNormalized=True ):
    NAt = hmoSol.NAt
    NP = hmoSol.NE // 2  # Number of electron pairs (for closed-shell systems)
    aoIx = []
    for i in range( 1, NAt+1 ):
        aoIx.append( np.array([ i ]) )

    CAONAO = np.eye( NAt )
    CLMO = hmoSol.C[ :, 0:NP ]
    ELMO = hmoSol.E[ 0:NP ]
    atIx1 = range( 1, NAt+1 )
    moIx1 = range( 1, NP+1 )
    raoType = 'uni'

    dmrt( aoIx, hmoSol, CAONAO, CLMO, ELMO, atIx1, moIx1, MAX_LP, \
          maxLSPercent, fixedLS, excludedLS, prec, degCri, raoType, \
          isNormalized )
    return
# enddef dmrt_hmo()



# DMRT analysis in the simple Hueckel Molecular Orbital (HMO) framework with 
# a specified set of Lewis structures:
#   hmoSol:      Solution of simple HMO theory
#   ixLS:         Lewis structures specified in either of the following ways:
#                 (1) List of indices of the specified Lewis structures
#                 This choice involves enumeration of all possbile LSs. If the
#                 enumeration has been done previously, file LEWIS_Nc_Me.dat
#                 will be read in.
#                 NOTE: 1) Indices start from 1 (not 0)
#                       2) The list is from ordering of projections, determined
#                          by calling the dmrt_hmo() function prior to calling 
#                          this function
#                 (2) Straightfoward specification given by LP[] and BD[],
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[] and BD[], atomic indices start from 1 
#                       (not 0), and they correspond to ACTUAL indices of 
#                       atoms, which may not be consecutive, not as those in 
#                       LP0[] and BD0[] (also see notes in body of this func)
def dmrt_hmo_spec( hmoSol, ixLS ):
    NAt = hmoSol.NAt
    NP = hmoSol.NE // 2  # Number of electron pairs (for closed-shell systems)
    aoIx = []
    for i in range( 1, NAt+1 ):
        aoIx.append( np.array([ i ]) )

    CAONAO = np.eye( NAt )
    CLMO = hmoSol.C[ :, 0:NP ]
    ELMO = hmoSol.E[ 0:NP ]
    atIx1 = range( 1, NAt+1 )
    moIx1 = range( 1, NP+1 )
    raoType = 'uni'

    dmrt_spec( aoIx, hmoSol, CAONAO, CLMO, ELMO, atIx1, moIx1, ixLS, raoType )
    return
# enddef dmrt_hmo_spec()
