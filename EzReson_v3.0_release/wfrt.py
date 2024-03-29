# wave-function based resonance theory (wfrt)
#
#  Created on Apr 7, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    May 16, 2022:
#    - In functions wfrt_spec(), wfrt_spec_spin(), wfrt_hmo_spec(),
#    wfrt_hmo_spec_spin(), added the `sciform` option
#
#    May 12, 2022:
#    - Added function wfrt_spec_spin() as an extension of wfrt_spec() for
#    spin-unrestricted cases
#
#    May 1, 2022:
#    - Added function proj_DM_WF_spin() as an extension of proj_DM_WF() for
#    spin-unrestricted cases
#
#    Feb 18, 2022:
#    - Changed PWSO scheme to P2WSO, for the WFRT weight is approx.
#    proportional to the square of WF-projection, not the projection.
#    More importantly, for open-shell systems, the projection of resonance 
#    structure onto the actual structure can be negative. Thus, the squared
#    form seems more reasonable and natural.
#
#    Feb 18, 2022:
#    - Extended functions wfrt_kekule() to radical systems also, and for this
#    purpose, a new function wfrt_spec_spin() was created
#    - Allowed to write RAOs to a fchk file when calling wfrt_kekule() by
#    assigning the keyword "writeRAOs = TRUE" in the input file
#
#    Dec 11, 2020:
#    - Fixed a bug in WFRT-HMO analysis for fragmented pi-conjugate systems.
#    The code in function rao_uni() fails. The solution is to directly set the
#    RAO as the identity matrix without the need of calling rao_uni().
#
#    Sep 1, 2020:
#    - Added function wfrt_hmo_kekule() to perform WFRT analysis in the HMO
#      framework using only Kekule structures
#    - Added in proj_DM_WF(), proj_DM_WF_kekule() the argument raoFlipAtoms 
#      to allow one to manually assign the RAO phases
#
#    Aug 28, 2020:
#    - Fixed a bug in proj_DM_WF() that when ixLS==-1, only 1 Lewis structure
#      is considered
#
#    Jul 29, 2020:
#    - Added an option in wfrt_spec() to manually assign RAO phases
#    - Added an option in wfrt_spec() to output RAOs into a *.fchk file
#
#    Jul 23, 2020:
#    - Added an option in wfrt() to manually assign RAO phases
#    - Added an option in wfrt() to output RAOs into a *.fchk file
#
#    From May 17, 2020:
#    - Modified function wfrt_spec() so that for cases where the specified
#      Lewis structures have linear dependence, the analysis is still carried
#      out, but the weights are not calculated.
#
#    From Apr 28, 2020:
#    - Added implementation using wave functions obtained from simple Hueckel  
#      Molecular Orbital theory
#
#    Apr 25, 2020:
#      - In function proj_DM_WF(), allow for straightforward specification
#        of Lewis structures parsed from strings, or for a given set of Kekule 
#        structures read from an external *.kek file geneated by enumKekule()
#
#    Apr 23, 2020:
#      - Allow for straightforward specification of Lewis structures parsed
#        from strings, in function wfrt_spec()
#
#      - WFRT analysis for a given set of Kekule structures read from an
#        external *.kek file, which is geneated by enumKekule() function
#
#    From Apr 18, 2020:
#      - Calculations of energies of resonance strcutures and their interaction
#        energies
#

import globals
from readNBOMat import *
from resonance import *
from weights import *
from kekule import *


# Lewis unit orbitals (LUOs) in AO basis (NAOs, Lowdin's AOs etc.):
#   LUOs are composed of RAOs and have two types:
#     1) LP-LUOs, corresponding to a lone pair localized at an atom, equal
#        to the RAO of that atom
#     2) BD-LUOs, corresponding to a pure covalent bond between two atoms,
#        being a linear combination (50-50) of the RAOs of both atoms
#
#  CRAO: Resonance atomic orbitals in orthonormal AO basis (NAOs, Lowdin's AOs)
#  FNAO: Fock matrix in NAO basis
#  atIx: Indices of all atoms involved in the resonance analysis (start from 1)
#  NAt:  Number of atoms/centers in the subsystem, to which WFRT is applied
#  NE:   Number of electrons in the subsystem, to which WFRT is applied
def luo( CRAO, FNAO, atIx, NAt, NE ):
    # Total number of LUOs == number of all combinations of 2e pairs 
    #                      == C_{NAt}^{2} + C_{NAt}^{1}
    N_LUO = NAt * (NAt+1) // 2 
    print( '%i Lewis unit orbitals (LUOs) in total' % N_LUO )

    # Numeration of all possible LPs and BDs associated with the LUOs
    atIx_LUO = np.zeros( ( N_LUO, 2 ), dtype='i' ) # Inidices start from ZERO
    # 2c BDs go first:
    k = 0
    for i in range( NAt ):
        for j in range( i+1, NAt ):
            atIx_LUO[ k, : ] = [ i, j ]
            k += 1
    # 1c LPs go second:
    for i in range( NAt ): 
        atIx_LUO[ k, : ] = [ i, -1 ]
        k += 1
    #print( atIx_LUO )

    # Coefficient matrix of LUOs in terms of RAO, i.e., RAO -> LUO trans. mat.
    Q_LUO = np.zeros( ( NAt, N_LUO ) )
    for i in range( N_LUO ):
        if atIx_LUO[ i, 1 ] == -1: # LP
            Q_LUO[ atIx_LUO[i,0], i ] = 1.
        else: # BD
            Q_LUO[ atIx_LUO[i,:], i ] = 1. / np.sqrt(2)
    #writeMat( Q_LUO )

    # Overlap matrix between LUOs:
    S_LUO = Q_LUO.T @ Q_LUO  # Note that RAOs are orthonormal to each other
    #writeMat( S_LUO )

    # LUO coefficients in NAO basis only for involved atoms, atIx[]:
    C_RAO = CRAO[ :, atIx-1 ] # NOTE: Indices in atIx[] start from 1 (not 0)
    C_LUO = C_RAO @ Q_LUO #

    # Fock matrix between LUOs:
    F_LUO = C_LUO.T @ FNAO @ C_LUO
    #writeMat( F_LUO )
    #input()

    return ( atIx_LUO, C_LUO, S_LUO, F_LUO )
# enddef luo()


# Overlap matrix between a given set of Lewis structures:
#  NE: Number of electrons
#  luoIX: Matrix storing the indices of LUOs of the Lewis structures
#  S_LUO: Overlap matrix between LUOs:
def overlapLewis( NE, luoIX, S_LUO ):
    NLS = len( luoIX ) # Number of Lewis structures
    NP = NE // 2 # Number of electron pairs
    S = np.eye( NLS ) # Initialize overlap matrix
    A = np.zeros( ( NP, NP ) ) # Initialize single-spin determinant
    for i in range( NLS ):
        for j in range( i+1, NLS ):
            for i1 in range( NP ): # Running over all alpha electrons
                ix1 = luoIX[ i, i1 ] # Index of LUO-1
                for i2 in range( NP ): # Running over all alpha electrons
                    ix2 = luoIX[ j, i2 ] # Index of LUO-2
                    A[ i1, i2 ] = S_LUO[ ix1, ix2 ]
            S[ i, j ] = np.linalg.det( A ) ** 2
            S[ j, i ] = S[ i, j ]

    # Determine the rank of S:
    _, E, _ = svd( S )
    rank = len( np.where( E > 1E-6 )[0] ) 
    return ( S, rank )
# enddef overlapLewis()


# Overlap & Hamiltonian/Fock matrices between a given set of Lewis structures:
# NOTE: Current implement ONLY FOR CLOSED-SHELL systems
#
#  NE: Number of electrons
#  luoIX: Matrix storing the indices of LUOs of the Lewis structures
#  S_LUO: Overlap matrix between LUOs:
#  F_LUO: Fock matrix between LUOs:
def overlapFockLewis( NE, luoIX, S_LUO, F_LUO ):
    NLS = len( luoIX ) # Number of Lewis structures
    NP = NE // 2 # Number of electron pairs
    S = np.eye( NLS ) # Initialize overlap matrix
    F = np.zeros( ( NLS, NLS ) ) # Initialize Fock matrix
    A = np.zeros( ( NP, NP ) ) # Initialize single-spin determinant
    for i in range( NLS ):
        for j in range( i+1, NLS ):
            # Overlap matrix:
            for i1 in range( NP ): # Running over all alpha electrons
                ix1 = luoIX[ i, i1 ] # Index of LUO-1
                for i2 in range( NP ): # Running over all alpha electrons
                    ix2 = luoIX[ j, i2 ] # Index of LUO-2
                    A[ i1, i2 ] = S_LUO[ ix1, ix2 ]
            detA = np.linalg.det( A )
            S[ i, j ] = detA * detA
            S[ j, i ] = S[ i, j ]
            # Fock matrix:
            for i1 in range( NP ): # Running over all rows of B
                B = A.copy() # Initialize B as A
                ix1 = luoIX[ i, i1 ] # Index of LUO-1
                # Change in the current row the elements to Fock integrals:
                for i2 in range( NP ): # Running over all alpha electrons
                    ix2 = luoIX[ j, i2 ] # Index of LUO-2
                    B[ i1, i2 ] = F_LUO[ ix1, ix2 ]
                detB = np.linalg.det( B )
                F[ i, j ] += detA * detB
            F[ j, i ] = F[ i, j ]
        # Speical case when i == j for Fock matrix (i.e., diagonal elements):
        for i1 in range( NP ): # Running over all alpha electrons
            ix1 = luoIX[ i, i1 ] # Index of LUO-1
            F[ i, i ] += F_LUO[ ix1, ix1 ]
    F = F * 2 # For closed-shell system, adding alpha+beta is doubling alpha

    # Determine the rank of S:
    _, E, _ = svd( S )
    rank = len( np.where( E > 1E-6 )[0] ) 
    return ( S, F, rank )
# enddef overlapFockLewis()


# Spin-unrestricted version of overlapFockLewis()
# Overlap & Hamiltonian/Fock matrices between a given set of Lewis structures:
#
#  NEa: Number of alpha electrons
#  NEb: Number of beta electrons
#  luoIXa: Matrix storing the indices of LUOs of alpha elec. for Lewis str.
#  luoIXb: Matrix storing the indices of LUOs of beta elec. for Lewis str.
#  S_LUOa: Overlap matrix between alpha LUOs
#  S_LUOb: Overlap matrix between beta LUOs
#  F_LUOa: Fock matrix between alpha LUOs
#  F_LUOb: Fock matrix between beta LUOs
def overlapFockLewis_spin( NEa, NEb, luoIXa, luoIXb, \
        S_LUOa, S_LUOb, F_LUOa, F_LUOb ):
    NLS = len( luoIXa ) # Number of Lewis structures
    S = np.eye( NLS ) # Initialize overlap matrix
    F = np.zeros( ( NLS, NLS ) ) # Initialize Fock matrix
    A = np.zeros( ( NEa, NEa ) ) # Initialize alpha determinant
    B = np.zeros( ( NEb, NEb ) ) # Initialize beta determinant
    for i in range( NLS ):
        for j in range( i+1, NLS ):
            # Overlap matrix:
            # 1) Alpha part:
            for i1 in range( NEa ): # Running over all alpha electrons
                ix1 = luoIXa[ i, i1 ] # Index of LUO-1
                for i2 in range( NEa ): # Running over all alpha electrons
                    ix2 = luoIXa[ j, i2 ] # Index of LUO-2
                    A[ i1, i2 ] = S_LUOa[ ix1, ix2 ]
            detA = np.linalg.det( A )
            # 2) Beta part:
            for i1 in range( NEb ): # Running over all beta electrons
                ix1 = luoIXb[ i, i1 ] # Index of LUO-1
                for i2 in range( NEb ): # Running over all beta electrons
                    ix2 = luoIXb[ j, i2 ] # Index of LUO-2
                    B[ i1, i2 ] = S_LUOb[ ix1, ix2 ]
            detB = np.linalg.det( B )
            S[ i, j ] = detA * detB
            S[ j, i ] = S[ i, j ]

            # Fock matrix:
            # 1) Alpha part:
            for i1 in range( NEa ): # Running over all rows of AF
                AF = A.copy() # Initialize AF as A
                ix1 = luoIXa[ i, i1 ] # Index of LUO-1
                # Change in the current row the elements to Fock integrals:
                for i2 in range( NEa ): # Running over all alpha electrons
                    ix2 = luoIXa[ j, i2 ] # Index of LUO-2
                    AF[ i1, i2 ] = F_LUOa[ ix1, ix2 ]
                detAF = np.linalg.det( AF )
                F[ i, j ] += detAF * detB
            # 2) Beta part:
            for i1 in range( NEb ): # Running over all rows of BF
                BF = B.copy() # Initialize BF as B
                ix1 = luoIXb[ i, i1 ] # Index of LUO-1
                # Change in the current row the elements to Fock integrals:
                for i2 in range( NEb ): # Running over all beta electrons
                    ix2 = luoIXb[ j, i2 ] # Index of LUO-2
                    BF[ i1, i2 ] = F_LUOb[ ix1, ix2 ]
                detBF = np.linalg.det( BF )
                F[ i, j ] += detA * detBF
            F[ j, i ] = F[ i, j ]

        # Speical case when i == j for Fock matrix (i.e., diagonal elements):
        # 1) Alpha part:
        for i1 in range( NEa ): # Running over all alpha electrons
            ix1 = luoIXa[ i, i1 ] # Index of LUO-1
            F[ i, i ] += F_LUOa[ ix1, ix1 ]
        # 2) Beta part:
        for i1 in range( NEb ): # Running over all beta electrons
            ix1 = luoIXb[ i, i1 ] # Index of LUO-1
            F[ i, i ] += F_LUOb[ ix1, ix1 ]

    # Determine the rank of S:
    _, E, _ = svd( S )
    rank = len( np.where( E > 1E-6 )[0] ) 
    return ( S, F, rank )
# enddef overlapFockLewis_spin()


# Open-shell (spin-restricted) version of overlapFockLewis()
# Overlap & Hamiltonian/Fock matrices between a given set of Lewis structures:
#
#  NEa: Number of alpha electrons
#  NEb: Number of beta electrons
#  luoIXa: Matrix storing the indices of LUOs of alpha elec. for Lewis str.
#  luoIXb: Matrix storing the indices of LUOs of beta elec. for Lewis str.
#  S_LUO: Overlap matrix between LUOs
#  F_LUO: Fock matrix between LUOs
def overlapFockLewis_open( NEa, NEb, luoIXa, luoIXb, S_LUO, F_LUO ):
    NLS = len( luoIXa ) # Number of Lewis structures
    S = np.eye( NLS ) # Initialize overlap matrix
    F = np.zeros( ( NLS, NLS ) ) # Initialize Fock matrix
    A = np.zeros( ( NEa, NEa ) ) # Initialize alpha determinant
    B = np.zeros( ( NEb, NEb ) ) # Initialize beta determinant
    for i in range( NLS ):
        for j in range( i+1, NLS ):
            # Overlap matrix:
            # 1) Alpha part:
            for i1 in range( NEa ): # Running over all alpha electrons
                ix1 = luoIXa[ i, i1 ] # Index of LUO-1
                for i2 in range( NEa ): # Running over all alpha electrons
                    ix2 = luoIXa[ j, i2 ] # Index of LUO-2
                    A[ i1, i2 ] = S_LUO[ ix1, ix2 ]
            detA = np.linalg.det( A )
            # 2) Beta part:
            for i1 in range( NEb ): # Running over all beta electrons
                ix1 = luoIXb[ i, i1 ] # Index of LUO-1
                for i2 in range( NEb ): # Running over all beta electrons
                    ix2 = luoIXb[ j, i2 ] # Index of LUO-2
                    B[ i1, i2 ] = S_LUO[ ix1, ix2 ]
            detB = np.linalg.det( B )
            S[ i, j ] = detA * detB
            S[ j, i ] = S[ i, j ]

            # Fock matrix:
            # 1) Alpha part:
            for i1 in range( NEa ): # Running over all rows of AF
                AF = A.copy() # Initialize AF as A
                ix1 = luoIXa[ i, i1 ] # Index of LUO-1
                # Change in the current row the elements to Fock integrals:
                for i2 in range( NEa ): # Running over all alpha electrons
                    ix2 = luoIXa[ j, i2 ] # Index of LUO-2
                    AF[ i1, i2 ] = F_LUO[ ix1, ix2 ]
                detAF = np.linalg.det( AF )
                F[ i, j ] += detAF * detB
            # 2) Beta part:
            for i1 in range( NEb ): # Running over all rows of BF
                BF = B.copy() # Initialize BF as B
                ix1 = luoIXb[ i, i1 ] # Index of LUO-1
                # Change in the current row the elements to Fock integrals:
                for i2 in range( NEb ): # Running over all beta electrons
                    ix2 = luoIXb[ j, i2 ] # Index of LUO-2
                    BF[ i1, i2 ] = F_LUO[ ix1, ix2 ]
                detBF = np.linalg.det( BF )
                F[ i, j ] += detA * detBF
            F[ j, i ] = F[ i, j ]

        # Speical case when i == j for Fock matrix (i.e., diagonal elements):
        # 1) Alpha part:
        for i1 in range( NEa ): # Running over all alpha electrons
            ix1 = luoIXa[ i, i1 ] # Index of LUO-1
            F[ i, i ] += F_LUO[ ix1, ix1 ]
        # 2) Beta part:
        for i1 in range( NEb ): # Running over all beta electrons
            ix1 = luoIXb[ i, i1 ] # Index of LUO-1
            F[ i, i ] += F_LUO[ ix1, ix1 ]

    # Determine the rank of S:
    _, E, _ = svd( S )
    rank = len( np.where( E > 1E-6 )[0] ) 
    return ( S, F, rank )
# enddef overlapFockLewis_open()



# Determine the indices of LUOs of a given Lewis structure:
#   lp0:   Indices of atoms in the LPs for a given Lewis structure
#   bd0:   Indices of atoms in the BDs for a given Lewis structure
def luoIx_one_Lewis( atIx_LUO, lp0, bd0 ):
    # Identify BD in terms of LUOs:
    nBD = len( bd0 ) # Number of BDs in this Lewis structure
    luoIx_BD = []
    for i in range( nBD ):
        ix_BD = bd0[i] - 1 # NOTE: Indices in BDO[] start from 1 (not 0)
        ix_luo = np.nonzero( np.all( atIx_LUO-ix_BD==0, axis=1 ) )[0]
        # Since ix_luo is a one-element array, put [0] to extract this number
        luoIx_BD.append( ix_luo[0] )
    #print( luoIx_BD )

    # Identify LP in terms of LUOs:
    nLP = len( lp0 ) # Number of LPs in this Lewis structure
    luoIx_LP = []
    for i in range( nLP ):
        ix_LP = [ lp0[i]-1, -1 ] # NOTE: Indices in LPO[] start from 1 (not 0)
        ix_luo = np.nonzero( np.all( atIx_LUO-ix_LP==0, axis=1 ) )[0]
        # Since ix_luo is a one-element array, put [0] to extract this number
        luoIx_LP.append( ix_luo[0] )
    #print( luoIx_LP )

    # Final indices of LUOs for all BDs, LPs in this Lewis structure:
    if nBD == 0:   # Only LPs 
        luoIx = np.array( luoIx_LP, dtype='i')
    elif nLP == 0: # Only BDs
        luoIx = np.array( luoIx_BD, dtype='i')
    else:          # Both BDs and LPs are present
        luoIx = np.hstack( ( np.array( luoIx_BD, dtype='i' ), \
                             np.array( luoIx_LP ,dtype='i' ) ) )

    return luoIx
# enddef luoIx_one_Lewis()


# Open-shell version of luoIx_one_Lewis()
# Determine the indices of LUOs of the same spin for a given Lewis structure:
#   lp0:   Indices of atoms in the LPs for a given Lewis structure
#   bd0:   Indices of atoms in the BDs for a given Lewis structure
#   r0:   Indices of atoms in the RAs/RBs for a given Lewis structure
def luoIx_one_Lewis_spin( atIx_LUO, lp0, bd0, r0 ):
    # Identify BD in terms of LUOs:
    nBD = len( bd0 ) # Number of BDs in this Lewis structure
    luoIx_BD = []
    for i in range( nBD ):
        ix_BD = bd0[i] - 1 # NOTE: Indices in BDO[] start from 1 (not 0)
        ix_luo = np.nonzero( np.all( atIx_LUO-ix_BD==0, axis=1 ) )[0]
        # Since ix_luo is a one-element array, put [0] to extract this number
        luoIx_BD.append( ix_luo[0] )
    #print( luoIx_BD )

    # Identify LP in terms of LUOs:
    nLP = len( lp0 ) # Number of LPs in this Lewis structure
    luoIx_LP = []
    for i in range( nLP ):
        ix_LP = [ lp0[i]-1, -1 ] # NOTE: Indices in LPO[] start from 1 (not 0)
        ix_luo = np.nonzero( np.all( atIx_LUO-ix_LP==0, axis=1 ) )[0]
        # Since ix_luo is a one-element array, put [0] to extract this number
        luoIx_LP.append( ix_luo[0] )
    #print( luoIx_LP )

    # Identify RA/RB (i.e., radical site) in terms of LUOs:
    nR = len( r0 ) # Number of radical sites in this Lewis structure
    luoIx_R = []
    for i in range( nR ):
        ix_R = [ r0[i]-1, -1 ] # NOTE: Indices in RAO[] start from 1 (not 0)
        ix_luo = np.nonzero( np.all( atIx_LUO-ix_R==0, axis=1 ) )[0]
        # Since ix_luo is a one-element array, put [0] to extract this number
        luoIx_R.append( ix_luo[0] )
    #print( luoIx_R )

    # Final indices of LUOs for all BDs, LPs, RAs/RBs in this Lewis structure:
    if nBD == 0 and nR == 0: # Only LPs 
        luoIx = np.array( luoIx_LP, dtype='i')
    elif nLP == 0 and nR ==0: # Only BDs
        luoIx = np.array( luoIx_BD, dtype='i')
    elif nBD == 0 and nLP ==0: # Only RAs/RBs
        luoIx = np.array( luoIx_R, dtype='i')
    elif nR == 0: # BDs and LPs 
        luoIx = np.hstack( ( np.array( luoIx_BD, dtype='i' ), \
                             np.array( luoIx_LP, dtype='i' ) ) )
    elif nLP == 0: # BDs and RAs/RBs 
        luoIx = np.hstack( ( np.array( luoIx_BD, dtype='i' ), \
                             np.array( luoIx_R, dtype='i' ) ) )
    elif nBD == 0: # LPs and RAs/RBs 
        luoIx = np.hstack( ( np.array( luoIx_LP, dtype='i' ), \
                             np.array( luoIx_R, dtype='i' ) ) )
    else: # BDs, LPs and RAs/RBs are all present
        luoIx = np.hstack( ( np.array( luoIx_BD, dtype='i' ), \
                             np.array( luoIx_LP, dtype='i' ), \
                             np.array( luoIx_R, dtype='i' ) ) )

    return luoIx
# enddef luoIx_one_Lewis_spin()

# Get matrix storing the indices of LUOs of a given set of Lewis structures:
#   atIx_LUO: Natural inidices of atoms in each LUO, being 0, 1, 2, ...
#   NE: Number of electrons
#   LP0: List of LPs for all considered Lewis structures
#   BD0: List of BDs for all considered Lewis structures
# Return:
#   luoIX: a matrix, each row corresponding to a Lewis structure (LS);
#          and the columns are the inidices of LUOs consitituting the LS
def luoIxLewis( atIx_LUO, NE, LP0, BD0 ):
    NLS = len( LP0 ) # Number of Lewis structures
    if NE % 2 == 1:
        raise NotImplementedError( 'Number of electrons must be even' )
    NP = NE // 2 # Number of electron pairs

    # Initialize the lists of LUO indices for all considered Lewis structures:
    luoIX = np.ones( ( NLS, NP ), dtype='i' ) * -1 
    for i in range( NLS ):
        luoIX[ i, : ] = luoIx_one_Lewis( atIx_LUO, LP0[i], BD0[i] )

    return luoIX
# enddef luoIxLewis()


# Open-shell version of luoIxLewis()
# Get matrix storing the indices of LUOs of a given set of Lewis structures:
#   atIx_LUO: Natural inidices of atoms in each LUO, being 0, 1, 2, ...
#   NE: Number of electrons
#   LP0: List of LPs for all considered Lewis structures
#   BD0: List of BDs for all considered Lewis structures
#   RA0: List of RAs for all considered Lewis structures
#   RB0: List of RBs for all considered Lewis structures
# Return:
#   luoIXa, luoIXb: Two matrices for alpha and beta electrons, respectively
#                   In each matrix, each row corresponding to a Lewis 
#                   structure (LS); and the columns are the inidices of LUOs 
#                   consitituting the LS
def luoIxLewis_open( atIx_LUO, NEa, NEb, LP0, BD0, RA0, RB0 ):
    NLS = len( LP0 ) # Number of Lewis structures

    # Initialize the lists of LUO indices for all considered Lewis structures:
    luoIXa = np.ones( ( NLS, NEa ), dtype='i' ) * -1 
    luoIXb = np.ones( ( NLS, NEb ), dtype='i' ) * -1 
    for i in range( NLS ):
        luoIXa[ i, : ] = luoIx_one_Lewis_spin( \
                atIx_LUO, LP0[i], BD0[i], RA0[i] )
        luoIXb[ i, : ] = luoIx_one_Lewis_spin( \
                atIx_LUO, LP0[i], BD0[i], RB0[i] )

    return ( luoIXa, luoIXb )
# enddef luoIxLewis_open()


# Determine which Lewis structures violate Rumer's rule
def antiRumerLewis( luoIX, atIx_LUO, NAt ):
    NLS = len( luoIX ) # Number of Lewis structures

    ix = []
    for i in range( NLS ):
        if not ifObeyRumer( luoIX[i], atIx_LUO, NAt ):
            ix.append( i )
    return np.array( ix, dtype='i' )
# enddef antiRumerLewis()

# Check if a given Lewis structure obeys Rumer's rule:
def ifObeyRumer( luo_ix, atIx_LUO, NAt ):
    N_LUO = len( luo_ix )
    # Get all bonds:
    bnd = []
    for i in range( N_LUO ):
        EP = atIx_LUO[ luo_ix[i], : ]
        if EP[1] != -1: # a BD
            bnd.append( EP )

    # Find if two bonds cross:
    Nbnd = len( bnd )
    for j in range( Nbnd ):
        for k in range( j+1, Nbnd ):
            if bondsCross( NAt, bnd[j], bnd[k] ):
                return False

    return True
# enddef ifObeyRumer()

# Determine if two bonds cross:
def bondsCross( N, bnd1, bnd2 ):
    # Generate coordinates of N-points on a circle:
    t = np.arange( N ) * 2*np.pi / N
    v = np.vstack( ( np.cos(t), np.sin(t) ) )

    x1 = v[ 0, bnd1[0] ]
    y1 = v[ 1, bnd1[0] ]
    x2 = v[ 0, bnd1[1] ]
    y2 = v[ 1, bnd1[1] ]

    u1 = v[ 0, bnd2[0] ]
    v1 = v[ 1, bnd2[0] ]
    u2 = v[ 0, bnd2[1] ]
    v2 = v[ 1, bnd2[1] ]

    #***** Cross point:  *****#
    # Solving equations:
    # (x2-x1)*a -(u2-u1)*b = u1-x1      (1)
    # (y2-y1)*a -(v2-v1)*b = v1-y1      (2)
    # Using Cramer's rule:
    D = np.linalg.det( np.array([ [x2-x1, u1-u2], [y2-y1, v1-v2] ]) )
    # Special case - two bonds are parallel:
    if abs(D) < 1E-10:
        return False
    a = np.linalg.det( np.array([ [u1-x1, u1-u2], [v1-y1, v1-v2] ]) ) / D
    b = np.linalg.det( np.array([ [x2-x1, u1-x1], [y2-y1, v1-y1] ]) ) / D
    # cross point:
    x = x1 + (x2-x1)*a
    y = y1 + (y2-y1)*a

    # Determine if the cross point is inside the unit circle:
    return x**2 + y**2 < 1
# enddef bondsCross()


# Wave function based resonance theory (WFRT) analysis:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   CLMO:        LMO coefficients in NAO basis
#   ELMO:        LMO energies
#   atIx1:         Indices of atoms in the resonance subsystem, start from 1
#   moIx1:         Indices of LMOs in the resonance subsystem, start from 1
#   MAX_LP:       Maximum number of lone pairs to be considered
#   projCut:      Cutoff value of projections relative to the principal LS
#   raoType:      'uni'--> universal RAOs;  'var'--> varying RAOs
#   inpBaseName:      Write RAOs into a *.fchk file; Default: '' --> Do not write
def wfrt( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
          MAX_LP=inf, projCut=0., raoType='uni', inpBaseName='', \
          raoFlipAtoms=None ):
    # Vectorize atIx and make moIx strat from ZERO (not 1):
    atIx = np.array( atIx1 ) # Starting from 1
    moIx = np.array( moIx1 ) - 1 # Starting from 0

    # Number of atoms:
    NAt = len( atIx )
    # Number of electrons:
    NE = len( moIx ) * 2  # Closed-shell system

    # Choose all LMOs involved in the resoncance subsystem:
    C_LMO = CLMO[:,moIx]

    if type( naoInfo ) is list: # Interface for wfrt_hmo()
        ifHMO = True
        aoIx = naoInfo
    else:                       # Normal usage
        ifHMO = False
        aoIx = naoInfo.naoIx

    NP = NE // 2 # Number of electron pairs

    # Get all resonance structures:
    #  LP0:   List of LPs for all considered Lewis structures
    #  BD0:   List of BDs for all considered Lewis structures
    #  Suffix 0 means the indices of atoms stored in LP0[] and BD0[] follow the
    #  order of appearance, i.e., 1, 2, 3, ..., which is NOT necessarily the
    #  actural indices of atoms appeared in the given molecule
    #  The latter indices are stored in LP[] and BD[]
    LP0, BD0 = lewis_read( NAt, NE, MAX_LP )

    # Set actual indices of atoms:
    atIx = np.array( atIx )
    LP = list( map(lambda x: atIx[x-1], LP0 ) )
    BD = []
    for arr in BD0:
        BD.append(  np.array( list( map(lambda x: atIx[x-1], arr) ) )  )

    NLS = len( LP )  # Total number of Lewis structures
    print( '%i Lewis structures in total to be considered' % NLS )

    # Density matrix from pi-LMOs in NAO basis:
    D0 = densityMatrixLMO( C_LMO )
    # Total energy of the LMOs for the resonance subsystem:
    Etot_LMO = LMOEnergy( ELMO, moIx )
    # Fock matrix in the original AO basis:
    F = FockMatrix( FchkInfo.C, FchkInfo.E )
    # Fock matrix in NAO basis:
    FNAO = CAONAO.T @ F @ CAONAO

    # Resonance atomic orbitals (RAOs) in AO basis (NAOs, Lowdin's AOs etc.):
    if raoType == 'uni': # Universal RAOs that are independent of LSs
        CRAO = rao_uni( D0, aoIx, raoFlipAtoms )
    elif raoType == 'uni-gen': # General universal RAOs
        raise NotImplementedError( 'General universal RAOs not supported yet '
                'in WFRT' )
    else:
        raise ValueError( 
                'Unrecognized value for argument raoType in function wfrt()' )

    # Write RAOs into a *.fchk file:
    if len(inpBaseName) > 0:
        raoOutpName = inpBaseName + '_RAO.fchk'
        CAORAO = CAONAO @ CRAO
        FRAO = CAORAO.T @ F @ CAORAO
        E_RAO = np.diag( FRAO )
        writeFchkOrb( inpBaseName+'.fchk', raoOutpName, CAORAO, E_RAO )
        print( 'RAOs written into ' + raoOutpName )


    # Lewis unit orbitals (LUOs) in AO basis (NAOs, Lowdin's AOs etc.):
    atIx_LUO, C_LUO, S_LUO, F_LUO = luo( CRAO, FNAO, atIx, NAt, NE )
    N_LUO = len( S_LUO )
    # Print the RAO indices in each LUO:
    print( '-'*25 )
    print( '%5s   %4s   %s' % ('LUO', 'Type', 'RAO(s)') )
    print( '-'*25 )
    k = 1
    for ix in atIx_LUO:
        print( '%5i' % k, end='' )
        if ix[1] == -1: # LP
            print( '  (LP):  %3i' % atIx[ix[0]] )
        else:           # BD
            print( '  (BD):  %3i %3i' % (atIx[ix[0]],atIx[ix[1]]) )
        k += 1
    print( '-'*25 )

    # Matrix storing the indices of LUOs of a given set of Lewis structures:
    luoIX = luoIxLewis( atIx_LUO, NE, LP0, BD0 )


    # Overlap matrix between LUOs and LMOs:
    S_LUO_LMO = np.zeros( ( N_LUO, NP ) )
    for i in range( N_LUO ): # Running over all LUOs
        for j in range( NP ): # Running over all LMOs
            S_LUO_LMO[ i, j ] = C_LUO[ :, i ].T @ C_LMO[ :, j ]
    #writeMat( S_LUO_LMO )


    # Projection of Lewis structure determinants onto that of Psi0:
    P = np.zeros( ( NLS, 1 ) ) # Initialize the projections
    A = np.zeros( ( NP, NP ) ) # Initialize single-spin determinant
    for i in range( NLS ):
        for j in range( NP ): # Running over all electron pairs in an LS
            iLUO = luoIX[ i, j ] # Index of LUO for this electron pair
            for iLMO in range( NP ): # Index of LMO
                A[ j, iLMO ] = S_LUO_LMO[ iLUO, iLMO ]
        P[i] = np.linalg.det( A ) ** 2
    #writeMat( P, '/dev/stdout', '%12.6f' )


    # Sort by P:
    ix_sort = np.argsort( -P[:,0] )
    P = P[ ix_sort ]

    # Update LP[], BD[], and luoIX:
    LP = [ LP[i] for i in ix_sort ]
    BD = [ BD[i] for i in ix_sort ]
    luoIX = luoIX[ ix_sort, : ]

    # Strings for Lewis structures:
    str_LS = []
    for i in range( NLS ):
        str_LS.append( lewis_str( LP[i], BD[i] ) )

    #------------------------------------------------------------------
    # Print out the reordered list of Lewis structures:
    #------------------------------------------------------------------
    print( '-'*80 )
    print( '%5s %15s    ' % ( 'No.', 'Projection'), end='' )
    print( ' '*NE, end='' )
    print( ' LUOs', end='' )
    print( ' '*NE, end='' )
    print( '  Lewis structure' )
    print( '-'*80 )
    for i in range( NLS ):
        print( '%5i %15.10f   ' % ( i+1, P[i] ), end='' )
        for ix in luoIX[i]:
            print( ' %4i' % ( ix+1 ), end='' )
        print( '     %s' % str_LS[i] )
    print( '-'*80 )
    #------------------------------------------------------------------

    ## Overlap matrix between determinants of all considered Lewis structures:
    #S, rank = overlapLewis( NE, luoIX, S_LUO )
    ##writeMat( S, '/dev/stdout', ' %4.2f' )
    #print( '%i linearly independent Lewis structures' % rank )
    
    # Apply Rumer's rule to pick up linearly independent Lewis structures:
    ix_depLS = antiRumerLewis( luoIX, atIx_LUO, NAt ) # Indices starting from 0
    ix_indepLS = np.setdiff1d( range(NLS), ix_depLS ) # Indices starting from 0
    N_depLS = len( ix_depLS )
    print( '%i Lewis structures violating Rumer\'s rule:' % N_depLS )
    for i in range( N_depLS ):
        print( ' %5i' % (ix_depLS[i]+1), end='' )
        if (i+1) % 10 == 0:
            print()
    if (i+1) % 10 != 0:
        print()


    #------------------------------------------------------------------
    # Solve linear equation for coeff. using ridge linear regression:
    #------------------------------------------------------------------
    print( 'Solving linear equation for WFRT coefficients ...' )

    # Exclude certain Lewis structures:
    #   1. Remove the anti-Rumer Lewis structures:
    ix_reducLS = ix_indepLS
    print( '%i Lewis structures violating Rumer\'s rule excluded' 
            % ( len( ix_depLS ) ) )
    #   2. Remove Lewis structures having small projections, if applicable:
    ix_minorLS = np.nonzero( P < projCut )[0]
    print( '%i Lewis structures of small projections ( < %.2E ) excluded' 
            % ( len( ix_minorLS), projCut ) )
    ix_reducLS = np.setdiff1d( ix_reducLS, ix_minorLS ) # Ix. starting from 0
    print( 'Using totally %i Lewis structures to expand the actual wave '
            'function' % ( len( ix_reducLS ) ) )

    # Reduce P[] and luoIX[]:
    Pr = P[ ix_reducLS ]
    luoIXr = luoIX[ ix_reducLS, : ]
    #Sr = S[np.ix_( ix_reducLS, ix_reducLS )]
    # Overlap matrix between determinants of all considered Lewis structures:
    #Sr, rank = overlapLewis( NE, luoIXr, S_LUO )
    Sr, Fr, rank = overlapFockLewis( NE, luoIXr, S_LUO, F_LUO )
    #writeMat( Sr, '/dev/stdout', ' %4.2f' )
    #print('*'*80)
    #writeMat( Fr, '/dev/stdout', ' %5.2f' )
    #print( '%i linearly independent Lewis structures' % rank )
    assert( len( ix_reducLS ) == rank )
    

    # Use ridge linear regression:
    a = 1E-10
    Cr = np.linalg.pinv( Sr@Sr + a*np.eye(len(Pr)) ) @ Sr @ Pr

    # Normalization:
    NormFactor = sqrt( Cr.T @ Sr @ Cr )  # <Psi|Psi>
    print( 'Normalizing the wave function by a factor of '
            '%.10f ...' % (NormFactor) )
    Cr /= NormFactor
    assert( abs( Cr.T @ Sr @ Cr - 1. ) < 1E-10 )


    # Reproducibility:
    R = Pr.T @ Cr


    # Weights:
    Wr_Mull = weightMulliken( Cr, Sr )
    Wr_Bick = weightBickelhaupt( Cr, Sr )
    Wr_RS = weightRos_Schuit( Cr )
    Wr_Lowd = weightLowdin( Cr, Sr )
    Wr_P2WSO = weightP2WSO( Cr, Sr, Pr )
    #writeMat( Wr_Lowd, '/dev/stdout', '%12.6f' )

    # Recover the full coefficient vector by setting zeros for dependent LSs:
    C = np.zeros( NLS ) # Initialize
    W_Mull = np.zeros( NLS ) # Initialize
    W_Bick = np.zeros( NLS ) # Initialize
    W_RS = np.zeros( NLS ) # Initialize
    W_Lowd = np.zeros( NLS ) # Initialize
    W_P2WSO = np.zeros( NLS ) # Initialize
    for i in range( NLS ):
        ix_reduc = np.nonzero( ix_reducLS == i )[0]
        if len( ix_reduc ) > 0:
            C[i] = Cr[ ix_reduc ]
            W_Mull[i] = Wr_Mull[ ix_reduc ]
            W_Bick[i] = Wr_Bick[ ix_reduc ]
            W_RS[i] = Wr_RS[ ix_reduc ]
            W_Lowd[i] = Wr_Lowd[ ix_reduc ]
            W_P2WSO[i] = Wr_P2WSO[ ix_reduc ]
    #assert( abs( C.T @ S @ C - 1. ) < 1E-10 )


    # Energies of non-interacting Lewis structures:
    N_reducLS = len( ix_reducLS )
    Er = np.diag( Fr ).reshape( N_reducLS, 1 )
    Erel_r = Er - Etot_LMO
    Er_min = np.min( Er )
    REr = Er - Er_min # Relative energies of LSs with respect to Er_min
    # The Er_min gives the Pauling--Wheland resonance energy:
    if ifHMO:
        Ereson = Er_min - Etot_LMO
    else:
        Ereson = au2kcal( Er_min - Etot_LMO )
    ix_min_LS = []
    for k in range( N_reducLS ):
        if ( not ifHMO and abs( Er[k] - Er_min ) < 0.001/au2kcal(1.) ) \
            or ( ifHMO and abs( Er[k] - Er_min ) < 0.001 ):
            ix_min_LS.append( k+1 )

    # Total energy of the expanded wave function:
    Esum = 0.
    for i in range( N_reducLS ):
        Esum += Cr[i] * Cr[i] * Fr[i,i]
        for j in range( i+1, N_reducLS ):
            if j == i:
                continue
            Esum += 2 * Cr[i] * Cr[j] * Fr[i,j]

    # Effective resonance energies of LSs using Mulliken-like partition scheme:
    REeff_r = Cr * Cr * REr
    for i in range( N_reducLS ):
        for j in range( N_reducLS ):
            if j == i:
                continue
            REeff_r[i] += Cr[i] * Cr[j] * ( Fr[i,j] - Er_min * Sr[i,j] )

    # Recover the full energy information by setting zeros for dependent LSs:
    Erel = np.zeros( NLS ) # Initialize
    RE = np.zeros( NLS ) # Initialize
    REeff = np.zeros( NLS ) # Initialize
    for i in range( NLS ):
        ix_reduc = np.nonzero( ix_reducLS == i )[0]
        if len( ix_reduc ) > 0:
            Erel[i] = Erel_r[ ix_reduc ]
            RE[i] = REr[ ix_reduc ]
            REeff[i] = REeff_r[ ix_reduc ]
    #------------------------------------------------------------------


    #------------------------------------------------------------------
    # Print out the final results:
    #------------------------------------------------------------------
    ## Rearrange in descending order of P2WSO weights:
    ix_sort = np.argsort( -W_P2WSO )
    # Rearrange in descending order of Lowdin weights:
    #ix_sort = np.argsort( -W_Lowd )
    print( '-'*105 )
    print( '%5s %10s %11s' % ( 'No.', 'Projection', 'Coefficient' ), end='' )
    print( ' %7s' % 'RE', end='' )
    print( ' %6s' % 'REeff', end='' )
    print( ' %7s' % 'Mulli.', end='' )
    print( ' %7s' % 'Bickel.', end='' )
    print( ' %7s' % 'Ros-Sc.', end='' )
    print( ' %7s' % 'Lowdin', end='' )
    print( ' %7s' % 'P2WSO', end='' )
    print( '  Lewis structure' )
    print( '-'*105 )
    for iLS in range( NLS ):
        i = ix_sort[ iLS ]
        # Only show Lewis structures in the reduced set:
        if np.count_nonzero( ix_reducLS == i ) == 0:
            continue
        print( '%5i %10.6f %11.7f' % ( i+1, P[i], C[i] ), end='' )
        if ifHMO:
            print( ' %7.3f' % ( Erel[i] ), end='' )
            print( ' %7.3f' % ( REeff[i] ), end='' )
        else:
            print( ' %7.2f' % ( au2kcal( Erel[i] ) ), end='' )
            print( ' %7.2f' % ( au2kcal( REeff[i] ) ), end='' )
        print( ' %6.2f%%' % ( W_Mull[i]*100 ), end='' )
        print( ' %6.2f%%' % ( W_Bick[i]*100 ), end='' )
        print( ' %6.2f%%' % ( W_RS[i]*100 ), end='' )
        print( ' %6.2f%%' % ( W_Lowd[i]*100 ), end='' )
        print( ' %6.2f%%' % ( W_P2WSO[i]*100 ), end='' )
        print( '  %s' % str_LS[i] )
    print( '-'*105 )
    print( 'Reproducibility = %10.6f%%' % ( R*100 ) )

    if ifHMO:
        print( 'Total HF energy of reference structure is %.4f |beta|' 
                % Etot_LMO )
        print( 'Total energy of expanded wave function is %.4f |beta|' 
            % Esum )
    else:
        print( 'Total HF energy of reference structure is %.8E a.u.' 
                % Etot_LMO )
        print( 'Total energy of expanded wave function is %.8E a.u.' 
            % Esum )
    print( '  Percent error = %.6f%% ' 
            % ( ( Esum/Etot_LMO - 1 ) * 100 * np.sign( Etot_LMO ) ) )
    if ifHMO:
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.4f |beta|' % Ereson )
        Ereson_expan =  -REeff.sum()
        print( 'Pauling--Wheland resonance energy given by expanded wave '
               'function is %.4f |beta|' % Ereson_expan )
    else:
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.2f kcal/mol' % Ereson )
        Ereson_expan =  -au2kcal( REeff.sum() )
        print( 'Pauling--Wheland resonance energy given by expanded wave '
               'function is %.2f kcal/mol' % Ereson_expan )
    print( '  Percent error = %.6f%% ' % ( ( 1 - Ereson_expan/Ereson )*100 ) )
    if len( ix_min_LS ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_LS:
        print( ' #%i' % i, end='' )
    print()

    #------------------------------------------------------------------

    return
# enddef wfrt()


# WFRT analysis with a specified set of Lewis structures:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   CLMO:        LMO coefficients in NAO basis
#   ELMO:        LMO energies
#   atIx1:         Indices of atoms in the resonance subsystem, start from 1
#   moIx1:         Indices of LMOs in the resonance subsystem, start from 1
#   ixLS:         Lewis structures specified in either of the following ways:
#                 (1) List of indices of the specified Lewis structures
#                 This choice involves enumeration of all possbile LSs. If the
#                 enumeration has been done previously, file LEWIS_Nc_Ne.dat
#                 will be read in.
#                 NOTE: 1) Indices start from 1 (not 0)
#                       2) The list is from ordering of projections, determined
#                          by calling the wfrt() function prior to calling this
#                          function
#                 (2) Straightfoward specification given by LP[] and BD[],
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[] and BD[], atomic indices start from 1 
#                       (not 0), and they correspond to ACTUAL indices of 
#                       atoms, which may not be consecutive, not as those in 
#                       LP0[] and BD0[] (also see notes in body of this func)
#   sciform:      whether showing numerical results in scientific format
#   raoType:      'uni'--> universal RAOs;  'var'--> varying RAOs
#   inpBaseName:      Write RAOs into a *.fchk file; Default: '' --> Do not write
def wfrt_spec( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, ixLS, \
               sciform=False, raoType='uni', inpBaseName='', \
               raoFlipAtoms=None ):
    # Vectorize atIx and make moIx strat from ZERO (not 1):
    atIx = np.array( atIx1 ) # Starting from 1
    moIx = np.array( moIx1 ) - 1 # Starting from 0

    # Number of atoms:
    NAt = len( atIx )
    # Number of electrons:
    NE = len( moIx ) * 2  # Closed-shell system

    # Choose all LMOs involved in the resoncance subsystem:
    C_LMO = CLMO[:,moIx]

    if type( naoInfo ) is list: # Interface for wfrt_hmo()
        ifHMO = True
        aoIx = naoInfo
    else:                       # Normal usage
        ifHMO = False
        aoIx = naoInfo.naoIx

    NP = NE // 2 # Number of electron pairs

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
        atIx = np.array( atIx, dtype='i' )
        LP = list( map(lambda x: atIx[x-1], LP0 ) )
        BD = []
        for arr in BD0:
            BD.append(  np.array( list( map(lambda x: atIx[x-1], arr) ), \
                                  dtype='i' )  )

    #-- OPTION 2: Straightforward specification ----------
    # CLOSED-SHELL CASE:
    elif len( ixLS ) == 2:
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

    # OPEN-SHELL CASE:
    else:
        ifEnume = False
        LP, BD, RA, RB, isRadical = ixLS

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

        RA0 = list( map(lambda x: sqIx[x-1]+1, RA ) )
        RB0 = list( map(lambda x: sqIx[x-1]+1, RB ) )
    #------------------------------------------------------------------------

    NLS = len( LP )  # Total number of Lewis structures
    print( '%i Lewis structures in total to be considered' % NLS )

    # Density matrix from pi-LMOs in NAO basis:
    D0 = densityMatrixLMO( C_LMO )
    # Total energy of the LMOs for the resonance subsystem:
    Etot_LMO = LMOEnergy( ELMO, moIx )
    # Fock matrix in the original AO basis:
    F = FockMatrix( FchkInfo.C, FchkInfo.E )
    # Fock matrix in NAO basis:
    FNAO = CAONAO.T @ F @ CAONAO

    # Resonance atomic orbitals (RAOs) in AO basis (NAOs, Lowdin's AOs etc.):
    if raoType == 'uni': # Universal RAOs that are independent of LSs
        if ifHMO:
            CRAO = np.identity( NAt )
        else:
            CRAO = rao_uni( D0, aoIx, raoFlipAtoms )
    elif raoType == 'uni-gen': # General universal RAOs
        raise NotImplementedError( 'General universal RAOs not supported yet '
                'in WFRT' )
    else:
        raise ValueError( 
                'Unrecognized value for argument raoType in function wfrt()' )

    # Write RAOs into a *.fchk file:
    if len(inpBaseName) > 0:
        raoOutpName = inpBaseName + '_RAO.fchk'
        CAORAO = CAONAO @ CRAO
        FRAO = CAORAO.T @ F @ CAORAO
        E_RAO = np.diag( FRAO )
        writeFchkOrb( inpBaseName+'.fchk', raoOutpName, CAORAO, E_RAO )
        print( 'RAOs written into ' + raoOutpName )

    # Lewis unit orbitals (LUOs) in AO basis (NAOs, Lowdin's AOs etc.):
    atIx_LUO, C_LUO, S_LUO, F_LUO = luo( CRAO, FNAO, atIx, NAt, NE )
    N_LUO = len( S_LUO )
    # Print the RAO indices in each LUO:
    print( '-'*25 )
    print( '%5s   %4s   %s' % ('LUO', 'Type', 'RAO(s)') )
    print( '-'*25 )
    k = 1
    for ix in atIx_LUO:
        print( '%5i' % k, end='' )
        if ix[1] == -1: # LP
            print( '  (LP):  %3i' % atIx[ix[0]] )
        else:           # BD
            print( '  (BD):  %3i %3i' % (atIx[ix[0]],atIx[ix[1]]) )
        k += 1
    print( '-'*25 )

    # Matrix storing the indices of LUOs of a given set of Lewis structures:
    luoIX = luoIxLewis( atIx_LUO, NE, LP0, BD0 )

    # Overlap matrix between LUOs and LMOs:
    S_LUO_LMO = np.zeros( ( N_LUO, NP ) )
    for i in range( N_LUO ): # Running over all LUOs
        for j in range( NP ): # Running over all LMOs
            S_LUO_LMO[ i, j ] = C_LUO[ :, i ].T @ C_LMO[ :, j ]
    #writeMat( S_LUO_LMO )

    # Projection of Lewis structure determinants onto that of Psi0:
    P = np.zeros( ( NLS, 1 ) ) # Initialize the projections
    A = np.zeros( ( NP, NP ) ) # Initialize single-spin determinant
    for i in range( NLS ):
        for j in range( NP ): # Running over all electron pairs in an LS
            iLUO = luoIX[ i, j ] # Index of LUO for this electron pair
            for iLMO in range( NP ): # Index of LMO
                A[ j, iLMO ] = S_LUO_LMO[ iLUO, iLMO ]
        P[i] = np.linalg.det( A ) ** 2
    #writeMat( P, '/dev/stdout', '%12.6f' )

# DO NOT SORT FOR wfrt_spec:
#    # Sort by P:
#    ix_sort = np.argsort( -P[:,0] )
#    P = P[ ix_sort ]

    # Do not sort:
    ix_sort = np.arange( NLS )
    P = P[ ix_sort ]

    # Update LP[], BD[], and luoIX:
    LP = [ LP[i] for i in ix_sort ]
    BD = [ BD[i] for i in ix_sort ]
    luoIX = luoIX[ ix_sort, : ]

    # Strings for Lewis structures:
    str_LS = []
    for i in range( NLS ):
        str_LS.append( lewis_str( LP[i], BD[i] ) )

    #------------------------------------------------------------------
    # Print out the reordered list of Lewis structures:
    #------------------------------------------------------------------
    print( '-'*80 )
    print( '%5s %15s    ' % ( 'No.', 'Projection'), end='' )
    print( ' '*NE, end='' )
    print( ' LUOs', end='' )
    print( ' '*NE, end='' )
    print( '  Lewis structure' )
    print( '-'*80 )
    for i in range( NLS ):
        print( '%5i %15.10f   ' % ( i+1, P[i] ), end='' )
        for ix in luoIX[i]:
            print( ' %4i' % ( ix+1 ), end='' )
        print( '     %s' % str_LS[i] )
    print( '-'*80 )
    #------------------------------------------------------------------

    # NOTE: Since the user specified set is usually small and not complete,
    #       we DO NOT apply Rumer's rule to deal with linearly dependence

    #------------------------------------------------------------------
    # Solve linear equation for coeff. using ridge linear regression:
    #------------------------------------------------------------------
    print( 'Solving linear equation for WFRT coefficients ...' )

    # Only choose the specified Lewis structures:
    if ifEnume:
        ix_reducLS = np.array( ixLS, dtype='i' ) - 1
        print( 'Using the %i specified Lewis structures to expand the '
               'actual wave function' % ( len( ix_reducLS ) ) )
    else:
        ix_reducLS = np.arange( NLS )
        print( 'Using %i straightforwardly specified Lewis structures to '
               'expand the actual wave function' % ( len( ix_reducLS ) ) )

    # Reduce P[] and luoIX[]:
    Pr = P[ ix_reducLS ]
    luoIXr = luoIX[ ix_reducLS, : ]

    # Overlap matrix between determinants of all considered Lewis structures:
    #Sr, rank = overlapLewis( NE, luoIXr, S_LUO )
    Sr, Fr, rank = overlapFockLewis( NE, luoIXr, S_LUO, F_LUO )
    #writeMat( Sr, '/dev/stdout', ' %4.2f' )
    #print( '%i linearly independent Lewis structures' % rank )
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
    

    # Use ridge linear regression:
    a = 1E-10
    Cr = np.linalg.pinv( Sr@Sr + a*np.eye(len(Pr)) ) @ Sr @ Pr

    # Normalization:
    NormFactor = sqrt( Cr.T @ Sr @ Cr )  # <Psi|Psi>
    print( 'Normalizing the wave function by a factor of '
            '%.10f ...' % (NormFactor) )
    Cr /= NormFactor
    assert( abs( Cr.T @ Sr @ Cr - 1. ) < 1E-10 )


    # Reproducibility:
    R = Pr.T @ Cr

    # Weights:
    if rank_defi == 0:
        Wr_Mull = weightMulliken( Cr, Sr )
        Wr_Bick = weightBickelhaupt( Cr, Sr )
        Wr_RS = weightRos_Schuit( Cr )
        Wr_Lowd = weightLowdin( Cr, Sr )
        Wr_P2WSO = weightP2WSO( Cr, Sr, Pr )
        #writeMat( Wr_Lowd, '/dev/stdout', '%12.6f' )

    # Recover the full coefficient vector by setting zeros for unchosen LSs:
    C = np.zeros( NLS ) # Initialize
    if rank_defi == 0:
        W_Mull = np.zeros( NLS ) # Initialize
        W_Bick = np.zeros( NLS ) # Initialize
        W_RS = np.zeros( NLS ) # Initialize
        W_Lowd = np.zeros( NLS ) # Initialize
        W_P2WSO = np.zeros( NLS ) # Initialize
    for i in range( NLS ):
        ix_reduc = np.nonzero( ix_reducLS == i )[0]
        if len( ix_reduc ) > 0:
            C[i] = Cr[ ix_reduc ]
            if rank_defi == 0:
                W_Mull[i] = Wr_Mull[ ix_reduc ]
                W_Bick[i] = Wr_Bick[ ix_reduc ]
                W_RS[i] = Wr_RS[ ix_reduc ]
                W_Lowd[i] = Wr_Lowd[ ix_reduc ]
                W_P2WSO[i] = Wr_P2WSO[ ix_reduc ]
    #assert( abs( C.T @ S @ C - 1. ) < 1E-10 )


    # Energies of non-interacting Lewis structures:
    N_reducLS = len( ix_reducLS )
    Er = np.diag( Fr ).reshape( N_reducLS, 1 )
    Erel_r = Er - Etot_LMO
    Er_min = np.min( Er )
    REr = Er - Er_min # Relative energies of LSs with respect to Er_min
    # The Er_min gives the Pauling--Wheland resonance energy:
    if ifHMO:
        Ereson = Er_min - Etot_LMO
    else:
        Ereson = au2kcal( Er_min - Etot_LMO )
    ix_min_LS = []
    for k in range( N_reducLS ):
        if ( not ifHMO and abs( Er[k] - Er_min ) < 0.001/au2kcal(1.) ) \
            or ( ifHMO and abs( Er[k] - Er_min ) < 0.001 ):
            ix_min_LS.append( k+1 )

    # Total energy of the expanded wave function:
    Esum = 0.
    for i in range( N_reducLS ):
        Esum += Cr[i] * Cr[i] * Fr[i,i]
        for j in range( i+1, N_reducLS ):
            if j == i:
                continue
            Esum += 2 * Cr[i] * Cr[j] * Fr[i,j]

    # Effective resonance energies of LSs using Mulliken-like partition scheme:
    REeff_r = Cr * Cr * REr
    for i in range( N_reducLS ):
        for j in range( N_reducLS ):
            if j == i:
                continue
            REeff_r[i] += Cr[i] * Cr[j] * ( Fr[i,j] - Er_min * Sr[i,j] )

    # Recover the full energy information by setting zeros for dependent LSs:
    Erel = np.zeros( NLS ) # Initialize
    RE = np.zeros( NLS ) # Initialize
    REeff = np.zeros( NLS ) # Initialize
    for i in range( NLS ):
        ix_reduc = np.nonzero( ix_reducLS == i )[0]
        if len( ix_reduc ) > 0:
            Erel[i] = Erel_r[ ix_reduc ]
            RE[i] = REr[ ix_reduc ]
            REeff[i] = REeff_r[ ix_reduc ]
    #------------------------------------------------------------------


    #------------------------------------------------------------------
    # Print out the final results:
    #------------------------------------------------------------------
    ## Rearrange in descending order of P2WSO weights:
    # Rearrange in descending order of Lowdin weights:
    if rank_defi == 0:
        ix_sort = np.argsort( -W_P2WSO )
        #ix_sort = np.argsort( -W_Lowd )
    else:
        ix_sort = np.argsort( Erel )
    print( '-'*105 )
    if sciform:
        print( '%5s %15s %11s' % ('No.', 'Projection', 'Coefficient'), end='' )
    else:
        print( '%5s %10s %11s' % ('No.', 'Projection', 'Coefficient'), end='' )
    print( ' %7s' % 'Erel', end='' )
    print( ' %7s' % 'REeff', end='' )
    if rank_defi == 0:
        print( ' %7s' % 'Mulli.', end='' )
        print( ' %7s' % 'Bickel.', end='' )
        print( ' %7s' % 'Ros-Sc.', end='' )
        print( ' %7s' % 'Lowdin', end='' )
        if sciform:
            print( ' %15s' % 'P2WSO', end='' )
        else:
            print( ' %7s' % 'P2WSO', end='' )
    print( '  Lewis structure' )
    print( '-'*105 )
    for iLS in range( NLS ):
        i = ix_sort[ iLS ]
        # Only show Lewis structures in the reduced set:
        if np.count_nonzero( ix_reducLS == i ) == 0:
            continue
        if sciform:
            print( '%5i %15.8E %11.7f' % ( i+1, P[i], C[i] ), end='' )
        else:
            print( '%5i %10.6f %11.7f' % ( i+1, P[i], C[i] ), end='' )
        if ifHMO:
            print( ' %7.3f' % ( Erel[i] ), end='' )
            print( ' %7.3f' % ( REeff[i] ), end='' )
        else:
            print( ' %7.2f' % ( au2kcal( Erel[i] ) ), end='' )
            print( ' %7.2f' % ( au2kcal( REeff[i] ) ), end='' )
        if rank_defi == 0:
            print( ' %6.2f%%' % ( W_Mull[i]*100 ), end='' )
            print( ' %6.2f%%' % ( W_Bick[i]*100 ), end='' )
            print( ' %6.2f%%' % ( W_RS[i]*100 ), end='' )
            print( ' %6.2f%%' % ( W_Lowd[i]*100 ), end='' )
            if sciform:
                print( ' %15.8E%%' % ( W_P2WSO[i]*100 ), end='' )
            else:
                print( ' %6.2f%%' % ( W_P2WSO[i]*100 ), end='' )
        print( '  %s' % str_LS[i] )
    print( '-'*105 )
    print( 'Reproducibility = %10.6f%%' % ( R*100 ) )

    if ifHMO:
        print( 'Total HF energy of reference structure is %.4f |beta|' 
                % Etot_LMO )
        print( 'Total energy of expanded wave function is %.4f |beta|' 
            % Esum )
    else:
        print( 'Total HF energy of reference structure is %.8E a.u.' 
                % Etot_LMO )
        print( 'Total energy of expanded wave function is %.8E a.u.' 
            % Esum )
    print( '  Percent error = %.6f%% ' 
            % ( ( Esum/Etot_LMO - 1 ) * 100 * np.sign( Etot_LMO ) ) )
    if ifHMO:
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.4f |beta|' % Ereson )
        Ereson_expan =  -REeff.sum()
        print( 'Pauling--Wheland resonance energy given by expanded wave '
               'function is %.4f |beta|' % Ereson_expan )
    else:
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.2f kcal/mol' % Ereson )
        Ereson_expan =  -au2kcal( REeff.sum() )
        print( 'Pauling--Wheland resonance energy given by expanded wave '
               'function is %.2f kcal/mol' % Ereson_expan )
    print( '  Percent error = %.6f%% ' % ( ( 1 - Ereson_expan/Ereson )*100 ) )
    if len( ix_min_LS ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_LS:
        print( ' #%i' % i, end='' )
    print()
    #------------------------------------------------------------------


    # Write overlap matrix to external file:
    if globals.ifWrite_Overlap:
        # Write the unsorted OVERLAP:
        overlapFile = globals.basename + '.OVERLAP'
        print( 'Writing overlap matrix to file %s ...' % overlapFile )
        np.savetxt( overlapFile, Sr, fmt='%24.15E' )
        print( 'File %s written' % overlapFile )

        # Write the sorted OVERLAP:
        overlapFile = globals.basename + '.OVERLAP_SORT'
        print( 'Writing overlap matrix to file %s ...' % overlapFile )
        Sr_sortByW = np.zeros( ( NLS, NLS ) )
        for iLS in range( NLS ):
            i = ix_sort[ iLS ]
            for jLS in range( NLS ):
                j = ix_sort[ jLS ]
                Sr_sortByW[iLS][jLS] = Sr[i][j]
        np.savetxt( overlapFile, Sr_sortByW, fmt='%24.15E' )
        print( 'File %s written' % overlapFile )

    # Write Hamiltonian matrix to external file:
    if globals.ifWrite_Hamilt:
        # Write the unsorted HAMILTONIAN:
        HamiltonianFile = globals.basename + '.HAMILT'
        print( 'Writing Hamiltonian matrix to file %s ...' % HamiltonianFile )
        np.savetxt( HamiltonianFile, Fr, fmt='%24.15E' )
        print( 'File %s written' % HamiltonianFile )

        # Write the sorted HAMILTONIAN:
        HamiltonianFile = globals.basename + '.HAMILT_SORT'
        print( 'Writing Hamiltonian matrix to file %s ...' % HamiltonianFile )
        Fr_sortByW = np.zeros( ( NLS, NLS ) )
        for iLS in range( NLS ):
            i = ix_sort[ iLS ]
            for jLS in range( NLS ):
                j = ix_sort[ jLS ]
                Fr_sortByW[iLS][jLS] = Fr[i][j]
        np.savetxt( HamiltonianFile, Fr_sortByW, fmt='%24.15E' )
        print( 'File %s written' % HamiltonianFile )

    return
# enddef wfrt_spec()


# Spin-unrestricted version of wfrt_spec():
# WFRT analysis with a specified set of Lewis structures:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   CLMO:        LMO coefficients in NAO basis
#   ELMO:        LMO energies
#   atIx1:         Indices of atoms in the resonance subsystem, start from 1
#   moIx1:         Indices of LMOs in the resonance subsystem, start from 1
#   ixLS:         Lewis structures specified in either of the following ways:
#                 (1) List of indices of the specified Lewis structures
#                 This choice involves enumeration of all possbile LSs. If the
#                 enumeration has been done previously, file LEWIS_Nc_Ne.dat
#                 will be read in.
#                 NOTE: 1) Indices start from 1 (not 0)
#                       2) The list is from ordering of projections, determined
#                          by calling the wfrt() function prior to calling this
#                          function
#                 (2) Straightfoward specification given by LP[] and BD[],
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[] and BD[], atomic indices start from 1 
#                       (not 0), and they correspond to ACTUAL indices of 
#                       atoms, which may not be consecutive, not as those in 
#                       LP0[] and BD0[] (also see notes in body of this func)
#   sciform:      whether showing numerical results in scientific format
#   raoType:      'uni'--> universal RAOs;  'var'--> varying RAOs
#   inpBaseName:      Write RAOs into a *.fchk file; Default: '' --> Do not write
def wfrt_spec_spin( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
        ixLS, sciform=False, raoType='uni', inpBaseName='', \
        raoFlipAtoms=None ):
    # Treat spin-unrestricted wave functions:
    CLMOa, CLMOb = CLMO
    ELMOa, ELMOb = ELMO
    moIx1_a, moIx1_b = moIx1

    # Vectorize atIx and make moIx strat from ZERO (not 1):
    atIx = np.array( atIx1 ) # Starting from 1
    moIx_a = np.array( moIx1_a ) - 1 # Starting from 0
    moIx_b = np.array( moIx1_b ) - 1 # Starting from 0

    # Number of atoms:
    NAt = len( atIx )
    # Number of electrons:
    NEa = len( moIx_a )
    NEb = len( moIx_b )
    NE = NEa + NEb

    # Choose all LMOs involved in the resoncance subsystem:
    C_LMOa = CLMOa[:,moIx_a]
    C_LMOb = CLMOb[:,moIx_b]
    # Make alpha-LMOs and beta-LMOs match each other:
    #C_LMOa, C_LMOb = matchAlphaBetaLMOs( C_LMOa, C_LMOb )
    #print( C_LMOa.T @ C_LMOb )
    #input()

    if type( naoInfo ) is list: # Interface for wfrt_hmo()
        ifHMO = True
        aoIx = naoInfo
    else:                       # Normal usage
        ifHMO = False
        aoIx = naoInfo.naoIx

    #NP = NE // 2 # Number of electron pairs

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
        atIx = np.array( atIx, dtype='i' )
        LP = list( map(lambda x: atIx[x-1], LP0 ) )
        BD = []
        for arr in BD0:
            BD.append(  np.array( list( map(lambda x: atIx[x-1], arr) ), \
                                  dtype='i' )  )

    #-- OPTION 2: Straightforward specification ----------
    # OPEN-SHELL CASE:
    else:
        ifEnume = False
        LP, BD, RA, RB = ixLS

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

        RA0 = list( map(lambda x: sqIx[x-1]+1, RA ) )
        RB0 = list( map(lambda x: sqIx[x-1]+1, RB ) )
    #------------------------------------------------------------------------

    NLS = len( LP )  # Total number of Lewis structures
    print( '%i Lewis structures in total to be considered' % NLS )

    # Density matrix from pi-LMOs in NAO basis:
    D0a, D0b = densityMatrixLMO( ( C_LMOa, C_LMOb ) )
    # Total energy of the LMOs for the resonance subsystem:
    Etot_LMOa = LMOEnergy( ELMOa, moIx_a )
    Etot_LMOb = LMOEnergy( ELMOb, moIx_b )
    Etot_LMO = Etot_LMOa + Etot_LMOb
    # Fock matrix in the original AO basis:
    Fa, Fb = FockMatrix( FchkInfo.C, FchkInfo.E )
    # Fock matrix in NAO basis:
    FNAOa = CAONAO.T @ Fa @ CAONAO
    FNAOb = CAONAO.T @ Fb @ CAONAO

    # Resonance atomic orbitals (RAOs) in AO basis (NAOs, Lowdin's AOs etc.):
    if raoType == 'uni': # Universal RAOs that are independent of LSs
        if ifHMO:
            CRAOa = np.identity( NAt )
            CRAOb = np.identity( NAt )
        else:
            CRAOa = rao_uni( D0a, aoIx, raoFlipAtoms )
            CRAOb = rao_uni( D0b, aoIx, raoFlipAtoms )
    elif raoType == 'uni-gen': # General universal RAOs
        raise NotImplementedError( 'General universal RAOs not supported yet '
                'in WFRT' )
    else:
        raise ValueError( 
                'Unrecognized value for argument raoType in function wfrt()' )

    # Write RAOs into a *.fchk file:
    if len( inpBaseName ) > 0:
        raoOutpName = inpBaseName + '_RAO.fchk'
        CAORAOa = CAONAO @ CRAOa
        CAORAOb = CAONAO @ CRAOb
        FRAOa = CAORAOa.T @ Fa @ CAORAOa
        FRAOb = CAORAOb.T @ Fb @ CAORAOb
        E_RAOa = np.diag( FRAOa )
        E_RAOb = np.diag( FRAOb )
        writeFchkOrb( inpBaseName+'.fchk', raoOutpName, \
                ( CAORAOa, CAORAOb ), ( E_RAOa, E_RAOb ) )
        print( 'RAOs written into ' + raoOutpName )

    # Lewis unit orbitals (LUOs) in AO basis (NAOs, Lowdin's AOs etc.):
    print( 'Constructing alpha LUOs ...' )
    atIx_LUO, C_LUOa, S_LUOa, F_LUOa = luo( CRAOa, FNAOa, atIx, NAt, NEa )
    print( 'Constructing beta LUOs ...' )
    atIx_LUO, C_LUOb, S_LUOb, F_LUOb = luo( CRAOb, FNAOb, atIx, NAt, NEb )
    print( 'Summary of LUOs and their Wiberg bond orders and spin densities:' )
    #print( atIx_LUO )
    #print( C_LUOa )
    #print( C_LUOb )
    #print( np.array_equal( C_LUOa, C_LUOb ) )
    N_LUO = len( S_LUOa )

    # Get spin density for 1c-LUOs and bond order for 2c-LUOs:
    _, _, _, SD, BO, _, _, = popDM_spin( D0a, D0b, aoIx )
    #writeMat( s, '/dev/stdout', '%12.6f' )


    # Print the RAO indices in each LUO:
    print( '-'*33 )
    print( '%5s   %4s   %s    BO/Spin' % ('LUO', 'Type', 'RAO(s)') )
    print( '-'*33 )
    k = 1
    for ix in atIx_LUO:
        print( '%5i' % k, end='' )
        if ix[1] == -1: # LP
            print( '  (LP):  %3i        %7.4f' % ( atIx[ix[0]], SD[ix[0]] ) )
        else:           # BD
            print( '  (BD):  %3i %3i    %7.4f' % 
                    ( atIx[ix[0]], atIx[ix[1]], BO[ix[0],ix[1]] ) )
        k += 1
    print( '-'*33 )

    # Matrix storing the indices of LUOs of a given set of Lewis structures:
    #luoIX = luoIxLewis( atIx_LUO, NE, LP0, BD0 )
    luoIXa, luoIXb = luoIxLewis_open( atIx_LUO, NEa, NEb, LP0, BD0, RA0, RB0 )
    #print( luoIXa )
    #print( luoIXb )
    #input()

    # Overlap matrix between LUOs and LMOs:
    S_LUO_LMOa = np.zeros( ( N_LUO, NEa ) )
    S_LUO_LMOb = np.zeros( ( N_LUO, NEb ) )
    for i in range( N_LUO ): # Running over all LUOs
        for j in range( NEa ): # Running over all alpha LMOs
            S_LUO_LMOa[ i, j ] = C_LUOa[ :, i ].T @ C_LMOa[ :, j ]
        for j in range( NEb ): # Running over all beta LMOs
            S_LUO_LMOb[ i, j ] = C_LUOb[ :, i ].T @ C_LMOb[ :, j ]
    #writeMat( S_LUO_LMOa )
    #writeMat( S_LUO_LMOb )
    #print( S_LUO_LMOa.shape )
    #print( S_LUO_LMOb.shape )
    #input()

    # Projection of Lewis structure determinants onto that of Psi0:
    P = np.zeros( ( NLS, 1 ) ) # Initialize the projections
    Aaa = np.zeros( ( NEa, NEa ) ) # Initialize alpha-alpha determinant
    Abb = np.zeros( ( NEb, NEb ) ) # Initialize beta-beta determinant
    #print( S_LUO_LMOa )
    #print( S_LUO_LMOb )
    #print('x'*20)
    for i in range( NLS ):
        for j in range( NEa ): # Running over all alpha electrons in an LS
            iLUO = luoIXa[ i, j ] # Index of LUOa for this alpha electron
            for iLMO in range( NEa ): # Index of LMOa
                Aaa[ j, iLMO ] = S_LUO_LMOa[ iLUO, iLMO ]
        for j in range( NEb ): # Running over all beta electrons in an LS
            iLUO = luoIXb[ i, j ] # Index of LUOb for this beta electron
            for iLMO in range( NEb ): # Index of LMOb
                Abb[ j, iLMO ] = S_LUO_LMOb[ iLUO, iLMO ]
        #print( i, np.linalg.det( Aaa ), np.linalg.det( Abb ) )
        #input()
        P[i] = np.linalg.det( Aaa ) * np.linalg.det( Abb )
    #writeMat( P, '/dev/stdout', '%12.6f' )

# DO NOT SORT FOR wfrt_spec:
#    # Sort by P:
#    ix_sort = np.argsort( -P[:,0] )
#    P = P[ ix_sort ]

    # Do not sort:
    ix_sort = np.arange( NLS )
    P = P[ ix_sort ]

    # Update LP[], BD[], and luoIX:
    LP = [ LP[i] for i in ix_sort ]
    BD = [ BD[i] for i in ix_sort ]
    RA = [ RA[i] for i in ix_sort ]
    RB = [ RB[i] for i in ix_sort ]
    luoIXa = luoIXa[ ix_sort, : ]
    luoIXb = luoIXb[ ix_sort, : ]

    # Strings for Lewis structures:
    str_LS = []
    for i in range( NLS ):
        str_LS.append( lewis_str( LP[i], BD[i], RA[i], RB[i] ) )
    #print( str_LS )
    #input()

    #------------------------------------------------------------------
    # Print out the reordered list of Lewis structures:
    #------------------------------------------------------------------
    print( '-'*80 )
    print( '%5s %15s    ' % ( 'No.', 'Projection'), end='' )
    print( ' '*NE, end='' )
    print( ' LUOs', end='' )
    print( ' '*NE, end='' )
    print( '  Lewis structure' )
    print( '-'*80 )
    for i in range( NLS ):
        print( '%5i %15.10f   ' % ( i+1, P[i] ), end='' )
        for ix in luoIXa[i]:
            print( ' %4i' % ( ix+1 ), end='' )
        print( '     %s' % str_LS[i] )
    print( '-'*80 )
    #------------------------------------------------------------------

    # NOTE: Since the user specified set is usually small and not complete,
    #       we DO NOT apply Rumer's rule to deal with linearly dependence

    #------------------------------------------------------------------
    # Solve linear equation for coeff. using ridge linear regression:
    #------------------------------------------------------------------
    print( 'Solving linear equation for WFRT coefficients ...' )

    # Only choose the specified Lewis structures:
    if ifEnume:
        ix_reducLS = np.array( ixLS, dtype='i' ) - 1
        print( 'Using the %i specified Lewis structures to expand the '
               'actual wave function' % ( len( ix_reducLS ) ) )
    else:
        ix_reducLS = np.arange( NLS )
        print( 'Using %i straightforwardly specified Lewis structures to '
               'expand the actual wave function' % ( len( ix_reducLS ) ) )

    # Reduce P[] and luoIX[]:
    Pr = P[ ix_reducLS ]
    luoIXr_a = luoIXa[ ix_reducLS, : ]
    luoIXr_b = luoIXb[ ix_reducLS, : ]

    # Overlap matrix between determinants of all considered Lewis structures:
    #Sr, rank = overlapLewis( NE, luoIXr, S_LUO )
    Sr, Fr, rank = overlapFockLewis_spin( NEa, NEb, luoIXr_a, luoIXr_b, \
            S_LUOa, S_LUOb, F_LUOa, F_LUOb )
    #writeMat( Sr, '/dev/stdout', ' %4.2f' )
    print( '%i linearly independent Lewis structures' % rank )
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
    

    # Use ridge linear regression:
    a = 1E-10
    Cr = np.linalg.pinv( Sr@Sr + a*np.eye(len(Pr)) ) @ Sr @ Pr

    # Normalization:
    NormFactor = sqrt( Cr.T @ Sr @ Cr )  # <Psi|Psi>
    print( 'Normalizing the wave function by a factor of '
            '%.10f ...' % (NormFactor) )
    Cr /= NormFactor
    assert( abs( Cr.T @ Sr @ Cr - 1. ) < 1E-10 )


    # Reproducibility:
    R = Pr.T @ Cr

    # Weights:
    if rank_defi == 0:
        Wr_Mull = weightMulliken( Cr, Sr )
        Wr_Bick = weightBickelhaupt( Cr, Sr )
        Wr_RS = weightRos_Schuit( Cr )
        Wr_Lowd = weightLowdin( Cr, Sr )
        Wr_P2WSO = weightP2WSO( Cr, Sr, Pr )
        #writeMat( Wr_Lowd, '/dev/stdout', '%12.6f' )

    # Recover the full coefficient vector by setting zeros for unchosen LSs:
    C = np.zeros( NLS ) # Initialize
    if rank_defi == 0:
        W_Mull = np.zeros( NLS ) # Initialize
        W_Bick = np.zeros( NLS ) # Initialize
        W_RS = np.zeros( NLS ) # Initialize
        W_Lowd = np.zeros( NLS ) # Initialize
        W_P2WSO = np.zeros( NLS ) # Initialize
    for i in range( NLS ):
        ix_reduc = np.nonzero( ix_reducLS == i )[0]
        if len( ix_reduc ) > 0:
            C[i] = Cr[ ix_reduc ]
            if rank_defi == 0:
                W_Mull[i] = Wr_Mull[ ix_reduc ]
                W_Bick[i] = Wr_Bick[ ix_reduc ]
                W_RS[i] = Wr_RS[ ix_reduc ]
                W_Lowd[i] = Wr_Lowd[ ix_reduc ]
                W_P2WSO[i] = Wr_P2WSO[ ix_reduc ]
    #assert( abs( C.T @ S @ C - 1. ) < 1E-10 )


    # Energies of non-interacting Lewis structures:
    N_reducLS = len( ix_reducLS )
    Er = np.diag( Fr ).reshape( N_reducLS, 1 )
    Erel_r = Er - Etot_LMO
    Er_min = np.min( Er )
    REr = Er - Er_min # Relative energies of LSs with respect to Er_min
    # The Er_min gives the Pauling--Wheland resonance energy:
    if ifHMO:
        Ereson = Er_min - Etot_LMO
    else:
        Ereson = au2kcal( Er_min - Etot_LMO )
    ix_min_LS = []
    for k in range( N_reducLS ):
        if ( not ifHMO and abs( Er[k] - Er_min ) < 0.001/au2kcal(1.) ) \
            or ( ifHMO and abs( Er[k] - Er_min ) < 0.001 ):
            ix_min_LS.append( k+1 )

    # Total energy of the expanded wave function:
    Esum = 0.
    for i in range( N_reducLS ):
        Esum += Cr[i] * Cr[i] * Fr[i,i]
        for j in range( i+1, N_reducLS ):
            if j == i:
                continue
            Esum += 2 * Cr[i] * Cr[j] * Fr[i,j]

    # Effective resonance energies of LSs using Mulliken-like partition scheme:
    REeff_r = Cr * Cr * REr
    for i in range( N_reducLS ):
        for j in range( N_reducLS ):
            if j == i:
                continue
            REeff_r[i] += Cr[i] * Cr[j] * ( Fr[i,j] - Er_min * Sr[i,j] )

    # Recover the full energy information by setting zeros for dependent LSs:
    Erel = np.zeros( NLS ) # Initialize
    RE = np.zeros( NLS ) # Initialize
    REeff = np.zeros( NLS ) # Initialize
    for i in range( NLS ):
        ix_reduc = np.nonzero( ix_reducLS == i )[0]
        if len( ix_reduc ) > 0:
            Erel[i] = Erel_r[ ix_reduc ]
            RE[i] = REr[ ix_reduc ]
            REeff[i] = REeff_r[ ix_reduc ]
    #------------------------------------------------------------------


    #------------------------------------------------------------------
    # Print out the final results:
    #------------------------------------------------------------------
    ## Rearrange in descending order of P2WSO weights:
    # Rearrange in descending order of Lowdin weights:
    if rank_defi == 0:
        ix_sort = np.argsort( -W_P2WSO )
        #ix_sort = np.argsort( -W_Lowd )
    else:
        ix_sort = np.argsort( Erel )
    print( '-'*105 )
    if sciform:
        print( '%5s %15s %11s' % ('No.', 'Projection', 'Coefficient'), end='' )
    else:
        print( '%5s %10s %11s' % ('No.', 'Projection', 'Coefficient'), end='' )
    print( ' %7s' % 'Erel', end='' )
    print( ' %7s' % 'REeff', end='' )
    if rank_defi == 0:
        print( ' %7s' % 'Mulli.', end='' )
        print( ' %7s' % 'Bickel.', end='' )
        print( ' %7s' % 'Ros-Sc.', end='' )
        print( ' %7s' % 'Lowdin', end='' )
        if sciform:
            print( ' %15s' % 'P2WSO', end='' )
        else:
            print( ' %7s' % 'P2WSO', end='' )
    print( '  Lewis structure' )
    print( '-'*105 )
    for iLS in range( NLS ):
        i = ix_sort[ iLS ]
        # Only show Lewis structures in the reduced set:
        if np.count_nonzero( ix_reducLS == i ) == 0:
            continue
        if sciform:
            print( '%5i %15.8E %11.7f' % ( i+1, P[i], C[i] ), end='' )
        else:
            print( '%5i %10.6f %11.7f' % ( i+1, P[i], C[i] ), end='' )
        if ifHMO:
            print( ' %7.3f' % ( Erel[i] ), end='' )
            print( ' %7.3f' % ( REeff[i] ), end='' )
        else:
            print( ' %7.2f' % ( au2kcal( Erel[i] ) ), end='' )
            print( ' %7.2f' % ( au2kcal( REeff[i] ) ), end='' )
        if rank_defi == 0:
            print( ' %6.2f%%' % ( W_Mull[i]*100 ), end='' )
            print( ' %6.2f%%' % ( W_Bick[i]*100 ), end='' )
            print( ' %6.2f%%' % ( W_RS[i]*100 ), end='' )
            print( ' %6.2f%%' % ( W_Lowd[i]*100 ), end='' )
            if sciform:
                print( ' %15.8E%%' % ( W_P2WSO[i]*100 ), end='' )
            else:
                print( ' %6.2f%%' % ( W_P2WSO[i]*100 ), end='' )
        print( '  %s' % str_LS[i] )
    print( '-'*105 )
    print( 'Reproducibility = %10.6f%%' % ( R*100 ) )

    if ifHMO:
        print( 'Total HF energy of reference structure is %.4f |beta|' 
                % Etot_LMO )
        print( 'Total energy of expanded wave function is %.4f |beta|' 
            % Esum )
    else:
        print( 'Total HF energy of reference structure is %.8E a.u.' 
                % Etot_LMO )
        print( 'Total energy of expanded wave function is %.8E a.u.' 
            % Esum )
    print( '  Percent error = %.6f%% ' 
            % ( ( Esum/Etot_LMO - 1 ) * 100 * np.sign( Etot_LMO ) ) )
    if ifHMO:
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.4f |beta|' % Ereson )
        Ereson_expan =  -REeff.sum()
        print( 'Pauling--Wheland resonance energy given by expanded wave '
               'function is %.4f |beta|' % Ereson_expan )
    else:
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.2f kcal/mol' % Ereson )
        Ereson_expan =  -au2kcal( REeff.sum() )
        print( 'Pauling--Wheland resonance energy given by expanded wave '
               'function is %.2f kcal/mol' % Ereson_expan )
    print( '  Percent error = %.6f%% ' % ( ( 1 - Ereson_expan/Ereson )*100 ) )
    if len( ix_min_LS ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_LS:
        print( ' #%i' % i, end='' )
    print()
    #------------------------------------------------------------------


    # Write overlap matrix to external file:
    if globals.ifWrite_Overlap:
        # Write the unsorted OVERLAP:
        overlapFile = globals.basename + '.OVERLAP'
        print( 'Writing overlap matrix to file %s ...' % overlapFile )
        np.savetxt( overlapFile, Sr, fmt='%24.15E' )
        print( 'File %s written' % overlapFile )

        # Write the sorted OVERLAP:
        overlapFile = globals.basename + '.OVERLAP_SORT'
        print( 'Writing overlap matrix to file %s ...' % overlapFile )
        Sr_sortByW = np.zeros( ( NLS, NLS ) )
        for iLS in range( NLS ):
            i = ix_sort[ iLS ]
            for jLS in range( NLS ):
                j = ix_sort[ jLS ]
                Sr_sortByW[iLS][jLS] = Sr[i][j]
        np.savetxt( overlapFile, Sr_sortByW, fmt='%24.15E' )
        print( 'File %s written' % overlapFile )

    # Write Hamiltonian matrix to external file:
    if globals.ifWrite_Hamilt:
        # Write the unsorted HAMILTONIAN:
        HamiltonianFile = globals.basename + '.HAMILT'
        print( 'Writing Hamiltonian matrix to file %s ...' % HamiltonianFile )
        np.savetxt( HamiltonianFile, Fr, fmt='%24.15E' )
        print( 'File %s written' % HamiltonianFile )

        # Write the sorted HAMILTONIAN:
        HamiltonianFile = globals.basename + '.HAMILT_SORT'
        print( 'Writing Hamiltonian matrix to file %s ...' % HamiltonianFile )
        Fr_sortByW = np.zeros( ( NLS, NLS ) )
        for iLS in range( NLS ):
            i = ix_sort[ iLS ]
            for jLS in range( NLS ):
                j = ix_sort[ jLS ]
                Fr_sortByW[iLS][jLS] = Fr[i][j]
        np.savetxt( HamiltonianFile, Fr_sortByW, fmt='%24.15E' )
        print( 'File %s written' % HamiltonianFile )

    return
# enddef wfrt_spec_spin()


# Calculate the projections of Lewis structures using both DMRT and WFRT:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   CLMO:        LMO coefficients in NAO basis
#   ELMO:        LMO energies
#   atIx1:         Indices of atoms in the resonance subsystem, start from 1
#   moIx1:         Indices of LMOs in the resonance subsystem, start from 1
#   ixLS:         Lewis structures specified in either of the following ways:
#                 (1) List of indices of the specified Lewis structures
#                 This choice involves enumeration of all possbile LSs. If the
#                 enumeration has been done previously, file LEWIS_Nc_Ne.dat
#                 will be read in.
#                 NOTE: 1) Indices start from 1 (not 0)
#                       2) The list is from ordering of projections, determined
#                          by calling the wfrt() function prior to calling this
#                          function
#                 (2) Straightfoward specification given by LP[] and BD[],
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[] and BD[], atomic indices start from 1 
#                       (not 0), and they correspond to ACTUAL indices of 
#                       atoms, which may not be consecutive, not as those in 
#                       LP0[] and BD0[] (also see notes in body of this func)
#                 (3) If ixLS == -1 (default), no LSs are specified
#   inpBaseName:      Write RAOs into a *.fchk file; Default: '' --> Do not write
def proj_DM_WF( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                ixLS=-1, raoType='uni', inpBaseName='', raoFlipAtoms=None ):
    # Vectorize atIx and make moIx strat from ZERO (not 1):
    atIx = np.array( atIx1 ) # Starting from 1
    moIx = np.array( moIx1 ) - 1 # Starting from 0

    # Number of atoms:
    NAt = len( atIx )
    # Number of electrons:
    NE = len( moIx ) * 2  # Closed-shell system

    # Choose all LMOs involved in the resoncance subsystem:
    C_LMO = CLMO[:,moIx]
#    writeMat( C_LMO, '/dev/stdout', '%12.6f' )

    # Density matrix of the actual system:
    D0 = densityMatrixLMO( C_LMO )
    #writeMat( D0, 'D0', ' %20.10E' )
    #writeMat( D0, '/dev/stdout', ' %12.6f' )
    # Total energy of the LMOs for the resonance subsystem:
    Etot_LMO = LMOEnergy( ELMO, moIx )
    # Fock matrix in the original AO basis:
    F = FockMatrix( FchkInfo.C, FchkInfo.E )
    # Fock matrix in NAO basis:
    FNAO = CAONAO.T @ F @ CAONAO
    #writeMat( FNAO, '/dev/stdout', ' %12.6f' )

    if type( naoInfo ) is list: # Interface for wfrt_hmo()
        ifHMO = True
        aoIx = naoInfo
    else:                       # Normal usage
        ifHMO = False
        aoIx = naoInfo.naoIx
    NP = NE // 2 # Number of electron pairs

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
        atIx = np.array( atIx, dtype='i' )
        LP = list( map(lambda x: atIx[x-1], LP0 ) )
        BD = []
        for arr in BD0:
            BD.append(  np.array( list( map(lambda x: atIx[x-1], arr) ), \
                                  dtype='i' )  )

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

    # Only choose the specified Lewis structures:
    NLS = len( LP )  # Total number of Lewis structures
    if ifEnume:
        print( 'Totally %i possible Lewis structures' % NLS )
        if type(ixLS) is int:
            if ixLS == -1:
                ixLS = [ -1 ]
            else:
                ixLS = [ -ixLS ]
        if ixLS[0] <= 0:
            ix_reducLS = np.arange( len( LP ) )
        else:
            ix_reducLS = np.array( ixLS, dtype='i' ) - 1
        NLS = len( ix_reducLS ) # Number of Lewis structures to be considered
        # Reduce LP[] and BD[]:
        LP = [ LP[i] for i in ix_reducLS ] # P is a list
        BD = [ BD[i] for i in ix_reducLS ] # BD is a list

    print( 'Using %i straightforwardly specified Lewis structures' % NLS )

    # Resonance atomic orbitals (RAOs) in AO basis (NAOs, Lowdin's AOs etc.):
    if raoType == 'uni': # Universal RAOs that are independent of LSs
        CRAO = rao_uni( D0, aoIx, raoFlipAtoms )
    elif raoType == 'uni-gen': # General universal RAOs
        raise NotImplementedError( 'General universal RAOs not supported yet '
                'in WFRT' )
    else:
        raise ValueError( 
                'Unrecognized value for argument raoType in function '
                'proj_DM_WF()' )

    # Write RAOs into a *.fchk file:
    if len(inpBaseName) > 0:
        raoOutpName = inpBaseName + '_RAO.fchk'
        CAORAO = CAONAO @ CRAO
        FRAO = CAORAO.T @ F @ CAORAO
        E_RAO = np.diag( FRAO )
        writeFchkOrb( inpBaseName+'.fchk', raoOutpName, CAORAO, E_RAO )
        print( 'RAOs written into ' + raoOutpName )

    # Lewis unit orbitals (LUOs) in AO basis (NAOs, Lowdin's AOs etc.):

    #print( '-'*80 )
    #writeMat( CRAO, '/dev/stdout', ' %12.6f' )
    #input()


    # Density matrix of RAOs in NAO (or Lowdin etc.) basis:
    print( 'Calculating densitry matrices for the Lewis strctures ...' )
    D = [] 
    for k in range( NLS ):
        D.append( DM_RAO_uni( CRAO, LP[k], BD[k] ) )

    #writeMat( D[0], '/dev/stdout', ' %12.6f' )
    #print( '-'*60 )
    #writeMat( D[1], '/dev/stdout', ' %12.6f' )
    #input()

    # Energies of Lewis structures:
    print( 'Calculating energies for the Lewis strctures ...' )
    E = [] 
    for k in range( NLS ):
        E.append( np.trace( D[k] @ FNAO ) )
    E = np.c_[ E ]


    #------------------------- WFRT calculations -------------------------#
    print( '='*80 )
    print( 'Performing Wave Function based Resonance Theory (WFRT) '
            'calculations ...' )
    print( '='*80 )
    # Lewis unit orbitals (LUOs) in AO basis (NAOs, Lowdin's AOs etc.):
    atIx_LUO, C_LUO, S_LUO, F_LUO = luo( CRAO, FNAO, atIx, NAt, NE )
    N_LUO = len( S_LUO )
    # Print the RAO indices in each LUO:
    print( '-'*25 )
    print( '%5s   %4s   %s' % ('LUO', 'Type', 'RAO(s)') )
    print( '-'*25 )
    k = 1
    for ix in atIx_LUO:
        print( '%5i' % k, end='' )
        if ix[1] == -1: # LP
            print( '  (LP):  %3i' % atIx[ix[0]] )
        else:           # BD
            print( '  (BD):  %3i %3i' % (atIx[ix[0]],atIx[ix[1]]) )
        k += 1
    print( '-'*25 )

    # Matrix storing the indices of LUOs of a given set of Lewis structures:
    luoIX = luoIxLewis( atIx_LUO, NE, LP0, BD0 )

    # Overlap matrix between LUOs and LMOs:
    print( 'Calculating overlap matrix between LUOs and LMOs ...' )
    S_LUO_LMO = np.zeros( ( N_LUO, NP ) )
    for i in range( N_LUO ): # Running over all LUOs
        for j in range( NP ): # Running over all LMOs
            S_LUO_LMO[ i, j ] = C_LUO[ :, i ].T @ C_LMO[ :, j ]
    #writeMat( S_LUO_LMO )

    # Projection of Lewis structure determinants onto that of Psi0:
    print( 'Calculating projections of Lewis structure wave functions onto '
           'the actual wave function ...' )
    Pwf = np.zeros( ( NLS, 1 ) ) # Initialize the projections
    A = np.zeros( ( NP, NP ) ) # Initialize single-spin determinant
    for i in range( NLS ):
        for j in range( NP ): # Running over all electron pairs in an LS
            iLUO = luoIX[ i, j ] # Index of LUO for this electron pair
            for iLMO in range( NP ): # Index of LMO
                A[ j, iLMO ] = S_LUO_LMO[ iLUO, iLMO ]
        Pwf[i] = np.linalg.det( A ) ** 2
    #writeMat( Pwf, '/dev/stdout', '%12.6f' )

    #----------------------------------------------------------------------

    #------------------------- DMRT calculations -------------------------#
    print( '='*80 )
    print( 'Performing Density Matrix based Resonance Theory (DMRT) '
            'calculations ...' )
    print( '='*80 )

    # Projections onto the actual DM:
    print( 'Calculating projections of Lewis structure density matrices onto '
           'the actual density matrix ...' )
    Pdm = projectDM( D, D0 )
    Pdm /= 2*NE  # Normalized
    #----------------------------------------------------------------------


    # Sort:
    if ifHMO:  # Sort by DM projection of Lewis structures
        ix_sort = np.argsort( -Pdm[:,0] )
    else:      # Sort by energy of Lewis structures
        ix_sort = np.argsort( E[:,0] )
    E = E[ ix_sort ]
    Pwf = Pwf[ ix_sort ]
    Pdm = Pdm[ ix_sort ]
    LP = [ LP[i] for i in ix_sort ]
    BD = [ BD[i] for i in ix_sort ]
    luoIX = luoIX[ ix_sort, : ]
    D = [ D[i] for i in ix_sort ]

    # Lewis structure energies relative to the LMO total energy:
    if ifHMO:
        Etot_LMO = Etot_LMO
        Erel = E - Etot_LMO
    else:
        Etot_LMO = au2kcal( Etot_LMO )
        Erel = au2kcal( E ) - Etot_LMO
    # The lowest Erel is the Pauling--Wheland resonance energy
    Ereson = min( Erel )
    ix_min_LS = []
    for k in range( NLS ):
        if abs( Erel[k] - Ereson ) < 0.001:
            ix_min_LS.append( k+1 )

    if ifHMO:
        print( 'HF Energy of reference structure is %.2f |beta|' % Etot_LMO )
        print( 'Pauling--Wheland resonance energy is %.2f |beta|' % Ereson )
    else:
        print( 'HF Energy of reference structure is %.2f kcal/mol' % Etot_LMO )
        print( 'Pauling--Wheland resonance energy is %.2f kcal/mol' % Ereson )
    if len( ix_min_LS ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_LS:
        print( ' #%i' % i, end='' )
    print()

    # Strings for Lewis structures:
    str_LS = []
    for i in range( NLS ):
        str_LS.append( lewis_str( LP[i], BD[i] ) )



    #------------------------------------------------------------------
    # Print out the reordered list of Lewis structures:
    #------------------------------------------------------------------
    print( '-'*80 )
    print( '%5s' % 'No.', end='' )
    print( ' %12s' % 'Proj (WF)', end='' )
    print( ' %12s' % 'Proj (DM)', end='' )
    if ifHMO:
        print( ' %10s' % 'E (|beta|)', end='' )
    else:
        print( ' %10s' % 'E (kcal)', end='' )
    print( '    Lewis structure' )
    print( '-'*80 )
    for i in range( NLS ):
        print( '%5i' % ( i+1 ), end='' )
        print( ' %12.10f' % Pwf[i], end='' )
        print( ' %12.10f' % Pdm[i], end='' )
        if ifHMO:
            print( ' %10.4f' % Erel[i], end='' )
        else:
            print( ' %10.2f' % Erel[i], end='' )
        print( '    %s' % str_LS[i] )
    print( '-'*80 )
    #------------------------------------------------------------------

    return
# enddef proj_DM_WF()


# Spin-unrestricted version of proj_DM_WF():
# Calculate the projections of Lewis structures using both DMRT and WFRT:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   CLMO:        LMO coefficients in NAO basis
#   ELMO:        LMO energies
#   atIx1:         Indices of atoms in the resonance subsystem, start from 1
#   moIx1:         Indices of LMOs in the resonance subsystem, start from 1
#   ixLS:         Lewis structures specified in either of the following ways:
#                 (1) List of indices of the specified Lewis structures
#                 This choice involves enumeration of all possbile LSs. If the
#                 enumeration has been done previously, file LEWIS_Nc_Ne.dat
#                 will be read in.
#                 NOTE: 1) Indices start from 1 (not 0)
#                       2) The list is from ordering of projections, determined
#                          by calling the wfrt() function prior to calling this
#                          function
#                 (2) Straightfoward specification given by LP[] and BD[],
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[] and BD[], atomic indices start from 1 
#                       (not 0), and they correspond to ACTUAL indices of 
#                       atoms, which may not be consecutive, not as those in 
#                       LP0[] and BD0[] (also see notes in body of this func)
#                 (3) If ixLS == -1 (default), no LSs are specified
#   inpBaseName:      Write RAOs into a *.fchk file; Default: '' --> Do not write
def proj_DM_WF_spin( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                ixLS=-1, raoType='uni', inpBaseName='', raoFlipAtoms=None ):
    # Treat spin-unrestricted wave functions:
    CLMOa, CLMOb = CLMO
    ELMOa, ELMOb = ELMO
    moIx1_a, moIx1_b = moIx1

    # Vectorize atIx and make moIx strat from ZERO (not 1):
    atIx = np.array( atIx1 ) # Starting from 1
    moIx_a = np.array( moIx1_a ) - 1 # Starting from 0
    moIx_b = np.array( moIx1_b ) - 1 # Starting from 0

    # Number of atoms:
    NAt = len( atIx )
    # Number of electrons:
    NEa = len( moIx_a )
    NEb = len( moIx_b )
    NE = NEa + NEb

    # Choose all LMOs involved in the resoncance subsystem:
    C_LMOa = CLMOa[:,moIx_a]
    C_LMOb = CLMOb[:,moIx_b]
    # Make alpha-LMOs and beta-LMOs match each other:
    #C_LMOa, C_LMOb = matchAlphaBetaLMOs( C_LMOa, C_LMOb )
    #print( C_LMOa.T @ C_LMOb )
    #input()

    # Density matrix of the actual system in NAO basis:
    D0a, D0b = densityMatrixLMO( ( C_LMOa, C_LMOb ) )
    D0 = D0a + D0b
    # Total energy of the LMOs for the resonance subsystem:
    Etot_LMOa = LMOEnergy( ELMOa, moIx_a )
    Etot_LMOb = LMOEnergy( ELMOb, moIx_b )
    Etot_LMO = Etot_LMOa + Etot_LMOb
    # Fock matrix in the original AO basis:
    Fa, Fb = FockMatrix( FchkInfo.C, FchkInfo.E )
    # Fock matrix in NAO basis:
    FNAOa = CAONAO.T @ Fa @ CAONAO
    FNAOb = CAONAO.T @ Fb @ CAONAO

    if type( naoInfo ) is list: # Interface for wfrt_hmo()
        ifHMO = True
        aoIx = naoInfo
    else:                       # Normal usage
        ifHMO = False
        aoIx = naoInfo.naoIx
    NP = NE // 2 # Number of electron pairs

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
        atIx = np.array( atIx, dtype='i' )
        LP = list( map(lambda x: atIx[x-1], LP0 ) )
        BD = []
        for arr in BD0:
            BD.append(  np.array( list( map(lambda x: atIx[x-1], arr) ), \
                                  dtype='i' )  )

    #-- OPTION 2: Straightforward specification ----------
    # OPEN-SHELL CASE:
    else:
        ifEnume = False
        LP, BD, RA, RB = ixLS

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

        RA0 = list( map(lambda x: sqIx[x-1]+1, RA ) )
        RB0 = list( map(lambda x: sqIx[x-1]+1, RB ) )
    #------------------------------------------------------------------------

    # Only choose the specified Lewis structures:
    NLS = len( LP )  # Total number of Lewis structures
    if ifEnume:
        print( 'Totally %i possible Lewis structures' % NLS )
        if type(ixLS) is int:
            if ixLS == -1:
                ixLS = [ -1 ]
            else:
                ixLS = [ -ixLS ]
        if ixLS[0] <= 0:
            ix_reducLS = np.arange( len( LP ) )
        else:
            ix_reducLS = np.array( ixLS, dtype='i' ) - 1
        NLS = len( ix_reducLS ) # Number of Lewis structures to be considered
        # Reduce LP[] and BD[]:
        LP = [ LP[i] for i in ix_reducLS ] # P is a list
        BD = [ BD[i] for i in ix_reducLS ] # BD is a list

    print( 'Using %i straightforwardly specified Lewis structures' % NLS )

    # Resonance atomic orbitals (RAOs) in AO basis (NAOs, Lowdin's AOs etc.):
    if raoType == 'uni': # Universal RAOs that are independent of LSs
        if ifHMO:
            CRAO = np.identity( NAt )
        else:
            CRAOa = rao_uni( D0a, aoIx, raoFlipAtoms )
            CRAOb = rao_uni( D0b, aoIx, raoFlipAtoms )
    elif raoType == 'uni-gen': # General universal RAOs
        raise NotImplementedError( 'General universal RAOs not supported yet '
                'in WFRT' )
    else:
        raise ValueError( 
                'Unrecognized value for argument raoType in function '
                'proj_DM_WF_spin()' )

    # Write RAOs into a *.fchk file:
    if len(inpBaseName) > 0:
        raoOutpName = inpBaseName + '_RAO.fchk'
        CAORAOa = CAONAO @ CRAOa
        CAORAOb = CAONAO @ CRAOb
        FRAOa = CAORAOa.T @ Fa @ CAORAOa
        FRAOb = CAORAOb.T @ Fb @ CAORAOb
        E_RAOa = np.diag( FRAOa )
        E_RAOb = np.diag( FRAOb )
        writeFchkOrb( inpBaseName+'.fchk', raoOutpName, \
                ( CAORAOa, CAORAOb ), ( E_RAOa, E_RAOb ) )
        print( 'RAOs written into ' + raoOutpName )

    #print( '-'*80 )
    #writeMat( CRAO, '/dev/stdout', ' %12.6f' )
    #input()


    # Density matrix of RAOs in NAO (or Lowdin etc.) basis:
    print( 'Calculating densitry matrices for the Lewis strctures ...' )
    Da = [] 
    Db = [] 
    D = []
    for k in range( NLS ):
        Da_k = DM_RAO_uni_spin( CRAOa, LP[k], BD[k], RA[k] )
        Db_k = DM_RAO_uni_spin( CRAOb, LP[k], BD[k], RB[k] )
        Da.append( Da_k )
        Db.append( Db_k )
        D.append( Da_k + Db_k )

    #writeMat( D[0], '/dev/stdout', ' %12.6f' )
    #print( '-'*60 )
    #writeMat( D[1], '/dev/stdout', ' %12.6f' )
    #input()

    # Energies of Lewis structures:
    print( 'Calculating energies for the Lewis strctures ...' )
    E = [] 
    for k in range( NLS ):
        Ea = np.trace( Da[k] @ FNAOa )
        Eb = np.trace( Db[k] @ FNAOb )
        E.append( Ea + Eb )
    E = np.c_[ E ]


    #------------------------- WFRT calculations -------------------------#
    print( '='*80 )
    print( 'Performing Wave Function based Resonance Theory (WFRT) '
            'calculations ...' )
    print( '='*80 )
    # Lewis unit orbitals (LUOs) in AO basis (NAOs, Lowdin's AOs etc.):
    print( 'Constructing alpha LUOs ...' )
    atIx_LUO, C_LUOa, S_LUOa, F_LUOa = luo( CRAOa, FNAOa, atIx, NAt, NEa )
    print( 'Constructing beta LUOs ...' )
    atIx_LUO, C_LUOb, S_LUOb, F_LUOb = luo( CRAOb, FNAOb, atIx, NAt, NEb )
    print( 'Summary of LUOs and their Wiberg bond orders and spin densities:' )
    #print( atIx_LUO )
    #print( C_LUOa )
    #print( C_LUOb )
    #print( np.array_equal( C_LUOa, C_LUOb ) )
    N_LUO = len( S_LUOa )

    # Get spin density for 1c-LUOs and bond order for 2c-LUOs:
    _, _, _, SD, BO, _, _, = popDM_spin( D0a, D0b, aoIx )
    #writeMat( s, '/dev/stdout', '%12.6f' )


    # Print the RAO indices in each LUO:
    print( '-'*33 )
    print( '%5s   %4s   %s    BO/Spin' % ('LUO', 'Type', 'RAO(s)') )
    print( '-'*33 )
    k = 1
    for ix in atIx_LUO:
        print( '%5i' % k, end='' )
        if ix[1] == -1: # LP
            print( '  (LP):  %3i        %7.4f' % ( atIx[ix[0]], SD[ix[0]] ) )
        else:           # BD
            print( '  (BD):  %3i %3i    %7.4f' %
                    ( atIx[ix[0]], atIx[ix[1]], BO[ix[0],ix[1]] ) )
        k += 1
    print( '-'*33 )

    # Matrix storing the indices of LUOs of a given set of Lewis structures:
    #luoIX = luoIxLewis( atIx_LUO, NE, LP0, BD0 )
    luoIXa, luoIXb = luoIxLewis_open( atIx_LUO, NEa, NEb, LP0, BD0, RA0, RB0 )
    #print( luoIXa )
    #print( luoIXb )
    #input()

    # Overlap matrix between LUOs and LMOs:
    print( 'Calculating overlap matrix between LUOs and LMOs ...' )
    S_LUO_LMOa = np.zeros( ( N_LUO, NEa ) )
    S_LUO_LMOb = np.zeros( ( N_LUO, NEb ) )
    for i in range( N_LUO ): # Running over all LUOs
        for j in range( NEa ): # Running over all alpha LMOs
            S_LUO_LMOa[ i, j ] = C_LUOa[ :, i ].T @ C_LMOa[ :, j ]
        for j in range( NEb ): # Running over all beta LMOs
            S_LUO_LMOb[ i, j ] = C_LUOb[ :, i ].T @ C_LMOb[ :, j ]
    #writeMat( S_LUO_LMOa )
    #writeMat( S_LUO_LMOb )
    #print( S_LUO_LMOa.shape )
    #print( S_LUO_LMOb.shape )
    #input()

    # Projection of Lewis structure determinants onto that of Psi0:
    print( 'Calculating projections of Lewis structure wave functions onto '
           'the actual wave function ...' )
    Pwf = np.zeros( ( NLS, 1 ) ) # Initialize the projections
    Aaa = np.zeros( ( NEa, NEa ) ) # Initialize alpha-alpha determinant
    Abb = np.zeros( ( NEb, NEb ) ) # Initialize beta-beta determinant
    for i in range( NLS ):
        for j in range( NEa ): # Running over all alpha electrons in an LS
            iLUO = luoIXa[ i, j ] # Index of LUOa for this alpha electron
            for iLMO in range( NEa ): # Index of LMOa
                Aaa[ j, iLMO ] = S_LUO_LMOa[ iLUO, iLMO ]
        for j in range( NEb ): # Running over all beta electrons in an LS
            iLUO = luoIXb[ i, j ] # Index of LUOb for this beta electron
            for iLMO in range( NEb ): # Index of LMOb
                Abb[ j, iLMO ] = S_LUO_LMOb[ iLUO, iLMO ]
        #print( i, np.linalg.det( Aaa ), np.linalg.det( Abb ) )
        #input()
        Pwf[i] = np.linalg.det( Aaa ) * np.linalg.det( Abb )
    #writeMat( Pwf, '/dev/stdout', '%12.6f' )
    #input()
    #----------------------------------------------------------------------

    #------------------------- DMRT calculations -------------------------#
    print( '='*80 )
    print( 'Performing Density Matrix based Resonance Theory (DMRT) '
            'calculations ...' )
    print( '='*80 )

    # Projections onto the actual DM:
    print( 'Calculating projections of Lewis structure density matrices onto '
           'the actual density matrix ...' )
    Pdm_a = projectDM( Da, D0a )
    Pdm_a /= NEa  # Normalized
    Pdm_b = projectDM( Db, D0b )
    Pdm_b /= NEb  # Normalized
    Pdm = np.sqrt( Pdm_a * Pdm_b )
    # writeMat( Pdm, '/dev/stdout', '%12.6f' )
    #print( '-'*40 )
    #writeMat( Pdm_a, '/dev/stdout', '%12.6f' )
    #print( '-'*40 )
    #writeMat( Pdm_b, '/dev/stdout', '%12.6f' )
    #print(2*NE)
    #print( innerDM( D[0], D[0] ) )
    #print( innerDM( D0, D0 ) )
    #input()
    #----------------------------------------------------------------------


    #==================================================================
    # 1. Do not sort:
    #==================================================================
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

    # Strings for Lewis structures:
    str_LS = []
    for i in range( NLS ):
        str_LS.append( lewis_str( LP[i], BD[i], RA[i], RB[i] ) )


    print()
    print( '*'*80 )
    print( 'Projections of resonance structures in the original order' )
    #------------------------------------------------------------------
    # Print out the reordered list of Lewis structures:
    #------------------------------------------------------------------
    print( '-'*80 )
    print( '%5s' % 'No.', end='' )
    print( ' %13s' % 'Proj (WF)', end='' )
    print( ' %13s' % 'Proj (DM)', end='' )
    if ifHMO:
        print( ' %10s' % 'E (|beta|)', end='' )
    else:
        print( ' %10s' % 'E (kcal)', end='' )
    print( '    Lewis structure' )
    print( '-'*80 )
    for i in range( NLS ):
        print( '%5i' % ( i+1 ), end='' )
        print( ' %13.10f' % Pwf[i], end='' )
        print( ' %13.10f' % Pdm[i], end='' )
        if ifHMO:
            print( ' %10.4f' % Erel[i], end='' )
        else:
            print( ' %10.2f' % Erel[i], end='' )
        print( '    %s' % str_LS[i] )
    print( '-'*80 )

    if ifHMO:
        print( 'HF Energy of reference structure is %.2f |beta|' % Etot_LMO )
        print( 'Pauling--Wheland resonance energy is %.2f |beta|' % Ereson )
    else:
        print( 'HF Energy of reference structure is %.8f a.u.' % Etot_LMO )
        print( 'Pauling--Wheland resonance energy is %.2f kcal/mol' % Ereson )
    if len( ix_min_LS ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_LS:
        print( ' #%i' % i, end='' )
    print()

    #------------------------------------------------------------------
    #==================================================================



    #==================================================================
    # 2. Sort:
    #==================================================================
    if ifHMO:  # Sort by DM projection of Lewis structures
        ix_sort = np.argsort( -Pdm[:,0] )
    else:      # Sort by energy of Lewis structures
        #ix_sort = np.argsort( E[:,0] )
        ix_sort = np.argsort( -np.abs( Pwf[:,0] ) )
    E = E[ ix_sort ]
    Pwf = Pwf[ ix_sort ]
    Pdm = Pdm[ ix_sort ]
    LP = [ LP[i] for i in ix_sort ]
    BD = [ BD[i] for i in ix_sort ]
    RA = [ RA[i] for i in ix_sort ]
    RB = [ RB[i] for i in ix_sort ]
    #luoIXa = luoIXa[ ix_sort, : ]
    #luoIXb = luoIXb[ ix_sort, : ]

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

    # Strings for Lewis structures:
    str_LS = []
    for i in range( NLS ):
        str_LS.append( lewis_str( LP[i], BD[i], RA[i], RB[i] ) )


    print()
    print( '*'*80 )
    print( 'Projections of resonance structures in descending order '
            'of wave funtion projection:' )
    #------------------------------------------------------------------
    # Print out the reordered list of Lewis structures:
    #------------------------------------------------------------------
    print( '-'*80 )
    print( '%5s' % 'No.', end='' )
    print( ' %20s' % 'Proj (WF)', end='' )
    print( ' %20s' % 'Proj (DM)', end='' )
    if ifHMO:
        print( ' %20s' % 'E (|beta|)', end='' )
    else:
        print( ' %20s' % 'E (kcal)', end='' )
    print( '    Lewis structure' )
    print( '-'*80 )
    for i in range( NLS ):
        print( '%5i' % ( i+1 ), end='' )
        print( ' %20.10E' % Pwf[i], end='' )
        print( ' %20.10E' % Pdm[i], end='' )
        if ifHMO:
            print( ' %20.10E' % Erel[i], end='' )
        else:
            print( ' %20.10E' % Erel[i], end='' )
        print( '    %s' % str_LS[i] )
    print( '-'*80 )
    #------------------------------------------------------------------
    #==================================================================

    if ifHMO:
        print( 'HF Energy of reference structure is %.8f |beta|' % Etot_LMO )
        print( 'Pauling--Wheland resonance energy is %.8f |beta|' % Ereson )
    else:
        print( 'HF Energy of reference structure is %.8E a.u.' % Etot_LMO )
        print( 'Pauling--Wheland resonance energy is %.4f kcal/mol' % Ereson )
    if len( ix_min_LS ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_LS:
        print( ' #%i' % i, end='' )
    print()

    return
# enddef proj_DM_WF_spin()




# Calculate the projections of Lewis structures using WFRT:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   CLMO:        LMO coefficients in NAO basis
#   ELMO:        LMO energies
#   atIx1:         Indices of atoms in the resonance subsystem, start from 1
#   moIx1:         Indices of LMOs in the resonance subsystem, start from 1
#   ixLS:         Lewis structures specified in either of the following ways:
#                 (1) List of indices of the specified Lewis structures
#                 This choice involves enumeration of all possbile LSs. If the
#                 enumeration has been done previously, file LEWIS_Nc_Ne.dat
#                 will be read in.
#                 NOTE: 1) Indices start from 1 (not 0)
#                       2) The list is from ordering of projections, determined
#                          by calling the wfrt() function prior to calling this
#                          function
#                 (2) Straightfoward specification given by LP[] and BD[],
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[] and BD[], atomic indices start from 1 
#                       (not 0), and they correspond to ACTUAL indices of 
#                       atoms, which may not be consecutive, not as those in 
#                       LP0[] and BD0[] (also see notes in body of this func)
#                 (3) If ixLS == -1 (default), no LSs are specified
#   inpBaseName:      Write RAOs into a *.fchk file; Default: '' --> Do not write
def proj_WF( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                ixLS=-1, raoType='uni', inpBaseName='', raoFlipAtoms=None ):
    # Vectorize atIx and make moIx strat from ZERO (not 1):
    atIx = np.array( atIx1 ) # Starting from 1
    moIx = np.array( moIx1 ) - 1 # Starting from 0

    # Number of atoms:
    NAt = len( atIx )
    # Number of electrons:
    NE = len( moIx ) * 2  # Closed-shell system

    # Choose all LMOs involved in the resoncance subsystem:
    C_LMO = CLMO[:,moIx]
#    writeMat( C_LMO, '/dev/stdout', '%12.6f' )

    # Density matrix of the actual system:
    D0 = densityMatrixLMO( C_LMO )
    #writeMat( D0, 'D0', ' %20.10E' )
    #writeMat( D0, '/dev/stdout', ' %12.6f' )
    # Total energy of the LMOs for the resonance subsystem:
    Etot_LMO = LMOEnergy( ELMO, moIx )
    # Fock matrix in the original AO basis:
    F = FockMatrix( FchkInfo.C, FchkInfo.E )
    # Fock matrix in NAO basis:
    FNAO = CAONAO.T @ F @ CAONAO
    #writeMat( FNAO, '/dev/stdout', ' %12.6f' )

    if type( naoInfo ) is list: # Interface for wfrt_hmo()
        ifHMO = True
        aoIx = naoInfo
    else:                       # Normal usage
        ifHMO = False
        aoIx = naoInfo.naoIx
    NP = NE // 2 # Number of electron pairs

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
        maxNLP = ixLS
        if ixLS == -1:
            maxNLP = inf
        print( "Enumerating Lewis structures with at most", maxNLP, \
               "lone pairs"  )
        LP0, BD0 = lewis_read( NAt, NE, maxNLP )

        # Set actual indices of atoms:
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

    # Only choose the specified Lewis structures:
    NLS = len( LP )  # Total number of Lewis structures
    if ifEnume:
        print( 'Totally %i possible Lewis structures' % NLS )
        if type(ixLS) is int:
            if ixLS == -1:
                ixLS = [ -1 ]
            else:
                ixLS = [ -ixLS ]
        if ixLS[0] <= 0:
            ix_reducLS = np.arange( len( LP ) )
        else:
            ix_reducLS = np.array( ixLS, dtype='i' ) - 1
        NLS = len( ix_reducLS ) # Number of Lewis structures to be considered
        # Reduce LP[] and BD[]:
        LP = [ LP[i] for i in ix_reducLS ] # P is a list
        BD = [ BD[i] for i in ix_reducLS ] # BD is a list

    print( 'Using %i straightforwardly specified Lewis structures' % NLS )

    # Resonance atomic orbitals (RAOs) in AO basis (NAOs, Lowdin's AOs etc.):
    if raoType == 'uni': # Universal RAOs that are independent of LSs
        CRAO = rao_uni( D0, aoIx, raoFlipAtoms )
    elif raoType == 'uni-gen': # General universal RAOs
        raise NotImplementedError( 'General universal RAOs not supported yet '
                'in WFRT' )
    else:
        raise ValueError( 
                'Unrecognized value for argument raoType in function '
                'proj_DM_WF()' )

    # Write RAOs into a *.fchk file:
    if len(inpBaseName) > 0:
        raoOutpName = inpBaseName + '_RAO.fchk'
        CAORAO = CAONAO @ CRAO
        FRAO = CAORAO.T @ F @ CAORAO
        E_RAO = np.diag( FRAO )
        writeFchkOrb( inpBaseName+'.fchk', raoOutpName, CAORAO, E_RAO )
        print( 'RAOs written into ' + raoOutpName )

    # Lewis unit orbitals (LUOs) in AO basis (NAOs, Lowdin's AOs etc.):

    #print( '-'*80 )
    #writeMat( CRAO, '/dev/stdout', ' %12.6f' )
    #input()


    # Density matrix of RAOs in NAO (or Lowdin etc.) basis:
    print( 'Calculating densitry matrices for the Lewis strctures ...' )
    D = [] 
    for k in range( NLS ):
        D.append( DM_RAO_uni( CRAO, LP[k], BD[k] ) )

    #writeMat( D[0], '/dev/stdout', ' %12.6f' )
    #print( '-'*60 )
    #writeMat( D[1], '/dev/stdout', ' %12.6f' )
    #input()

    # Energies of Lewis structures:
    print( 'Calculating energies for the Lewis strctures ...' )
    E = [] 
    for k in range( NLS ):
        E.append( np.trace( D[k] @ FNAO ) )
    E = np.c_[ E ]


    #------------------------- WFRT calculations -------------------------#
    print( '='*80 )
    print( 'Performing Wave Function based Resonance Theory (WFRT) '
            'calculations ...' )
    print( '='*80 )
    # Lewis unit orbitals (LUOs) in AO basis (NAOs, Lowdin's AOs etc.):
    atIx_LUO, C_LUO, S_LUO, F_LUO = luo( CRAO, FNAO, atIx, NAt, NE )
    N_LUO = len( S_LUO )
    # Print the RAO indices in each LUO:
    print( '-'*25 )
    print( '%5s   %4s   %s' % ('LUO', 'Type', 'RAO(s)') )
    print( '-'*25 )
    k = 1
    for ix in atIx_LUO:
        print( '%5i' % k, end='' )
        if ix[1] == -1: # LP
            print( '  (LP):  %3i' % atIx[ix[0]] )
        else:           # BD
            print( '  (BD):  %3i %3i' % (atIx[ix[0]],atIx[ix[1]]) )
        k += 1
    print( '-'*25 )

    # Matrix storing the indices of LUOs of a given set of Lewis structures:
    luoIX = luoIxLewis( atIx_LUO, NE, LP0, BD0 )

    # Overlap matrix between LUOs and LMOs:
    print( 'Calculating overlap matrix between LUOs and LMOs ...' )
    S_LUO_LMO = np.zeros( ( N_LUO, NP ) )
    for i in range( N_LUO ): # Running over all LUOs
        for j in range( NP ): # Running over all LMOs
            S_LUO_LMO[ i, j ] = C_LUO[ :, i ].T @ C_LMO[ :, j ]
    #writeMat( S_LUO_LMO )

    # Projection of Lewis structure determinants onto that of Psi0:
    print( 'Calculating projections of Lewis structure wave functions onto '
           'the actual wave function ...' )
    Pwf = np.zeros( ( NLS, 1 ) ) # Initialize the projections
    A = np.zeros( ( NP, NP ) ) # Initialize single-spin determinant
    for i in range( NLS ):
        for j in range( NP ): # Running over all electron pairs in an LS
            iLUO = luoIX[ i, j ] # Index of LUO for this electron pair
            for iLMO in range( NP ): # Index of LMO
                A[ j, iLMO ] = S_LUO_LMO[ iLUO, iLMO ]
        Pwf[i] = np.linalg.det( A ) ** 2
    #writeMat( Pwf, '/dev/stdout', '%12.6f' )

    #----------------------------------------------------------------------


    # Strings for Lewis structures:
    str_LS = []
    for i in range( NLS ):
        str_LS.append( lewis_str( LP[i], BD[i] ) )



    #------------------------------------------------------------------
    # Print out the list of Lewis structures:
    #------------------------------------------------------------------
    print( '-'*80 )
    print( '%5s' % 'No.', end='' )
    print( ' %12s' % 'Proj (WF)', end='' )
    print( '    Lewis structure' )
    print( '-'*80 )
    for i in range( NLS ):
        print( '%5i' % ( i+1 ), end='' )
        print( ' %12.10f' % Pwf[i], end='' )
        print( '    %s' % str_LS[i] )
    print( '-'*80 )
    #------------------------------------------------------------------

    return
# enddef proj_WF()



# WFRT analysis only with all possible Kekule structures:
#   naoInfo:       Information from NAO's output
#   FchkInfo:      Informatoin of fchk output
#   CAONAO:        AO-->NAO transformation matrix
#   CLMO:          LMO coefficients in NAO basis
#   ELMO:          LMO energies
#   atIx1:         Indices of atoms in the resonance subsystem, start from 1
#   moIx1:         Indices of LMOs in the resonance subsystem, start from 1
#   kekFileName:   Name of *.kek file where Kekule structures are stored
#   sciform:      whether showing numerical results in scientific format
#   raoType:       'uni'--> universal RAOs;  'var'--> varying RAOs
def wfrt_kekule( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                 kekFileName, sciform=False, raoType='uni', inpBaseName='' ):
    # Read Kekule structures from external file:
    print( 'Reading Kekule structures from file %s ...' % kekFileName )
    LP, BD, RA, RB, isRadical = readKekule( kekFileName )

    # Check consistency
    if isRadical and ( not FchkInfo.spin ):
        raise ValueError( 'The Kekule structures are a radical, '
                'but the wavefunction is closed-shelled' )
    elif ( not isRadical ) and FchkInfo.spin:
        raise ValueError( 'The Kekule structures are closed-shelled, '
                'but the wavefunction is open-shelled' )

    # Print spin info:
    if isRadical:
        print( 'A spin-unrestricted wave function with spin multiplicity of %i'
                % FchkInfo.mult )
    else:
        print( 'A closed-shell wave function' )

    # Perform WFRT analysis:
    # 1. Closed-shell system:
    if not isRadical:
        wfrt_spec( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
            ( LP, BD ), sciform, raoType, inpBaseName )
    # 2. Spin-unrestricted:
    else:
        wfrt_spec_spin( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
            ( LP, BD, RA, RB ), sciform, raoType, inpBaseName )

    return
# enddef wfrt_kekule()


# Calculate DMRT and WFRT projections only with all possible Kekule structures:
#   naoInfo:     Information from NAO's output
#   FchkInfo:    Informatoin of fchk output
#   CAONAO:      AO-->NAO transformation matrix
#   CLMO:        LMO coefficients in NAO basis
#   ELMO:        LMO energies
#   atIx1:       Indices of atoms in the resonance subsystem, start from 1
#   moIx1:       Indices of LMOs in the resonance subsystem, start from 1
#   kekFileName: Name of *.kek file where Kekule structures are stored
#   ndiv:        Divide file in n parts lest a huge file runs out of memory
def proj_DM_WF_kekule( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                 kekFileName, ndiv=1, raoType='uni', inpBaseName='', \
                 raoFlipAtoms=None ):
    if ndiv > 1:
        print( 'Dividing file %s into %i parts ...' % ( kekFileName, ndiv ) )
    elif ndiv <=0 :
        raise ValueError( 'Invalid value for argument ndiv (%i) in '
                'proj_DM_WF_kekule' % ndiv )

    for i in range( 1, ndiv+1 ):
        if ndiv > 1:
            print( '########## PART %3i OUT OF %3i ##########' % ( i, ndiv ) )

        # Read Kekule structures from external file:
        print( 'Reading Kekule structures from file %s ...' % kekFileName )
        #LP, BD = readKekule( kekFileName, ( ndiv, i ) )
        LP, BD, RA, RB, isRadical = readKekule( kekFileName, ( ndiv, i ) )
        #isRadical = len( RA ) > 0 or len( RB ) > 0

        # Check consistency
        if isRadical and ( not FchkInfo.spin ):
            raise ValueError( 'The Kekule structures are a radical, '
                    'but the wavefunction is closed-shelled' )
        elif ( not isRadical ) and FchkInfo.spin:
            raise ValueError( 'The Kekule structures are closed-shelled, '
                    'but the wavefunction is open-shelled' )

        # Print spin info:
        if isRadical:
            print( 'A spin-unrestricted wave function with '
                    'spin multiplicity of %i' % FchkInfo.mult )
        else:
            print( 'A closed-shell wave function' )

        # Calculate DMRT and WFRT projections:
        # 1. Closed-shell system:
        if not isRadical:
            proj_DM_WF( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, \
                    moIx1, ( LP, BD ), raoType, inpBaseName, raoFlipAtoms )
            # 2. Spin-unrestricted:
        else:
            proj_DM_WF_spin( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, \
                    moIx1, ( LP, BD, RA, RB ), raoType, inpBaseName, \
                    raoFlipAtoms )

    return
# enddef proj_DM_WF_kekule()



# Calculate WFRT projections only with all possible Kekule structures:
#   naoInfo:     Information from NAO's output
#   FchkInfo:    Informatoin of fchk output
#   CAONAO:      AO-->NAO transformation matrix
#   CLMO:        LMO coefficients in NAO basis
#   ELMO:        LMO energies
#   atIx1:       Indices of atoms in the resonance subsystem, start from 1
#   moIx1:       Indices of LMOs in the resonance subsystem, start from 1
#   kekFileName: Name of *.kek file where Kekule structures are stored
#   ndiv:        Divide file in n parts lest a huge file runs out of memory
def proj_WF_kekule( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                 kekFileName, ndiv=1, raoType='uni', inpBaseName='', \
                 raoFlipAtoms=None ):
    if ndiv > 1:
        print( 'Dividing file %s into %i parts ...' % ( kekFileName, ndiv ) )
    elif ndiv <=0 :
        raise ValueError( 'Invalid value for argument ndiv (%i) in '
                'proj_WF_kekule' % ndiv )

    for i in range( 1, ndiv+1 ):
        if ndiv > 1:
            print( '########## PART %3i OUT OF %3i ##########' % ( i, ndiv ) )

        # Read Kekule structures from external file:
        print( 'Reading Kekule structures from file %s ...' % kekFileName )
        # LP, BD = readKekule( kekFileName, ( ndiv, i ) )
        LP, BD, RA, RB, isRadical = readKekule( kekFileName )

        # Calculate DMRT and WFRT projections:
        proj_WF( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                    ( LP, BD ), raoType, inpBaseName, raoFlipAtoms )
    return
# enddef proj_WF_kekule()


# WFRT analysis in the simple Hueckel Molecular Orbital (HMO) framework:
#   hmoSol:      Solution of simple HMO theory
#   atIx1:         Indices of atoms in the resonance subsystem, start from 1
#   MAX_LP:       Maximum number of lone pairs to be considered
#   projCut:      Cutoff value of projections relative to the principal LS
def wfrt_hmo( hmoSol, atIx1, MAX_LP=inf, projCut=0. ):
    NAt = hmoSol.NAt
    NP = hmoSol.NE // 2  # Number of electron pairs (for closed-shell systems)
    aoIx = []
    for i in range( 1, NAt+1 ):
        aoIx.append( np.array([ i ]) )

    CAONAO = np.eye( NAt )
    CLMO = hmoSol.C[ :, 0:NP ]
    ELMO = hmoSol.E[ 0:NP ]
    moIx1 = range( 1, NP+1 )

    wfrt( aoIx, hmoSol, CAONAO, CLMO, ELMO, atIx1, moIx1, \
            MAX_LP, projCut )
    return
# enddef wfrt_hmo()


# WFRT analysis in the simple Hueckel Molecular Orbital (HMO) framework with 
# a specified set of Lewis structures:
#   hmoSol:      Solution of simple HMO theory
#   ixLS:         Lewis structures specified in either of the following ways:
#                 (1) List of indices of the specified Lewis structures
#                 This choice involves enumeration of all possbile LSs. If the
#                 enumeration has been done previously, file LEWIS_Nc_Ne.dat
#                 will be read in.
#                 NOTE: 1) Indices start from 1 (not 0)
#                       2) The list is from ordering of projections, determined
#                          by calling the wfrt_hmo() function prior to calling 
#                          this function
#                 (2) Straightfoward specification given by LP[] and BD[],
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[] and BD[], atomic indices start from 1 
#                       (not 0), and they correspond to ACTUAL indices of 
#                       atoms, which may not be consecutive, not as those in 
#                       LP0[] and BD0[] (also see notes in body of this func)
#   sciform:      whether showing numerical results in scientific format
def wfrt_hmo_spec( hmoSol, ixLS, sciform=False ):
    LP, BD = ixLS

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

    wfrt_spec( aoIx, hmoSol, CAONAO, CLMO, ELMO, atIx1, moIx1, ixLS, sciform )
    return
# enddef wfrt_hmo_spec()


# Spin-polarized version of wfrt_hmo_spec():
# WFRT analysis in the simple Hueckel Molecular Orbital (HMO) framework with 
# a specified set of Lewis structures:
#   hmoSol:      Solution of simple HMO theory
#   ixLS:         Lewis structures specified in either of the following ways:
#                 (1) List of indices of the specified Lewis structures
#                 This choice involves enumeration of all possbile LSs. If the
#                 enumeration has been done previously, file LEWIS_Nc_Ne.dat
#                 will be read in.
#                 NOTE: 1) Indices start from 1 (not 0)
#                       2) The list is from ordering of projections, determined
#                          by calling the wfrt_hmo() function prior to calling 
#                          this function
#                 (2) Straightfoward specification given by LP[] and BD[],
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[] and BD[], atomic indices start from 1 
#                       (not 0), and they correspond to ACTUAL indices of 
#                       atoms, which may not be consecutive, not as those in 
#                       LP0[] and BD0[] (also see notes in body of this func)
#   sciform:      whether showing numerical results in scientific format
def wfrt_hmo_spec_spin( hmoSol, ixLS, sciform=False ):
    LP, BD, RA, RB, isRadical = ixLS

    NAt = hmoSol.NAt
    NEa = hmoSol.NEa # Number of alpha electrons
    NEb = hmoSol.NEb # Number of alpha electrons
    aoIx = []
    for i in range( 1, NAt+1 ):
        aoIx.append( np.array([ i ]) )

    CAONAO = np.eye( NAt )
    CLMO = ( hmoSol.C[ :, 0:NEa ], hmoSol.C[ :, 0:NEb ] )
    ELMO = ( hmoSol.E[ 0:NEa ], hmoSol.E[ 0:NEb ] )
    atIx1 = range( 1, NAt+1 )
    moIx1 = ( range( 1, NEa+1 ), range( 1, NEb+1 ) )

    hmoSol.C = ( hmoSol.C, hmoSol.C )
    hmoSol.E = ( hmoSol.E, hmoSol.E )
    wfrt_spec_spin( aoIx, hmoSol, CAONAO, CLMO, ELMO, atIx1, moIx1, \
            ( LP, BD, RA, RB ), sciform )

    return
# enddef wfrt_hmo_spec_spin()



# Energy Minimization Resonance Theory (EMRT) calculation with the simple 
# Hueckel approximation with a specified set of Lewis structures:
#   hmoSol:      Solution of simple HMO theory
#   ixLS:         Lewis structures specified in the following way:
#                 (1) Straightfoward specification given by LP[], BD[], RA[],
#                 RB[] (Lone-Pairs, BonDs, Radical-Alpha & Radical-Beta sites)
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[], BD[] RA[], RB[], atomic indices start 
#                       from 1 (not 0), and they correspond to ACTUAL indices 
#                       of atoms, which may not be consecutive, not as those in
#                       LP0[] and BD0[]
def emrt_spec( hmoSol, ixLS ):
    NAt = hmoSol.NAt
    atIx = range( 1, NAt+1 )
    atIx = np.array( atIx, dtype='i' )

    ifOpenShell = ( len(ixLS) == 4 )
    if ifOpenShell:
        LP, BD, RA, RB = ixLS
    else:
        LP, BD = ixLS
    # Set sequential indices of appearance of atoms:
    #  Suffix 0 means the indices of atoms stored in LP0[] and BD0[] follow the
    #  order of appearance, i.e., 1, 2, 3, ..., which is NOT necessarily the
    #  actural indices of atoms appeared in the given molecule
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

    if ifOpenShell:
        RA0 = list( map(lambda x: sqIx[x-1]+1, RA ) )
        RB0 = list( map(lambda x: sqIx[x-1]+1, RB ) )
        NEa = len(LP[0]) + len(BD[0]) + len( RA[0] ) # Num. of alpha electrons
        NEb = len(LP[0]) + len(BD[0]) + len( RB[0] ) # Num. of beta electrons
    else:
        NEa = len(LP[0]) + len(BD[0])
        NEb = NEa
    NE = NEa + NEb # Total number of electrons
    print( '%i alpha electron(s), %i beta electron(s), %i electron(s) in total'
            % ( NEa, NEb, NE ) )

    aoIx = []
    for i in range( 1, NAt+1 ):
        aoIx.append( np.array([ i ]) )

    NLS = len( LP )  # Total number of Lewis structures
    print( 'A total of %i Lewis structures to be considered' % NLS )

    CAONAO = np.eye( NAt )
    CRAO = np.identity( NAt )
    # Fock matrix in the AO basis:
    FAO = -hmoSol.A # Setting alpha=0, beta=-1
    #print( FAO )
    #input()

    # Lewis unit orbitals (LUOs) in AO basis (NAOs, Lowdin's AOs etc.):
    atIx_LUO, C_LUO, S_LUO, F_LUO = luo( CRAO, FAO, atIx, NAt, NE )
    N_LUO = len( S_LUO )
    # Print the RAO indices in each LUO:
    print( '-'*25 )
    print( '%5s   %4s   %s' % ('LUO', 'Type', 'RAO(s)') )
    print( '-'*25 )
    k = 1
    for ix in atIx_LUO:
        print( '%5i' % k, end='' )
        if ix[1] == -1: # LP
            if ifOpenShell:
                print( '  (AT):  %3i' % atIx[ix[0]] )
            else:
                print( '  (LP):  %3i' % atIx[ix[0]] )
        else:           # BD
            print( '  (BD):  %3i %3i' % (atIx[ix[0]],atIx[ix[1]]) )
        k += 1
    print( '-'*25 )

    # Matrix storing the indices of LUOs of a given set of Lewis structures:
    if ifOpenShell:
        luoIXa, luoIXb = luoIxLewis_open( atIx_LUO, NEa, NEb, \
                LP0, BD0, RA0, RB0 )
    else:
        luoIX = luoIxLewis( atIx_LUO, NE, LP0, BD0 )

    # Overlap matrix between determinants of all considered Lewis structures:
    if ifOpenShell:
        S, F, rank = overlapFockLewis_open( NEa, NEb, luoIXa, luoIXb, \
                S_LUO, F_LUO )
    else:
        S, F, rank = overlapFockLewis( NE, luoIX, S_LUO, F_LUO )

    #------------------------------------------------------------------
    # Solve linear equation for coeff. using ridge linear regression:
    #------------------------------------------------------------------
    print( 'Solving secular equation for EMRT coefficients ...' )

    rank_defi = NLS - rank
    if rank_defi > 0:
        print( 'Warning: The input Lewis strucutres are not all linearly'
                ' independent\nRank deficiency is %i' % rank_defi )

    # Lowdin's symmetric orthogonalization:
    Sd, U = np.linalg.eigh( S )
    Q = U @ np.diag( 1/np.sqrt(Sd) ) @ U.T
    F_ = Q.T @ F @ Q
    E, V = np.linalg.eigh( F_ )
    Emin = np.min( E ) # Lowest energy as the solution
    ixmin = np.argmin( E )
    C_ = V[ :, ixmin ]
    C = Q @ C_ # Coefficients for the original Lewis structure basis
    W_Lowd = np.square( C_ )
    print( Emin )
    print( C )
    print( W_Lowd )


    return
# enddef emrt_spec()



# WFRT analysis in the simple Hueckel Molecular Orbital (HMO) framework using
# only all possible Kekule structures:
#   hmoSol:      Solution of simple HMO theory
#   kekFileName:   Name of *.kek file where Kekule structures are stored
#   sciform:      whether showing numerical results in scientific format
def wfrt_hmo_kekule( hmoSol, kekFileName, sciform=False ):
    # Read Kekule structures from external file:
    print( 'Reading Kekule structures from file %s ...' % kekFileName )
    LP, BD, RA, RB, isRadical = readKekule( kekFileName )

    # Print spin info:
    if isRadical:
        print( 'An open-shell wave function with spin multiplicity of %i'
                % hmoSol.mult )
    else:
        print( 'A closed-shell wave function' )

    # Perform WFRT analysis:
    # 1. Closed-shell system:
    if not isRadical:
        wfrt_hmo_spec( hmoSol, ( LP, BD ), sciform )
    # 2. Spin-unrestricted:
    else:
        wfrt_hmo_spec_spin( hmoSol, ( LP, BD, RA, RB, isRadical ), sciform )
    return
# enddef wfrt_hmo_kekule()


# Calculate the projections of Lewis structures using both DMRT and WFRT, in
# the simple Hueckel Molecular Orbital (HMO) framework:
#   hmoSol:      Solution of simple HMO theory
#   ixLS:         Lewis structures specified in either of the following ways:
#                 (1) List of indices of the specified Lewis structures
#                 This choice involves enumeration of all possbile LSs. If the
#                 enumeration has been done previously, file LEWIS_Nc_Ne.dat
#                 will be read in.
#                 NOTE: 1) Indices start from 1 (not 0)
#                       2) The list is from ordering of projections, determined
#                          by calling the wfrt() function prior to calling this
#                          function
#                 (2) Straightfoward specification given by LP[] and BD[],
#                 which are passed as a packed tuple and can be obtained by
#                 parsing a sequence of strings of LSs
#                 NOTE: In lists LP[] and BD[], atomic indices start from 1 
#                       (not 0), and they correspond to ACTUAL indices of 
#                       atoms, which may not be consecutive, not as those in 
#                       LP0[] and BD0[] (also see notes in body of this func)
#                 (3) If ixLS == -1 (default), no LSs are specified
def proj_DM_WF_hmo( hmoSol, ixLS=-1 ):
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

    proj_DM_WF( aoIx, hmoSol, CAONAO, CLMO, ELMO, atIx1, moIx1, ixLS )
    return
# enddef proj_DM_WF_hmo()



# Calculate DMRT and WFRT projections only with all possible Kekule structures
# in the simple Hueckel Molecular Orbital (HMO) framework:
#   hmoSol:      Solution of simple HMO theory
#   kekFileName: Name of *.kek file where Kekule structures are stored
#   ndiv:        Divide file in n parts lest a huge file runs out of memory
def proj_DM_WF_kekule_hmo( hmoSol, kekFileName, ndiv=1 ):
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

    print( '--In the framework of simple Hueckel Molecular Orbital theory--' )

    proj_DM_WF_kekule( aoIx, hmoSol, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                       kekFileName, ndiv )
    return
# enddef proj_DM_WF_kekule_hmo()
