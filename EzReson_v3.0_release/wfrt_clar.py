# Wave-Function based Resonance Theory (WFRT) for Clar resonators
#
#  Created on Jun 22, 2021
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#  Sep 1:
#  - Refined the outputs of the CUOs: Only real C=C bonds are shown in the
#  BD-CUOs list, and only the C=C bonds or hexagonal rings involved in all Clar
#  resonators are listed; The ring indices (as in the 'Clar resonator' column
#  in the last table and as in the *.clar file) are also indicated for
#  sextet-CUOs.
#
#  Aug 31:
#  - Added outputs of the energies of CUOs and Clar resonators, which would be
#  useful for energy partioning in Clar resonators
#
#

from wfrt import *
from clar import *

# Clar unit orbitals (CUOs) in AO basis (NAOs, Lowdin's AOs etc.):
#   CUOs are the counterpart of the CUOs, but for Clar resonators
#   CUOs are composed of RAOs and have five types:
#     1) LP-CUOs, corresponding to a lone pair localized at an atom, equal
#        to the RAO of that atom
#     2) BD-CUOs, corresponding to a pure covalent bond between two atoms,
#        being a linear combination (50-50) of the RAOs of both atoms
#     3) Circular-Sextet-CUOs, corresponding to the lowest-energy pi MO of
#        benzene, without node in the MO
#     4) Vertical-Sextet-CUOs, corresponding to one the two degenerate HOMOs 
#        of benzene, with a vertical node in the MO
#     5) Horizontal-Sextet-CUOs, corresponding to one the two degenerate HOMOs 
#        of benzene, with a horizontal node in the MO
#
#  CRAO: Resonance atomic orbitals in orthonormal AO basis (NAOs, Lowdin's AOs)
#  FNAO:    Fock matrix in NAO basis
#  atIx: Indices of all atoms involved in the resonance analysis (start from 1)
#  NAt:  Number of atoms/centers in the subsystem, to which WFRT is applied
#  RG:   Indices of atoms in all rings
#  BD: List of numpy arrays storing the bonds for all Clar resonators
#  SX: List of numpy arrays storing the sextets for all Clar resonators
def cuo( CRAO, FNAO, atIx, NAt, RG, BD, SX ):
    # NOTE: Currently, we only consider covalent Clar resonators, 
    # i.e., excluding all ionic forms, but including all biradical ones

    # Get the adjacent matrix:
    A = np.zeros( ( NAt, NAt ), dtype='i' ) 
    valNum = np.zeros( ( NAt, 1 ), dtype='i' )
    for i1 in range( NAt ):
        # Search in the sextets:
        flag = False
        for sx in SX:
            for r in sx:
                rg = RG[r -1] # NOTE: indices in SX start from 1 (not 0)
                for k in range( len(rg) ):
                    if rg[k] != i1 +1:
                        continue
                    
                    if k != len(rg)-1:
                        i2 = rg[k+1] -1
                    else:
                        i2 = rg[0] -1
                    A[i1][i2] = 1
                    A[i2][i1] = 1
                    
                    if k != 0:
                        i2 = rg[k-1] -1
                    else:
                        i2 = rg[-1] -1

                    if A[i1][i2] == 1:
                        continue

                    A[i1][i2] = 1
                    A[i2][i1] = 1

                    valNum[i1] += 2
                    valNum[i2] += 2

                    # Maximum: trivalent graph:
                    if valNum[i1] >= 3:
                        flag = True
                        break
                if flag:
                    break
            if flag:
                break

        # Search in the bonds:
        flag = False
        for bd in BD:
            for b in bd:
                if b[0] != (i1 +1) and b[1] != (i1 +1):
                    continue
                if b[0] == (i1 +1):
                    i2 = b[1] -1
                else:
                    i2 = b[0] -1

                if A[i1][i2] == 1:
                    continue

                A[i1][i2] = 1
                A[i2][i1] = 1
                valNum[i1] += 1
                valNum[i2] += 1
                # Maximum: trivalent graph:
                if valNum[i1] >= 3:
                    flag = True
                    break
            if flag:
                break

    #print( A )
    #input()
    
    # Get the bonds; NOTE: atomic indices start from 0 (not 1)
    BONDS = []
    for i1 in range( NAt ):
        for i2 in range( i1+1, NAt ):
            if A[i1][i2] == 1:
                BONDS.append( [ i1, i2 ] )

    NB = len( BONDS ) # Number of bonds
    #print( NB )
    #print( BONDS )
    #input()


    # Choose hexagonal rings only for sextets:
    RG6 = []
    for rg in RG:
        # Only consider hexagons for sextets:
        if len( rg ) == 6:
           RG6.append( rg ) 
    NR6 = len( RG6 ) # Number of rings
    N_CUO = NB + 3*NR6
    print( '%i Clar unit orbitals (CUOs) in total' % N_CUO )

    # Numeration of all possible LPs and BDs associated with the CUOs
    #atIx_CUO = np.zeros( ( N_CUO, 2 ), dtype='i' ) # Inidices start from ZERO
    atIx_CUO = [] # Inidices start from ZERO
    # 2c BDs go first:
#    for i in range( NAt ):
#        for j in range( i+1, NAt ):
#            atIx_CUO.append( [ i, j ] )
    for b in BONDS:
        atIx_CUO.append( b )
    # 6c CSs (Cicular-Sextets) go second:
    for rg in RG6:
        atIx_CUO.append( [ rg[i]-1 for i in range(len(rg)) ] ) 
    # 6c VSs (Vertical-Sextets) go second:
    for rg in RG6:
        atIx_CUO.append( [ rg[i]-1 for i in [0,1,3,4] ] ) 
    # 6c HSs (Horizontal-Sextets) go second:
    for rg in RG6:
        ix = [ rg[i]-1 for i in range(len(rg)) ]
        ix.append( 0 ) # Add a terminal 0 to indicate that this is of HS type
        #atIx_CUO.append( [ [ rg[i]-1 for i in range(len(rg)) ] ] ) 
        atIx_CUO.append( ix )


    # Coefficient matrix of CUOs in terms of RAO, i.e., RAO -> CUO trans. mat.
    inv_sqrt2 = 1. / np.sqrt(2)
    inv_sqrt6 = 1. / np.sqrt(6)
    inv_sqrt12 = 1. / np.sqrt(12)
    dbl_inv_sqrt12 = 2. / np.sqrt(12)
    Q_CUO = np.zeros( ( NAt, N_CUO ) )
    for i in range( N_CUO ):
        if i < NB: # BD
            Q_CUO[ atIx_CUO[i], i ] = inv_sqrt2
        elif i < NB + NR6: # CS
            Q_CUO[ atIx_CUO[i], i ] = inv_sqrt6
        elif i < NB + 2*NR6: # VS
            ix = atIx_CUO[i]
            Q_CUO[ ix[0], i ] =  0.5
            Q_CUO[ ix[1], i ] =  0.5
            Q_CUO[ ix[2], i ] = -0.5
            Q_CUO[ ix[3], i ] = -0.5
        else: # HS
            ix = atIx_CUO[i]
            Q_CUO[ ix[0], i ] =  inv_sqrt12
            Q_CUO[ ix[1], i ] = -inv_sqrt12
            Q_CUO[ ix[2], i ] = -dbl_inv_sqrt12
            Q_CUO[ ix[3], i ] = -inv_sqrt12
            Q_CUO[ ix[4], i ] =  inv_sqrt12
            Q_CUO[ ix[5], i ] =  dbl_inv_sqrt12
    #writeMat( Q_CUO )

    # Overlap matrix between CUOs:
    S_CUO = Q_CUO.T @ Q_CUO  # Note that RAOs are orthonormal to each other
    #writeMat( S_CUO )

    # LUO coefficients in NAO basis only for involved atoms, atIx[]:
    C_RAO = CRAO[ :, atIx-1 ] # NOTE: Indices in atIx[] start from 1 (not 0)
    C_CUO = C_RAO @ Q_CUO

    # Fock matrix between CUOs:
    F_CUO = C_CUO.T @ FNAO @ C_CUO
    #writeMat( F_CUO )

    return ( atIx_CUO, C_CUO, S_CUO, F_CUO, NR6 )
# enddef cuo()



# Determine the indices of CUOs of a given Clar resonator:
#   lp0:   Indices of atoms in the LPs for a given Clar resonator
#   bd0:   Indices of atoms in the BDs for a given Clar resonator
def cuoIx_one_Clar( atIx_CUO, lp0, bd0, sx, RG, NR6 ):
    N_CUO = len( atIx_CUO )
    # Identify BD in terms of CUOs:
    nBD = len( bd0 ) # Number of BDs in this Clar resonator
    cuoIx_BD = []
    for i in range( nBD ):
        # NOTE: Indices in bdO[] start from 1 (not 0)
        # np.array->list to allow comparison with atIx_CUO[k](of list type too)
        ix_BD = ( bd0[i] - 1 ).tolist()
        for k in range( N_CUO ):
            if atIx_CUO[k] == ix_BD:
                ix_cuo = k
                break
        cuoIx_BD.append( ix_cuo )
    #print( cuoIx_BD )
           
    # Identify LP in terms of CUOs:
    nLP = len( lp0 ) # Number of LPs in this Clar resonator
    cuoIx_LP = []
    for i in range( nLP ):
        ix_LP = [ lp0[i] - 1, -1 ] # NOTE: Indices in LPO[] start from 1 (not 0)
        ix_cuo = np.nonzero( np.all( atIx_CUO-ix_LP==0, axis=1 ) )[0]
        # Since ix_cuo is a one-element array, put [0] to extract this number
        cuoIx_LP.append( ix_cuo[0] )
    #print( cuoIx_LP )

    # Identify CS (Circular-Sextet) in terms of CUOs:
    nSX = len( sx ) # Number of sextets in this Clar resonator
    cuoIx_CS = []
    cuoIx_VS = []
    cuoIx_HS = []
    for i in range( nSX ):
        # NOTE: Indices in rg[] start from 1 (not 0)
        #  np.array->list to allow comparison with atIx_CUO[k](of list type too)
        #  The sextet indices in sx[1] start from 1 (not 0)
        ix_SX = ( np.array(RG[ sx[i]-1 ], dtype='i') - 1 ).tolist()
        for k in range( N_CUO ):
            if len( atIx_CUO[k] ) == 6:
                if atIx_CUO[k] == ix_SX:
                    ix_cuo = k
                    break
        cuoIx_CS.append( ix_cuo )
        cuoIx_VS.append( ix_cuo + NR6 )
        cuoIx_HS.append( ix_cuo + 2*NR6 )
    #print( cuoIx_CS )
    #print( cuoIx_VS )
    #print( cuoIx_HS )


    # Final indices of CUOs for all BDs, LPs, CSs, VSs, HSs in this Clar reson.
    if nBD == 0 and nSX == 0:   # Only LPs 
        cuoIx = np.array( cuoIx_LP, dtype='i')
    elif nLP == 0 and nSX == 0: # Only BDs
        cuoIx = np.array( cuoIx_BD, dtype='i')
    elif nLP == 0 and nBD == 0: # Only SXs
        cuoIx = np.hstack( ( np.array( cuoIx_CS, dtype='i' ), \
                             np.array( cuoIx_VS ,dtype='i' ), \
                             np.array( cuoIx_HS ,dtype='i' ) ) )
    elif nLP == 0 and nBD == 0: # Only LPs and BDs
        cuoIx = np.hstack( ( np.array( cuoIx_BD, dtype='i' ), \
                             np.array( cuoIx_LP ,dtype='i' ) ) )
    elif nBD == 0: # Only LPs and SXs
        cuoIx = np.hstack( ( np.array( cuoIx_CS, dtype='i' ), \
                             np.array( cuoIx_VS ,dtype='i' ), \
                             np.array( cuoIx_HS ,dtype='i' ), \
                             np.array( cuoIx_LP ,dtype='i' ) ) )
    elif nLP == 0: # Only BDs and SXs
        cuoIx = np.hstack( ( np.array( cuoIx_CS, dtype='i' ), \
                             np.array( cuoIx_VS ,dtype='i' ), \
                             np.array( cuoIx_HS ,dtype='i' ), \
                             np.array( cuoIx_BD ,dtype='i' ) ) )
    else:          # SXs, BDs and LPs are all present
        cuoIx = np.hstack( ( np.array( cuoIx_CS, dtype='i' ), \
                             np.array( cuoIx_VS ,dtype='i' ), \
                             np.array( cuoIx_HS ,dtype='i' ), \
                             np.array( cuoIx_BD ,dtype='i' ), \
                             np.array( cuoIx_LP ,dtype='i' ) ) )
    #print( cuoIx )
    return cuoIx
# enddef cuoIx_one_Clar



# Get matrix storing the indices of CUOs of a given set of Clar resonators:
#   atIx_CUO: Natural inidices of atoms in each CUO, being 0, 1, 2, ...
#   NE: Number of electrons
#   LP0: List of LPs for all considered Clar resonators
#   BD0: List of BDs for all considered Clar resonators
# Return:
#   cuoIX: a matrix, each row corresponding to a Clar resonator;
#          and the columns are the inidices of CUOs consitituting the Clar
#          resonators.
def cuoIxClar( atIx_CUO, NE, LP0, BD0, SX, RG, NR6 ):
    NClar = len( LP0 ) # Number of Clar resonators
    if NE % 2 == 1:
        raise NotImplementedError( 'Number of electrons must be even' )
    NP = NE // 2 # Number of electron pairs
    
    # Initialize the lists of CUO indices for all considered Clar resonators:
    cuoIX = np.ones( ( NClar, NP ), dtype='i' ) * -1 
    for i in range( NClar ):
        cuoIX[ i, : ] = cuoIx_one_Clar( atIx_CUO, LP0[i], BD0[i], \
                SX[i], RG, NR6 )
    
    return cuoIX
# enddef cuoIxClar()



# Calculate the projections of Clar resonators using WFRT:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   CLMO:         LMO coefficients in NAO basis
#   ELMO:         LMO energies
#   atIx1:        Indices of atoms in the resonance subsystem, start from 1
#   moIx1:        Indices of LMOs in the resonance subsystem, start from 1
#   ixClar:       Clar resonators specified in either of the following ways:
#   clarstyle:    'short' for printing ring indices; 'explicit' for printing
#                 atom indices for all rings
#   sciform:      whether showing numerical results in scientific format
#   raoType:      'uni'--> universal RAOs;  'var'--> varying RAOs
#   inpBaseName:  Write RAOs into a *.fchk file; Default: '' --> Do not write
def do_proj_WF_clar( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                ixClar, clarstyle='short', sciform=False, raoType='uni', \
                inpBaseName='', raoFlipAtoms=None ):
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

    #-------------------- Clar resonators ----------------------------------
    #---------- Straightforward specification ----------
    LP, BD, SX, RG = ixClar

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

    # Only choose the specified Clar resonators:
    NClar = len( LP )  # Total number of Clar resonators
    
    print( 'Using %i straightforwardly specified Clar resonators' % NClar )


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

    # Clar unit orbitals (CUOs) in AO basis (NAOs, Lowdin's AOs etc.):
 
    #print( '-'*80 )
    #writeMat( CRAO, '/dev/stdout', ' %12.6f' )
    #input()


    # Density matrices of Clar resonators in NAO (or Lowdin etc.) basis:
    #print( 'Calculating densitry matrices for the Lewis strctures ...' )
    #D = [] 
    #for k in range( NLS ):
    #    D.append( DM_RAO_uni( CRAO, LP[k], BD[k] ) )
     
    #writeMat( D[0], '/dev/stdout', ' %12.6f' )
    #print( '-'*60 )
    #writeMat( D[1], '/dev/stdout', ' %12.6f' )
    #input()
     
    ## Energies of Lewis structures:
    #print( 'Calculating energies for the Lewis strctures ...' )
    #E = [] 
    #for k in range( NLS ):
    #    E.append( np.trace( D[k] @ FNAO ) )
    #E = np.c_[ E ]


    #------------------------- WFRT projections -------------------------#
    print( '='*80 )
    print( 'Performing Wave Function based Resonance Theory (WFRT) '
            'calculations ...' )
    print( '='*80 )
    # Clar unit orbitals (CUOs) in AO basis (NAOs, Lowdin's AOs etc.):
    atIx_CUO, C_CUO, S_CUO, F_CUO, NR6 = \
            cuo( CRAO, FNAO, atIx, NAt, RG, BD, SX )
    N_CUO = len( S_CUO )

    # Print the RAO indices in each CUO:
    print( '-'*72 )
    if ifHMO:
        print( '%5s %4s %4s  %15s  %12s   %s' % 
                ('CUO', 'Type', 'Ring', 'Energy(|beta|)', 'RE(|beta|)', \
                        'RAOs') )
    else:
        print( '%5s %4s %4s  %15s  %12s   %s' % 
                ('CUO', 'Type', 'Ring', 'Energy(a.u.)', 'RE(kcal/mol)', \
                        'RAOs') )
    print( '-'*72 )
    Etot_LMO_per_ele = Etot_LMO / NE
    E_SX1 = []
    E_SX2 = []
    E_SX3 = []
    Erel_SX1 = []
    Erel_SX2 = []
    Erel_SX3 = []
    k = 0
    for ix in atIx_CUO:
        if ifHMO:
            Erel_CUO = F_CUO[k][k]  - Etot_LMO_per_ele
        else:
            Erel_CUO = au2kcal( F_CUO[k][k]  - Etot_LMO_per_ele )
        print( '%5i' % (k+1), end='' )  
        if len( ix ) == 2: # BD        
            print( '   BD %4s  %15.8f  %12.6f  %3i %3i' % 
                    ('--', F_CUO[k][k], Erel_CUO, atIx[ix[0]], atIx[ix[1]]) )
        elif len( ix ) == 6: # CS
            # Get the ring index:
            rgix = -1
            for r in range( len(RG) ):
                if len( RG[r] ) != 6:
                    continue
                rg = [ atIx[ix[0]], atIx[ix[1]], atIx[ix[2]], \
                    atIx[ix[3]], atIx[ix[4]], atIx[ix[5]] ]
                if rg == RG[r]:
                    rgix = r +1
                    break
            print( '   CS %4i  %15.8f  %12.6f  %3i %3i %3i %3i %3i %3i' % 
                    (rgix, F_CUO[k][k], Erel_CUO, atIx[ix[0]], atIx[ix[1]], \
                    atIx[ix[2]], atIx[ix[3]], atIx[ix[4]], atIx[ix[5]]) )
            E_SX1.append( F_CUO[k][k] )
            Erel_SX1.append( Erel_CUO )
        elif len( ix ) == 4: # VS
            # Get the ring index:
            rgix = -1
            for r in range( len(RG) ):
                if len( RG[r] ) != 6:
                    continue
                rg = [ atIx[ix[0]], atIx[ix[1]], atIx[ix[2]], \
                    atIx[ix[3]] ]
                x = np.array(RG[r])
                if np.array_equal( rg, np.array(RG[r])[[0,1,3,4]] ):
                    rgix = r +1
                    break
            print( '   VS %4i  %15.8f  %12.6f  %3i %3i %3i %3i' % 
                    (rgix, F_CUO[k][k], Erel_CUO, atIx[ix[0]], atIx[ix[1]], \
                    atIx[ix[2]], atIx[ix[3]]) )
            E_SX2.append( F_CUO[k][k] )
            Erel_SX2.append( Erel_CUO )
        elif len( ix ) == 7: # HS
            # Get the ring index:
            rgix = -1
            for r in range( len(RG) ):
                if len( RG[r] ) != 6:
                    continue
                rg = [ atIx[ix[0]], atIx[ix[1]], atIx[ix[2]], \
                    atIx[ix[3]], atIx[ix[4]], atIx[ix[5]] ]
                if rg == RG[r]:
                    rgix = r +1
                    break
            print( '   HS %4i  %15.8f  %12.6f  %3i %3i %3i %3i %3i %3i' % 
                    (rgix, F_CUO[k][k], Erel_CUO, atIx[ix[0]], atIx[ix[1]], \
                    atIx[ix[2]], atIx[ix[3]], atIx[ix[4]], atIx[ix[5]])  )
            E_SX3.append( F_CUO[k][k] )
            Erel_SX3.append( Erel_CUO )
        k += 1
    print( '-'*72 )

    # One-electron energy per electron (in a.u.) for the sextets:
    E_SX = ( np.array(E_SX1) + np.array(E_SX2) + np.array(E_SX3) ) / 3
    #print( E_SX )
    # Relative one-electron energy per electron (in a.u.) for the sextets:
    Erel_SX = ( np.array(Erel_SX1) + np.array(Erel_SX2) \
              + np.array(Erel_SX3) ) / 3
    #print( Erel_SX )

    # Print out the one-electron energies of all involved bonds and sextets:
    print()
    print( 'One-electron energies of all involved bonds and sextets:' )
    print( '-'*72 )
    if ifHMO:
        print( '%5s  %4s  %16s  %16s   %s' % 
                ('No.', 'Type', 'E(|beta|)', 'RE/ele(|beta|)', 'Atoms') )
    else:
        print( '%5s %4s %16s  %16s   %s' % 
                ('No.', 'Type', 'E/ele(a.u.)', 'RE/ele(kcal/mol)', 'Atoms') )
    print( '-'*72 )
    k = 0
    iSX = 0
    for ix in atIx_CUO:
        if ifHMO:
            Erel_CUO = F_CUO[k][k]  - Etot_LMO_per_ele
        else:
            Erel_CUO = au2kcal( F_CUO[k][k]  - Etot_LMO_per_ele )
        # Bonds:
        if len( ix ) == 2: # BD        
            print( '%5i' % (k+1), end='' )  
            print( ' BOND %16.9f  %16.8f  %3i %3i' % 
                    (F_CUO[k][k], Erel_CUO, atIx[ix[0]], atIx[ix[1]]) )
        # Sextets:
        elif len( ix ) == 6: # CS
            # Get the ring index:
            rgix = -1
            for r in range( len(RG) ):
                if len( RG[r] ) != 6:
                    continue
                rg = [ atIx[ix[0]], atIx[ix[1]], atIx[ix[2]], \
                    atIx[ix[3]], atIx[ix[4]], atIx[ix[5]] ]
                if rg == RG[r]:
                    rgix = r +1
                    break
            print( '%5i' % rgix, end='' )  
            print( ' SEXT %16.9f  %16.8f  %3i %3i %3i %3i %3i %3i' % 
                    (E_SX[iSX], Erel_SX[iSX], atIx[ix[0]], atIx[ix[1]], \
                    atIx[ix[2]], atIx[ix[3]], atIx[ix[4]], atIx[ix[5]]) )
            iSX += 1

        k += 1
    print( '-'*72 )
    print()

    # Matrix storing the indices of CUOs of a given set of Clar resonators:
    cuoIX = cuoIxClar( atIx_CUO, NE, LP0, BD0, SX, RG, NR6 )

    # Overlap matrix between CUOs and LMOs:
    print( 'Calculating overlap matrix between CUOs and LMOs ...' )
    S_CUO_LMO = np.zeros( ( N_CUO, NP ) )
    for i in range( N_CUO ): # Running over all CUOs
        for j in range( NP ): # Running over all LMOs
            S_CUO_LMO[ i, j ] = C_CUO[ :, i ].T @ C_LMO[ :, j ]
    #writeMat( S_CUO_LMO )

    # Projection of Clar resonator determinants onto that of Psi0:
    print( 'Calculating projections of Clar resonator wave functions onto '
           'the actual wave function ...' )
    P = np.zeros( ( NClar, 1 ) ) # Initialize the projections
    A = np.zeros( ( NP, NP ) ) # Initialize single-spin determinant
    for i in range( NClar ):
        for j in range( NP ): # Running over all electron pairs in a Clar reson
            iCUO = cuoIX[ i, j ] # Index of CUO for this electron pair
            for iLMO in range( NP ): # Index of LMO
                A[ j, iLMO ] = S_CUO_LMO[ iCUO, iLMO ]
        P[i] = np.linalg.det( A ) ** 2
    #writeMat( P, '/dev/stdout', '%12.6f' )
         
    #----------------------------------------------------------------------


    # Sort by P:
    ix_sort = np.argsort( -P[:,0] )
    P = P[ ix_sort ]

    # Update LP[], BD[], SX[] and cuoIX:
    LP = [ LP[i] for i in ix_sort ]
    BD = [ BD[i] for i in ix_sort ]
    SX = [ SX[i] for i in ix_sort ]
    cuoIX = cuoIX[ ix_sort, : ]

    # One-electron energy of Clar resonators relative to Etot_LMO:
    Erel = np.zeros( ( NClar, 1 ) )
    for i in range( NClar ):
        for ix in cuoIX[i]:
            Erel[i] += F_CUO[ix][ix]
        if ifHMO:
            Erel[i] = 2*Erel[i] - Etot_LMO
        else:
            Erel[i] = au2kcal( 2*Erel[i] - Etot_LMO )
#    print( 'Etot_LMO = ', Etot_LMO )


    # Strings for Clar resonators:
    str_clar = []
    for i in range( NClar ):
        #str_clar.append( clar_str( LP[i], BD[i], SX[i], RG, True ) )
        str_clar.append( clar_str( LP[i], BD[i], SX[i], RG, \
                clarstyle=='explicit' ) )
    #print( str_clar )


    #------------------------------------------------------------------
    # Print out the reordered list of Clar resonators:
    #------------------------------------------------------------------
    #------------------------------------------------------------------
    print( '-'*80 )
    print( '%5s' % 'No.', end='' )
    if sciform:
        print( ' %20s' % 'Proj (WF)', end='' )
        if ifHMO:
            print( ' %20s' % 'Erel (|beta|)', end='' )
        else:
            print( ' %20s' % 'Erel (kcal/mol)', end='' )
    else:
        print( ' %12s' % 'Proj (WF)', end='' )
        if ifHMO:
            print( '   %15s' % 'Erel (|beta|)', end='' )
        else:
            print( '   %15s' % 'Erel (kcal/mol)', end='' )
    print( '    Clar resonator' )
    print( '-'*80 )
    for i in range( NClar ):
        print( '%5i' % ( ix_sort[i]+1 ), end='' )
        if sciform:
            print( ' %20.12E' % P[i], end='' )
            print( ' %20.10E' % Erel[i], end='' )
        else:
            print( ' %12.10f' % P[i], end='' )
            print( '   %15.8f' % Erel[i], end='' )
        print( '    %s' % str_clar[i] )
    print( '-'*80 )
    #------------------------------------------------------------------

    if ifHMO:
        print( 'Total HF energy of reference structure is %.8E |beta|' % 
                Etot_LMO )
    else:
        print( 'Total HF energy of reference structure is %.8E a.u.' % 
                Etot_LMO )

    return
# enddef do_proj_WF_clar()



# Calculate WFRT projections only with Clar resonators, which are read from 
# *.clar file and will be automatically enumerated if the file does not exist:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   CLMO:         LMO coefficients in NAO basis
#   ELMO:         LMO energies
#   atIx1:        Indices of atoms in the resonance subsystem, start from 1
#   moIx1:        Indices of LMOs in the resonance subsystem, start from 1
#   clarFileName: Name of *.clar file where Clar resonators are stored
#   ndiv:         Divide file in n parts lest a huge file runs out of memory
#   clarstyle:    'short' for printing ring indices; 'explicit' for printing
#                 atom indices for all rings
#   sciform:      whether showing numerical results in scientific format
#   raoType:      'uni'--> universal RAOs;  'var'--> varying RAOs
#   inpBaseName:  Write RAOs into a *.fchk file; Default: '' --> Do not write
def proj_WF_clar( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                 clarFileName, ndiv=1, clarstyle='short', sciform=False, \
                 raoType='uni', inpBaseName='', raoFlipAtoms=None ):
    if ndiv > 1:
        print( 'Dividing file %s into %i parts ...' % ( clarFileName, ndiv ) )
    elif ndiv <=0 :
        raise ValueError( 'Invalid value for argument ndiv (%i) in '
                'proj_WF_clar' % ndiv )

    for i in range( 1, ndiv+1 ):
        if ndiv > 1:
            print( '########## PART %3i OUT OF %3i ##########' % ( i, ndiv ) )

        # Read Clar resonators from external file:
        print( 'Reading Clar resonators from file %s ...' % clarFileName )
        LP, BD, SX, RG = readClar( clarFileName, ( ndiv, i ) )

        # Calculate WFRT projections:
        do_proj_WF_clar( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                    ( LP, BD, SX, RG ), clarstyle, sciform, raoType, \
                    inpBaseName, raoFlipAtoms )
    return
# enddef proj_WF_clar()


# Calculate WFRT projections only with Clar resonators, which are read from 
# *.clar file and will be automatically enumerated if the file does not exist,
# in the simple Hueckel Molecular Orbital (HMO) framework:
#   hmoSol:      Solution of simple HMO theory
#   clarFileName: Name of *.clar file where Clar resonators are stored
#   ndiv:        Divide file in n parts lest a huge file runs out of memory
#   clarstyle:    'short' for printing ring indices; 'explicit' for printing
#                 atom indices for all rings
#   sciform:      whether showing numerical results in scientific format
def proj_WF_clar_hmo( hmoSol, clarFileName, ndiv=1, clarstyle='short', \
        sciform=False ):
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

    proj_WF_clar( aoIx, hmoSol, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                       clarFileName, ndiv, clarstyle, sciform )
    return
# enddef proj_WF_clar_hmo()



# Perform WFRT analysis for a given set of Clar resonators:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   CLMO:         LMO coefficients in NAO basis
#   ELMO:         LMO energies
#   atIx1:        Indices of atoms in the resonance subsystem, start from 1
#   moIx1:        Indices of LMOs in the resonance subsystem, start from 1
#   ixClar:       Clar resonators specified in either of the following ways:
#   projCut:      Cutoff value of projections relative to the principal LS
#   clarstyle:    'short' for printing ring indices; 'explicit' for printing
#                 atom indices for all rings
#   sciform:      whether showing numerical results in scientific format
#   raoType:      'uni'--> universal RAOs;  'var'--> varying RAOs
#   inpBaseName:  Write RAOs into a *.fchk file; Default: '' --> Do not write
def do_wfrt_clar( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                ixClar, projCut=0., clarstyle='short', sciform=False, \
                raoType='uni', inpBaseName='', raoFlipAtoms=None ):
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

    #-------------------- Clar resonators ----------------------------------
    #---------- Straightforward specification ----------
    LP, BD, SX, RG = ixClar

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

    # Only choose the specified Clar resonators:
    NClar = len( LP )  # Total number of Clar resonators
    
    print( 'Using %i straightforwardly specified Clar resonators' % NClar )


    # Resonance atomic orbitals (RAOs) in AO basis (NAOs, Lowdin's AOs etc.):
    if raoType == 'uni': # Universal RAOs that are independent of LSs
        CRAO = rao_uni( D0, aoIx, raoFlipAtoms )
    elif raoType == 'uni-gen': # General universal RAOs
        raise NotImplementedError( 'General universal RAOs not supported yet '
                'in WFRT' )
    else:
        raise ValueError( 
                'Unrecognized value for argument raoType in function '
                'do_wfrt_clar()' )
    
    # Write RAOs into a *.fchk file:
    if len(inpBaseName) > 0:
        raoOutpName = inpBaseName + '_RAO.fchk'
        CAORAO = CAONAO @ CRAO
        FRAO = CAORAO.T @ F @ CAORAO
        E_RAO = np.diag( FRAO )
        writeFchkOrb( inpBaseName+'.fchk', raoOutpName, CAORAO, E_RAO )
        print( 'RAOs written into ' + raoOutpName )

    # Clar unit orbitals (CUOs) in AO basis (NAOs, Lowdin's AOs etc.):
 
    #print( '-'*80 )
    #writeMat( CRAO, '/dev/stdout', ' %12.6f' )
    #input()


    # Density matrices of Clar resonators in NAO (or Lowdin etc.) basis:
    #print( 'Calculating densitry matrices for the Lewis strctures ...' )
    #D = [] 
    #for k in range( NLS ):
    #    D.append( DM_RAO_uni( CRAO, LP[k], BD[k] ) )
     
    #writeMat( D[0], '/dev/stdout', ' %12.6f' )
    #print( '-'*60 )
    #writeMat( D[1], '/dev/stdout', ' %12.6f' )
    #input()
     
    ## Energies of Lewis structures:
    #print( 'Calculating energies for the Lewis strctures ...' )
    #E = [] 
    #for k in range( NLS ):
    #    E.append( np.trace( D[k] @ FNAO ) )
    #E = np.c_[ E ]


    #------------------------- WFRT calculations -------------------------#
    print( '='*80 )
    print( 'Performing Wave Function based Resonance Theory (WFRT) '
            'calculations ...' )
    print( '='*80 )
    # Clar unit orbitals (CUOs) in AO basis (NAOs, Lowdin's AOs etc.):
    atIx_CUO, C_CUO, S_CUO, F_CUO, NR6 = \
            cuo( CRAO, FNAO, atIx, NAt, RG, BD, SX )
    N_CUO = len( S_CUO )
    # Print the RAO indices in each CUO:
    print( '-'*72 )
    if ifHMO:
        print( '%5s %4s %4s  %15s  %12s   %s' % 
                ('CUO', 'Type', 'Ring', 'Energy(|beta|)', 'RE(|beta|)', \
                        'RAOs') )
    else:
        print( '%5s %4s %4s  %15s  %12s   %s' % 
                ('CUO', 'Type', 'Ring', 'Energy(a.u.)', 'RE(kcal/mol)', \
                        'RAOs') )
    print( '-'*72 )
    Etot_LMO_per_ele = Etot_LMO / NE
    E_SX1 = []
    E_SX2 = []
    E_SX3 = []
    Erel_SX1 = []
    Erel_SX2 = []
    Erel_SX3 = []
    k = 0
    for ix in atIx_CUO:
        if ifHMO:
            Erel_CUO = F_CUO[k][k]  - Etot_LMO_per_ele
        else:
            Erel_CUO = au2kcal( F_CUO[k][k]  - Etot_LMO_per_ele )
        print( '%5i' % (k+1), end='' )  
        if len( ix ) == 2: # BD        
            print( '   BD %4s  %15.8f  %12.6f  %3i %3i' % 
                    ('--', F_CUO[k][k], Erel_CUO, atIx[ix[0]], atIx[ix[1]]) )
        elif len( ix ) == 6: # CS
            # Get the ring index:
            rgix = -1
            for r in range( len(RG) ):
                if len( RG[r] ) != 6:
                    continue
                rg = [ atIx[ix[0]], atIx[ix[1]], atIx[ix[2]], \
                    atIx[ix[3]], atIx[ix[4]], atIx[ix[5]] ]
                if rg == RG[r]:
                    rgix = r +1
                    break
            print( '   CS %4i  %15.8f  %12.6f  %3i %3i %3i %3i %3i %3i' % 
                    (rgix, F_CUO[k][k], Erel_CUO, atIx[ix[0]], atIx[ix[1]], \
                    atIx[ix[2]], atIx[ix[3]], atIx[ix[4]], atIx[ix[5]]) )
            E_SX1.append( F_CUO[k][k] )
            Erel_SX1.append( Erel_CUO )
        elif len( ix ) == 4: # VS
            # Get the ring index:
            rgix = -1
            for r in range( len(RG) ):
                if len( RG[r] ) != 6:
                    continue
                rg = [ atIx[ix[0]], atIx[ix[1]], atIx[ix[2]], \
                    atIx[ix[3]] ]
                x = np.array(RG[r])
                if np.array_equal( rg, np.array(RG[r])[[0,1,3,4]] ):
                    rgix = r +1
                    break
            print( '   VS %4i  %15.8f  %12.6f  %3i %3i %3i %3i' % 
                    (rgix, F_CUO[k][k], Erel_CUO, atIx[ix[0]], atIx[ix[1]], \
                    atIx[ix[2]], atIx[ix[3]]) )
            E_SX2.append( F_CUO[k][k] )
            Erel_SX2.append( Erel_CUO )
        elif len( ix ) == 7: # HS
            # Get the ring index:
            rgix = -1
            for r in range( len(RG) ):
                if len( RG[r] ) != 6:
                    continue
                rg = [ atIx[ix[0]], atIx[ix[1]], atIx[ix[2]], \
                    atIx[ix[3]], atIx[ix[4]], atIx[ix[5]] ]
                if rg == RG[r]:
                    rgix = r +1
                    break
            print( '   HS %4i  %15.8f  %12.6f  %3i %3i %3i %3i %3i %3i' % 
                    (rgix, F_CUO[k][k], Erel_CUO, atIx[ix[0]], atIx[ix[1]], \
                    atIx[ix[2]], atIx[ix[3]], atIx[ix[4]], atIx[ix[5]])  )
            E_SX3.append( F_CUO[k][k] )
            Erel_SX3.append( Erel_CUO )
        #else:           # BD        
        #    print( '  (BD):  %3i %3i' % (atIx[ix[0]],atIx[ix[1]]) )
        k += 1         
    print( '-'*72 )

    # One-electron energy per electron (in a.u.) for the sextets:
    E_SX = ( np.array(E_SX1) + np.array(E_SX2) + np.array(E_SX3) ) / 3
    #print( E_SX )
    # Relative one-electron energy per electron (in a.u.) for the sextets:
    Erel_SX = ( np.array(Erel_SX1) + np.array(Erel_SX2) \
              + np.array(Erel_SX3) ) / 3
    #print( Erel_SX )

    # Print out the one-electron energies of all involved bonds and sextets:
    print()
    print( 'One-electron energies of all involved bonds and sextets:' )
    print( '-'*72 )
    if ifHMO:
        print( '%5s  %4s  %16s  %16s   %s' % 
                ('No.', 'Type', 'E(|beta|)', 'RE/ele(|beta|)', 'Atoms') )
    else:
        print( '%5s %4s %16s  %16s   %s' % 
                ('No.', 'Type', 'E/ele(a.u.)', 'RE/ele(kcal/mol)', 'Atoms') )
    print( '-'*72 )
    k = 0
    iSX = 0
    for ix in atIx_CUO:
        if ifHMO:
            Erel_CUO = F_CUO[k][k]  - Etot_LMO_per_ele
        else:
            Erel_CUO = au2kcal( F_CUO[k][k]  - Etot_LMO_per_ele )
        # Bonds:
        if len( ix ) == 2: # BD        
            print( '%5i' % (k+1), end='' )  
            print( ' BOND %16.9f  %16.8f  %3i %3i' % 
                    (F_CUO[k][k], Erel_CUO, atIx[ix[0]], atIx[ix[1]]) )
        # Sextets:
        elif len( ix ) == 6: # CS
            # Get the ring index:
            rgix = -1
            for r in range( len(RG) ):
                if len( RG[r] ) != 6:
                    continue
                rg = [ atIx[ix[0]], atIx[ix[1]], atIx[ix[2]], \
                    atIx[ix[3]], atIx[ix[4]], atIx[ix[5]] ]
                if rg == RG[r]:
                    rgix = r +1
                    break
            print( '%5i' % rgix, end='' )  
            print( ' SEXT %16.9f  %16.8f  %3i %3i %3i %3i %3i %3i' % 
                    (E_SX[iSX], Erel_SX[iSX], atIx[ix[0]], atIx[ix[1]], \
                    atIx[ix[2]], atIx[ix[3]], atIx[ix[4]], atIx[ix[5]]) )
            iSX += 1

        k += 1
    print( '-'*72 )
    print()

    # Matrix storing the indices of CUOs of a given set of Clar resonators:
    cuoIX = cuoIxClar( atIx_CUO, NE, LP0, BD0, SX, RG, NR6 )
    #print( cuoIX )


    # Overlap matrix between CUOs and LMOs:
    print( 'Calculating overlap matrix between CUOs and LMOs ...' )
    S_CUO_LMO = np.zeros( ( N_CUO, NP ) )
    for i in range( N_CUO ): # Running over all CUOs
        for j in range( NP ): # Running over all LMOs
            S_CUO_LMO[ i, j ] = C_CUO[ :, i ].T @ C_LMO[ :, j ]
    #writeMat( S_CUO_LMO )

    # Projection of Clar resonator determinants onto that of Psi0:
    print( 'Calculating projections of Clar resonator wave functions onto '
           'the actual wave function ...' )
    P = np.zeros( ( NClar, 1 ) ) # Initialize the projections
    A = np.zeros( ( NP, NP ) ) # Initialize single-spin determinant
    for i in range( NClar ):
        for j in range( NP ): # Running over all electron pairs in a Clar reson
            iCUO = cuoIX[ i, j ] # Index of CUO for this electron pair
            #print( iCUO )
            for iLMO in range( NP ): # Index of LMO
                A[ j, iLMO ] = S_CUO_LMO[ iCUO, iLMO ]
        P[i] = np.linalg.det( A ) ** 2
    #writeMat( P, '/dev/stdout', '%12.6f' )

    # Sort by P:
    ix_sort = np.argsort( -P[:,0] )
    P = P[ ix_sort ]

    # Update LP[], BD[], SX[] and cuoIX:
    LP = [ LP[i] for i in ix_sort ]
    BD = [ BD[i] for i in ix_sort ]
    SX = [ SX[i] for i in ix_sort ]
    cuoIX = cuoIX[ ix_sort, : ]

    # Strings for Clar resonators:
    str_clar = []
    for i in range( NClar ):
        #str_clar.append( clar_str( LP[i], BD[i], SX[i], RG, True ) )
        str_clar.append( clar_str( LP[i], BD[i], SX[i], RG, \
                clarstyle=='explicit' ) )
    #print( str_clar )


    #------------------------------------------------------------------
    # Print out the reordered list of Clar resonators:
    #------------------------------------------------------------------
    print( '-'*80 )
    print( '%5s %15s    ' % ( 'No.', 'Projection'), end='' )
    print( ' '*NE, end='' )
    print( ' CUOs', end='' )
    print( ' '*NE, end='' )
    print( '    Clar resonator' )
    print( '-'*80 )
    for i in range( NClar ):
        if sciform:
            print( '%5i %15.8E   ' % ( ix_sort[i]+1, P[i] ), end='' )
        else:
            print( '%5i %15.10f   ' % ( ix_sort[i]+1, P[i] ), end='' )
        for ix in cuoIX[i]:
            print( ' %4i' % ( ix+1 ), end='' )
        print( '     %s' % str_clar[i] )
    print( '-'*80 )
    #------------------------------------------------------------------


    # Overlap matrix between determinants of all considered Clar resonators:
    S, rank = overlapLewis( NE, cuoIX, S_CUO )
    #writeMat( S, '/dev/stdout', ' %4.2f' )
    print( '%i linearly independent Clar resonators' % rank )


    #------------------------------------------------------------------
    # Solve linear equation for coeff. using ridge linear regression:
    #------------------------------------------------------------------
    print( 'Solving linear equation for WFRT coefficients ...' )

    # Exclude certain Lewis structures:
    #   1. Remove the anti-Rumer Lewis structures:
    #      NOTE: NO NEED to apply Rumer's rule for nonionic Clar resonators
    ix_indep = range( NClar ) 
    ix_reduc = ix_indep       
    #print( '%i Lewis structures violating Rumer\'s rule excluded' 
    #        % ( len( ix_dep ) ) )
    #   2. Remove Lewis structures having small projections, if applicable:
    ix_minor = np.nonzero( P < projCut )[0]
    print( '%i Clar resonators of small projections ( < %.2E ) excluded' 
            % ( len( ix_minor ), projCut ) )
    ix_reduc = np.setdiff1d( ix_reduc, ix_minor ) # Ix. starting from 0
    print( 'Using totally %i Clar resonators to expand the actual wave '
            'function' % ( len( ix_reduc ) ) )
    if len( ix_reduc ) == 0:
        raise ValueError( 'Aborted. Please increase the value of projCut' )


    # Reduce P[] and cuoIX[]:
    Pr = P[ ix_reduc ]
    cuoIXr = cuoIX[ ix_reduc, : ]
    #Sr = S[np.ix_( ix_reduc, ix_reduc )]
    # Overlap matrix between determinants of all considered Clar resonators:
    #Sr, rank = overlapLewis( NE, luoIXr, S_LUO )
    Sr, Fr, rank = overlapFockLewis( NE, cuoIXr, S_CUO, F_CUO )
    #writeMat( Sr, '/dev/stdout', ' %7.5f' )
    #print('*'*80)    
    #writeMat( Fr, '/dev/stdout', ' %7.4f' )
    #input()
    #print( '%i linearly independent Lewis structures' % rank )
    assert( len( ix_reduc ) == rank )


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
    Wr_PWSO = weightPWSO( Cr, Sr, Pr )
    #writeMat( Wr_Lowd, '/dev/stdout', '%12.6f' )

    # Recover the full coefficient vector by setting zeros for dependent LSs:
    C = np.zeros( NClar ) # Initialize
    W_Mull = np.zeros( NClar ) # Initialize
    W_Bick = np.zeros( NClar ) # Initialize
    W_RS = np.zeros( NClar ) # Initialize
    W_Lowd = np.zeros( NClar ) # Initialize
    W_PWSO = np.zeros( NClar ) # Initialize
    for i in range( NClar ):
        ix_r = np.nonzero( ix_reduc == i )[0]
        if len( ix_r ) > 0:
            C[i] = Cr[ ix_r ]
            W_Mull[i] = Wr_Mull[ ix_r ]
            W_Bick[i] = Wr_Bick[ ix_r ]
            W_RS[i] = Wr_RS[ ix_r ]
            W_Lowd[i] = Wr_Lowd[ ix_r ]
            W_PWSO[i] = Wr_PWSO[ ix_r ]
    #assert( abs( C.T @ S @ C - 1. ) < 1E-10 )


    # Energies of non-interacting Lewis structures:
    N_reduc = len( ix_reduc )
    Er = np.diag( Fr ).reshape( N_reduc, 1 )
    Erel_r = Er - Etot_LMO
    Er_min = np.min( Er )
    REr = Er - Er_min # Relative energies of Clar reson. with respect to Er_min
    # The Er_min gives the Pauling--Wheland resonance energy:
    if ifHMO:
        Ereson = Er_min - Etot_LMO
    else:
        Ereson = au2kcal( Er_min - Etot_LMO )
    ix_min_Clar = []
    for k in range( N_reduc ):
        if ( not ifHMO and abs( Er[k] - Er_min ) < 0.001/au2kcal(1.) ) \
            or ( ifHMO and abs( Er[k] - Er_min ) < 0.001 ):
            ix_min_Clar.append( k+1 )


    # Total energy of the expanded wave function:
    Esum = 0.
    for i in range( N_reduc ):
        Esum += Cr[i] * Cr[i] * Fr[i,i]
        for j in range( i+1, N_reduc ):
            if j == i:
                continue
            Esum += 2 * Cr[i] * Cr[j] * Fr[i,j]

    # Recover the full energy info. by setting zeros for dependent Clars:
    Erel = np.zeros( NClar ) # Initialize
    RE = np.zeros( NClar ) # Initialize
    for i in range( NClar ):
        ix_r = np.nonzero( ix_reduc == i )[0]
        if len( ix_r ) > 0:
            Erel[i] = Erel_r[ ix_r ]
            RE[i] = REr[ ix_r ]
    #------------------------------------------------------------------

    #------------------------------------------------------------------
    # Print out the final results:
    #------------------------------------------------------------------
    ## Rearrange in descending order of PWSO weights:
    ix_sort2 = np.argsort( -W_PWSO )
    # Rearrange in descending order of Lowdin weights:
    #ix_sort2 = np.argsort( -W_Lowd )
    print( '-'*105 )
    print( '%5s %15s %15s' % ( 'No.', 'Projection', 'Coefficient' ), end='' )
    print( ' %7s' % 'RE', end='' )
    print( ' %7s' % 'Mulli.', end='' )
    print( ' %7s' % 'Bickel.', end='' )
    print( ' %7s' % 'Ros-Sc.', end='' )
    print( ' %7s' % 'Lowdin', end='' )
    if sciform:
        print( ' %15s' % 'PWSO', end='' )
    else:
        print( ' %7s' % 'PWSO', end='' )
    print( '  Clar resonator' )
    print( '-'*105 )
    for iClar in range( NClar ):
        i = ix_sort2[ iClar ]
        # Only show Lewis structures in the reduced set:
        if np.count_nonzero( ix_reduc == i ) == 0:
            continue
        if sciform:
            print( '%5i %15.8E %15.7E' % (ix_sort[i]+1, P[i], C[i]), end='' )
        else:
            print( '%5i %15.12f %15.12f' % (ix_sort[i]+1, P[i],C[i]), end='' )
        if ifHMO:
            print( ' %7.3f' % ( Erel[i] ), end='' )
        else:
            print( ' %7.2f' % ( au2kcal( Erel[i] ) ), end='' )
        print( ' %6.2f%%' % ( W_Mull[i]*100 ), end='' )
        print( ' %6.2f%%' % ( W_Bick[i]*100 ), end='' )
        print( ' %6.2f%%' % ( W_RS[i]*100 ), end='' )
        print( ' %6.2f%%' % ( W_Lowd[i]*100 ), end='' )
        if sciform:
            print( ' %15.6E%%' % ( W_PWSO[i]*100 ), end='' )
        else:
            print( ' %6.2f%%' % ( W_PWSO[i]*100 ), end='' )
        print( '  %s' % str_clar[i] )
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
    else:
        print( 'Pauling--Wheland resonance energy computed by exact wave '
               'function is %.2f kcal/mol' % Ereson )
    if len( ix_min_Clar ) > 1:
        print( 'The corresponding resonance structures are', end='' )
    else:
        print( 'The corresponding resonance structure is', end='' )
    for i in ix_min_Clar:
        print( ' #%i' % i, end='' )
    print()
     
    #------------------------------------------------------------------


    # Write overlap matrix to external file:
    if globals.ifWrite_Overlap:
        # Write the sorted OVERLAP:
        overlapFile = globals.basename + '.CLAR_OVERLAP'
        print()
        print( 'Writing overlap matrix to file %s ...' % overlapFile )

        Sr_sortByW = np.zeros( ( NClar, NClar ) )
        for iClar in range( NClar ):
            i = ix_sort2[ iClar ]
            for jClar in range( NClar ):
                j = ix_sort2[ jClar ]
                Sr_sortByW[iClar][jClar] = Sr[i][j]
                np.savetxt( overlapFile, Sr_sortByW, fmt='%24.15E' )
        print( 'File %s written' % overlapFile )

    # Write Hamiltonian matrix to external file:
    if globals.ifWrite_Hamilt:
        # Write the sorted HAMILTONIAN:
        HamiltonianFile = globals.basename + '.CLAR_HAMILT'
        print()
        print( 'Writing Hamiltonian matrix to file %s ...' % HamiltonianFile )

        Fr_sortByW = np.zeros( ( NClar, NClar ) )
        for iClar in range( NClar ):
            i = ix_sort2[ iClar ]
            for jClar in range( NClar ):
                j = ix_sort2[ jClar ]
                Fr_sortByW[iClar][jClar] = Fr[i][j]
                np.savetxt( HamiltonianFile, Fr_sortByW, fmt='%24.15E' )
        print( 'File %s written' % HamiltonianFile )

    # Write the matrix of interaction energy to external file:
    # NOTE: The diagonal elements are the relative one-electron energy of
    # each resonator, while the off-diagonal elements are the corresponding
    # interaction energies between two resonators (not splitted).
    # All values are in kcal/mol
    if globals.ifWrite_Eint:
        inteneFile = globals.basename + '.CLAR_INTENE'
        print()
        print( 'Writing interaction energy matrix to file %s ...' % 
                inteneFile )

        EINTr = np.zeros( ( NClar, NClar ) )
        for iClar in range( NClar ):
            i = ix_sort2[ iClar ]
            #EINTr[iClar][iClar] = W_Mull[i] * Fr[i][i]
            ## Convert to relative energy in kcal/mol:
            #EINTr[iClar][iClar] = au2kcal( EINTr[i][i] - W_Mull[i]*Etot_LMO )
            EINTr[iClar][iClar] = au2kcal( Erel[i] )
            for jClar in range( iClar+1, NClar ):
                j = ix_sort2[ jClar ]
                EINTr[iClar][jClar] = C[i]*C[j]* \
                        ( 2*Fr[i][j] - (Fr[i][i]+Fr[j][j]) * S[i][j] )
                # Convert to kcal/mol:
                EINTr[iClar][jClar] = au2kcal( EINTr[iClar][jClar] )
                EINTr[jClar][iClar] = EINTr[iClar][jClar]
        np.savetxt( inteneFile, EINTr, fmt='%15.8f' )
        print( 'File %s written' % inteneFile )

#enddef do_wfrt_clar()



# WFRT analysis only with Clar resonators, which are read from *.clar file and
# will be automatically enumerated if the file does not exist
# automatically:
#   naoInfo:      Information from NAO's output
#   FchkInfo:     Informatoin of fchk output
#   CAONAO:       AO-->NAO transformation matrix
#   CLMO:         LMO coefficients in NAO basis
#   ELMO:         LMO energies
#   atIx1:        Indices of atoms in the resonance subsystem, start from 1
#   moIx1:        Indices of LMOs in the resonance subsystem, start from 1
#   clarFileName: Name of *.clar file where Clar resonators are stored
#   ndiv:         Divide file in n parts lest a huge file runs out of memory
#   projCut:      Cutoff value of projections relative to the principal Clar 
#   clarstyle:    'short' for printing ring indices; 'explicit' for printing
#                 atom indices for all rings
#   sciform:      whether showing numerical results in scientific format
#   raoType:      'uni'--> universal RAOs;  'var'--> varying RAOs
#   inpBaseName:  Write RAOs into a *.fchk file; Default: '' --> Do not write
def wfrt_clar( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                 clarFileName, ndiv=1, projCut=0., clarstyle='short', \
                 sciform=False, raoType='uni', inpBaseName='', \
                 raoFlipAtoms=None ):
    if ndiv > 1:
        print( 'Dividing file %s into %i parts ...' % ( clarFileName, ndiv ) )
    elif ndiv <=0 :
        raise ValueError( 'Invalid value for argument ndiv (%i) in '
                'wfrt_clar' % ndiv )

    for i in range( 1, ndiv+1 ):
        if ndiv > 1:
            print( '########## PART %3i OUT OF %3i ##########' % ( i, ndiv ) )

        # Read Clar resonators from external file:
        print( 'Reading Clar resonators from file %s ...' % clarFileName )
        LP, BD, SX, RG = readClar( clarFileName, ( ndiv, i ) )

        # WFRT analysis for Clar resonators:
        do_wfrt_clar( naoInfo, FchkInfo, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                    ( LP, BD, SX, RG ), projCut, clarstyle, sciform, raoType, \
                    inpBaseName, raoFlipAtoms )
    return   
# enddef wfrt_clar()



# WFRT analysis in the simple Hueckel Molecular Orbital (HMO) framework, for 
# Clar resonators, which are read from *.clar file and will be automatically
# enumerated if the file does not exist:
#   hmoSol:      Solution of simple HMO theory
#   clarFileName: Name of *.clar file where Clar resonators are stored
#   ndiv:        Divide file in n parts lest a huge file runs out of memory
#   projCut:      Cutoff value of projections relative to the principal Clar
def wfrt_clar_hmo( hmoSol, clarFileName, ndiv=1, projCut=0., \
        clarstyle='short' ):
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

    #print( '--In the framework of simple Hueckel Molecular Orbital theory--' )

    wfrt_clar( aoIx, hmoSol, CAONAO, CLMO, ELMO, atIx1, moIx1, \
                       clarFileName, ndiv, projCut, clarstyle )
    return   
# enddef wfrt_clar_hmo()
