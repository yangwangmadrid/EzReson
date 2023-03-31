# Enumeration of all Kekule structures of a given molecular graph
#
#  Created on Apr 23, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    Feb  17, 2022:
#    - Extended function kekule() for both closed-shell and monoradical
#    conjugated systems
#    - Extended function readKekule() for both closed-shell and monoradical
#    conjugated systems
#
#    May 4, 2021:
#    - Added function getSCENRE_penta(). The difference from getSCENRE() is
#    that it takes NEGATIVE weights for pentagons, which has been shown to be
#    more reasonable to explain the relative WFRT weights of Kekule structures
#    for a few tested cases of pentagon-containing PAHs
#
#    Feb 12, 2021:
#    - Added function getFNRE() to obtain Formal Numbers of Ring Electrons
#
#    Feb 4, 2021:
#    - Added function getSCENRE() to obtain Self-Consistent Effective Numbers 
#      of Ring Electrons (SCENRE) based solely on topology
#

from topo import *
from util import *
import copy as cp


# Iterative enumeration of all Kekule structures for closed-shell system:
def kekule_closedshell( writer, K, B, NAt, iat ):
    # Number of bonds in the list:
    NB = len( B )

    # Finished:
    if NB == 0:       
        Knew = K
        Bnew = B

        # If this is an invalid Kekulean path:
        NK = len( K )
        if NK != (NAt // 2):
            return Knew, Bnew

        #writer.write( 'Kekule: ' )
        for i in range( NK ):
            # NOTE: Output atomic indices starting from 1 (not 0):
            writer.write( ' %i %i' % ( K[i][0]+1, K[i][1]+1 ) )
        writer.write( '\n' )

        return Knew, Bnew

    # In case that atom iat is not present in the bond list B:
    for i in range( NB ):
        if B[i][0] == iat or B[i][1] == iat:
            break
        Knew, Bnew = kekule_closedshell( writer, K, B, NAt, iat+1 )
        return Knew, Bnew

    # Find all bonds involving atom iat:
    B_iat = []
    for i in range( NB ):
        if B[i][0] == iat or B[i][1] == iat:
            B_iat.append( B[i] )

    # Run over all bonds involving atom iat:
    for i in range( len(B_iat) ):
        bnd = B_iat[i]
        # Update list of Kekule structures:
        Knew = cp.copy(K)
        Knew.append( bnd )
        # Remove all bonds involving atom iat from bond list B:
        Bnew = removeBonds( B, bnd )
        # Recursive calling:
        Knew, Bnew = kekule_closedshell( writer, Knew, Bnew, NAt, iat+1 );
        if iat == 0:
            print( 'Progress %.2f%% ...' % ( i/len(B_iat)*100 ) )

    return Knew, Bnew
# enddef kekule_closedshell()


# Iterative enumeration of all Kekule structures for both closed-shell and
# monoradical conjugated-systems:
#
#  NOTE: In the following lists, all labels start from 0, not 1
#  K:   List storing pairs of atoms for the kekule structure to be enumerated
#  A:   List storing the remaining atoms during the enumeration
#  B:   List storing the remaining bonds during the enumeration
#  NAt: Number of atoms
#  iat: Label (starting from 0) of the starting atom for the enumeration 
#
def kekule_backup( writer, K, A, B, NAt, iat ):
    # Number of bonds in the list:
    NB = len( B )

    # Finished:
    if NB == 0:       
        Knew = K
        Bnew = B
        Anew = A

        # Check consistency of the algorithm:
        if len( A ) > 1:
            raise RuntimeError( '%i atoms unassigned. Something going wrong.'
                    % len( A ) )

        # If this is an invalid Kekulean path:
        NK = len( K )
        if NK != (NAt // 2):
            return Knew, Anew, Bnew

        #writer.write( 'Kekule: ' )
        # Write pairs of atoms (i.e., Kekule bonds):
        for i in range( NK ):
            # NOTE: Output atomic indices starting from 1 (not 0):
            writer.write( ' %i %i' % ( K[i][0]+1, K[i][1]+1 ) )
        # Write the single leftover atom (i.e., monoradical site):
        writer.write( ' -%i' % ( Anew[0]+1 ) )
        writer.write( '\n' )

        return Knew, Anew, Bnew

    # In case that atom iat is not present in the bond list B:
    for i in range( NB ):
        if B[i][0] == iat or B[i][1] == iat:
            break
        Knew, Anew, Bnew = kekule_closedshell( writer, K, A, B, NAt, iat+1 )
        return Knew, Anew, Bnew

    # Find all bonds involving atom iat:
    B_iat = []
    for i in range( NB ):
        if B[i][0] == iat or B[i][1] == iat:
            B_iat.append( B[i] )

    # Run over all bonds involving atom iat:
    for i in range( len(B_iat) ):
        bnd = B_iat[i]
        # Update list of Kekule structures:
        Knew = cp.copy(K)
        Knew.append( bnd )
        # Remove all bonds involving atom iat from bond list B:
        Bnew = removeBonds( B, bnd )
        # Remove from atom list A all atoms involved in bnd:
        Anew = removeAtoms_from_Bond( A, bnd )

        # Recursive calling:
        Knew, Anew, Bnew = kekule_closedshell( writer, Knew, Anew, Bnew, NAt, iat+1 );
        if iat == 0:
            print( 'Progress %.2f%% ...' % ( i/len(B_iat)*100 ) )

    return Knew, Anew, Bnew
# enddef kekule_backup()


# Iterative enumeration of all Kekule structures for both closed-shell and
# monoradical conjugated-systems:
#
#  NOTE: In the following lists, all labels start from 0, not 1
#  K:    List storing pairs of atoms for the kekule structure to be enumerated
#  A:    List storing the remaining atoms during the enumeration
#  B:    List storing the remaining bonds during the enumeration
#  NAt:  Number of atoms
#  iat:  Label (starting from 0) of the starting atom for the enumeration 
#  iRad: Label (starting from 0) of the radical site atom for the enumeration
#        If iRad is -1 (default), the system is recognized as a closed-shell
#
def kekule( writer, K, B, NAt, iat, iRad=-1 ):
    # Get radical type:
    if isinstance( iRad, list ):
        if len( iRad ) == 3:
            radType = 3 # Triradical
        elif len( iRad ) == 2:
            if iRad[0] >= 0 and iRad[1] >= 0:
                radType = 2 # Singlet biradical
            elif iRad[0] >= 0 and iRad[1] < 0:
                radType = -2 # Triplet biradical
    else:
        if iRad == -1:
            radType = 0 # Closed-shell
        else:
            radType = 1 # Monoradical

    # Number of bonds in the list:
    NB = len( B )

    # Finished:
    if NB == 0:       
        Knew = K
        Bnew = B

        # If this is an invalid Kekulean path:
        NK = len( K )
        if NK != (NAt // 2):
            return Knew, Bnew

        #writer.write( 'Kekule: ' )
        for i in range( NK ):
            # NOTE: Output atomic indices starting from 1 (not 0):
            writer.write( ' %i %i' % ( K[i][0]+1, K[i][1]+1 ) )
        # For the monoradical site (if applicable):
        if radType == 1:
            writer.write( ' -%i' % ( iRad+1 ) )
        # For the singlet biradical sites (if applicable):
        elif radType == 2:
            writer.write( ' -%i +%i' % ( iRad[0]+1, iRad[1]+1 ) )
        # For the triplet biradical sites (if applicable):
        elif radType == -2:
            writer.write( ' -%i -%i' % ( iRad[0]+1, -iRad[1]+1 ) )
        # For the triradical sites (if applicable):
        elif radType == 3:
            writer.write( ' -%i -%i +%i' % (iRad[0]+1, iRad[1]+1, iRad[2]+1) )

        writer.write( '\n' )
        return Knew, Bnew

    # In case that atom iat is not present in the bond list B:
    for i in range( NB ):
        if B[i][0] == iat or B[i][1] == iat:
            break
        Knew, Bnew = kekule( writer, K, B, NAt, iat+1, iRad )
        return Knew, Bnew

    # Find all bonds involving atom iat:
    B_iat = []
    for i in range( NB ):
        if B[i][0] == iat or B[i][1] == iat:
            B_iat.append( B[i] )

    # Run over all bonds involving atom iat:
    for i in range( len(B_iat) ):
        bnd = B_iat[i]
        # Update list of Kekule structures:
        Knew = cp.copy(K)
        Knew.append( bnd )
        # Remove all bonds involving atom iat from bond list B:
        Bnew = removeBonds( B, bnd )
        # Recursive calling:
        Knew, Bnew = kekule( writer, Knew, Bnew, NAt, iat+1, iRad );
        if iRad == -1 and iat == 0:
            print( 'Progress %.2f%% ...' % ( i/len(B_iat)*100 ) )

    return Knew, Bnew
# enddef kekule()


# Enumeration of all Kekule structures from Cartesian coordinates:
#
#  spin --> 0: restricted; 1: unrestricted
#  mult: spin-multiplicity
def enumKekule( writer, elem, xyz, spin=0, mult=1 ):

    xyz, at = xyzNonH( elem, xyz ) # For non-hydrogen atoms
    NAt = xyz.shape[0]

    # Adjacent matrix:
    M = adjMat( xyz )

    # Bond list:
    B = bonds( M )

    # Assign back the original atomic labels in B[]:
    B0 = []
    for bnd in B:
        B0.append( np.array([ at[bnd[0]], at[bnd[1]] ]) )

    ## Atom list (labels starting from 0):
    #A0 = np.arange( NAt )
    #kekule_backup( writer, [], A0, B0, NAt, 0 )

    # For a closed-shell conjugated system:
    if spin == 0:
        if NAt % 2 == 1:
            raise ValueError( 'A closed-shell wavefunction but with '
                    'an odd number of electrons' )

        print( 'This is a CLOSED-SHELL conjugated system.' )
        kekule( writer, [], B0, NAt, 0 )
    # For an open-shell biradical conjugated system:
    elif spin == 1 and mult == 1 and NAt % 2 == 0:
        print( 'Considering NONRADICAL resonators ...' )
        kekule( writer, [], B0, NAt, 0 )
        print( 'Considering SINGLET BIRADICAL resonators ...' )
        # Run over all pairs of atoms, each of which is taken as the biradical 
        # sites:
        # Forcing that iRad1 is alpha and iRad2 is beta:
        for iRad1 in range( NAt ):
            for iRad2 in range( NAt ):
                # Assure that iRad2 is different from iRad1:
                if iRad2 == iRad1:
                    continue
                # Skip if iRad1--iRad2 form a bond:
                if M[ iRad1, iRad2 ]:
                    continue
                #print( 'Radical site on atom %i, %i' % ( iRad1, iRad2 ) )
                # Remove all bonds invovling atom iRad:
                B0_iRad = removeBonds_by_Atoms( B0, [iRad1, iRad2] )
                #print( B0_iRad )
                kekule( writer, [], B0_iRad, NAt-2, 0, [iRad1, iRad2] )
    # For an open-shell triplet conjugated system:
    elif spin == 1 and mult == 3 and NAt % 2 == 0:
        print( 'Considering TRIPLET BIRADICAL resonators ...' )
        # Run over all pairs of atoms, each of which is taken as the biradical 
        # sites:
        # Forcing that both iRad1 and iRad2 are alpha:
        for iRad1 in range( NAt ):
            for iRad2 in range( iRad1+1, NAt ):
                #print( 'Radical site on atom %i, %i' % ( iRad1, iRad2 ) )
                # Remove all bonds invovling atom iRad:
                B0_iRad = removeBonds_by_Atoms( B0, [iRad1, iRad2] )
                #print( B0_iRad )
                kekule( writer, [], B0_iRad, NAt-2, 0, [iRad1, -iRad2] )
    # For an open-shell monoradical/triradical conjugated system:
    elif spin == 1 and NAt % 2 == 1:
        # 1. Monoradical models:
        print( 'Considering MONORADICAL resonators ...' )
        # Run over all atoms, each of which is taken as the radical site:
        for iRad in range( NAt ):
            #print( 'Radical site on atom %i' % ( iRad+1 ) )
            # Remove all bonds invovling atom iRad:
            B0_iRad = removeBonds_by_Atoms( B0, iRad )
            #print( B0_iRad )
            kekule( writer, [], B0_iRad, NAt-1, 0, iRad )
            print( 'Progress %.2f%% ...' % ( (iRad+1)/NAt*100 ) )

        # 2. Triradical models:
        print( 'Considering TRIRADICAL resonators ...' )
        # Run over all atoms, each of which is taken as the radical site:
        # Forcing that iRad1 and iRad2 are both alpha and iRad3 is beta:
        for iRad1 in range( NAt ):
            for iRad2 in range( iRad1+1, NAt ):
                for iRad3 in range( NAt ):
                    # Assure that iRad3 is different from iRad1 and iRad2:
                    if iRad3 == iRad1 or iRad3 == iRad2:
                        continue
                    # Skip if iRad1--iRad3 or iRad2--iRad3 form a bond:
                    if M[ iRad1, iRad3 ] or M[ iRad2, iRad3 ]:
                        continue
                    #print( 'Radical site on atom %i, %i, %i' 
                    #        % ( iRad1, iRad2, iRad3 ) )
                    # Remove all bonds invovling atom iRad:
                    B0_iRad = removeBonds_by_Atoms( B0, [iRad1, iRad2, iRad3] )
                    #print( B0_iRad )
                    kekule( writer, [], B0_iRad, NAt-3, 0, \
                            [iRad1, iRad2, iRad3] )
            print( 'Progress %.2f%% ...' % ( (iRad1+1)/NAt*100 ) )

    return
# enddef enumKekule()


# Read Kekule structures from *.kek file and convert them to LP[] and BD[]
#   kekFileName: Name of *.kek file where Kekule structures are stored
#   div:         If div = (n,m), then divide kekFileName in n parts and only
#                consider the m-th part of the whole set of Kekule structures
#  Return:
#          - ( LP, BD ) for closed-shell system
#          - ( LP, BD, RA, RB ) for open-shell system, where RA and RB are the
#           radical sites of alpha and beta spins, respectively
#
def readKekule( kekFileName, div=(1,1) ):
    # Division of the whole file:
    n, m = div  # Divide into n parts and take the m-th partition
    # Number of lines of the file:
    nline = sum( 1 for line in open( kekFileName ) )
    wid = nline // n
    istart = (m-1) * wid + 1  # NOTE: istart == 1 (not 0) when m == 1
    iend = istart + wid - 1
    # Take care of the residual lines:
    if m == n and iend < nline:
        iend = nline

    LP = []
    BD = []
    RA = []
    RB = []
    isRadical = False
    iline = 0
    with open( kekFileName ) as f:
        for line in f:
            iline += 1
            if iline < istart or iline > iend:
                continue
            bnd = []
            rad_a = [] # radical sites of alpha spin
            rad_b = [] # radical sites of beta spin
            i = 0
            for s in line.split():
                if len( rad_a ) == 0:
                    if i == 0:
                        i = 1
                        b1 = int(s)
                        # In case of monoradical:
                        if b1 < 0:
                            rad_a.append( -b1 )
                    else:
                        i = 0
                        b2 = int(s)
                        bnd.append( [ b1, b2 ] )
                else:
                    b = int(s)
                    if b < 0:
                        rad_a.append( -b )
                    else:
                        rad_b.append( b )
            BD.append( np.array( bnd, dtype='i' ) )
            LP.append( np.array( [], dtype='i' ) )
            if len( rad_a ) > 0 or len( rad_b ) > 0:  # Radical system
                RA.append( np.array( rad_a, dtype='i' ) )
                RB.append( np.array( rad_b, dtype='i' ) )
                isRadical = True
            else:
                RA.append( np.array( [], dtype='i' ) )
                RB.append( np.array( [], dtype='i' ) )
        
    return LP, BD, RA, RB, isRadical
# enddef readKekule()


# Get Formal Numbers of Ring Electrons for a given set of Kekule structures
#   RG: 2D array storing the rings of the given molecular graph
#   BD: List of numpy arrays storing the bonds for the given Kekule structures
# Return:
#   FNRE: List of Formal Numbers of Ring Electrons for the Kekule structures
def getFNRE( RG, BD ):
    NR = len( RG )

    FNRE = []
    for bnd in BD: # Run over Kekule structures
        N = np.max( bnd ) # Number of atoms
        # Get linkage matrix:
        lm = np.zeros( ( N, N ), dtype='i' )
        for b in bnd:
            lm[ b[0]-1 ][ b[1]-1 ] = 1
            lm[ b[1]-1 ][ b[0]-1 ] = 1
        # Formal Numbers of Ring Electrons:
        fnre = np.zeros( NR, dtype='i' )
        for j in range( NR ):
            rg = RG[j]
            rgsize = len( rg )
            for k in range( rgsize ):
                iat1 = rg[k] - 1
                if k == rgsize-1:
                    iat2 = rg[0] - 1
                else:
                    iat2 = rg[k+1] - 1
                if lm[iat1][iat2] == 1:
                    fnre[j] += 2
        FNRE.append( fnre ) 
    return np.array( FNRE )
# enddef getFNRE()


# Get Self-Consistent Effective Numbers of Ring Electrons (SCENRE) for a given 
# set of Kekule structures
#   RG: 2D array storing the rings of the given molecular graph
#   BD: List of numpy arrays storing the bonds for the given Kekule structures
#    n: power law for evaluation of weights of rings (n = 2 is recommended)
# Return:
#   SCENRE: List of Formal Numbers of Ring Electrons for the Kekule structures
def getSCENRE( RG, BD, n ):
    MAX_ITER = 100 # Maximum number of iterations
    TOL = 1E-8 # Tolerence in weights

    NR = len( RG ) # Number of rings
    NK = len( BD ) # Number of Kekule structures

    ringSize = [] # ring sizes

    SCENRE = []
    FNRE = [] # Formal Numbers of Ring Electrons
    for bnd in BD: # Run over Kekule structures
        N = np.max( bnd ) # Number of atoms
        # Get linkage matrix:
        lm = np.zeros( ( N, N ), dtype='i' )
        for b in bnd:
            lm[ b[0]-1 ][ b[1]-1 ] = 1
            lm[ b[1]-1 ][ b[0]-1 ] = 1
        # Formal Numbers of Ring Electrons:
        fnre = np.zeros( NR, dtype='i' )
        for j in range( NR ):
            rg = RG[j]
            rgsize = len( rg )
            ringSize.append( rgsize )
            for k in range( rgsize ):
                iat1 = rg[k] - 1
                if k == rgsize-1:
                    iat2 = rg[0] - 1
                else:
                    iat2 = rg[k+1] - 1
                if lm[iat1][iat2] == 1:
                    fnre[j] += 2
        FNRE.append( fnre )

    FNRE = np.array( FNRE, dtype='f' )
    # Self-consistently prediction of Effective Numbers of Ring Electrons:
    # Initial guess of weights of Kekule structures based on FNREs:
    p = FNRE ** n
    # Total FNRE summed over all rings for each Kekule structure:
    p_kek = p.sum( axis=1 )
    W = p_kek / p_kek.sum()

    p_ring = np.zeros( NR )
    W0 = W
    for loop in range( MAX_ITER ):
        # Get the ENRE of each ring:
        for r in range( NR ):
            #p_ring[r] = sum( ( FNRE[:,r] ** n ) * W )
            #p_ring[r] = sum( ( FNRE[:,r] ** n ) * W ) / ringSize[r]
            p_ring[r] = sum( ( (FNRE[:,r]/ringSize[r]) ** n ) * W )
        W_ring = p_ring / p_ring.sum()
        #print( W_ring )
        # Evaluate the weight of each Kekule structure:
        for k in range( NK ):
            p_kek[k] = sum( ( FNRE[k,:] ** n ) * W_ring )
        W = p_kek / p_kek.sum()
        #print( W )

        # Check convergence:
        if max( abs( W - W0 ) ) < TOL:
            print( 'Self-consistent weights converged after %i iterations' %
                    ( loop+2 ) )
            break

        W0 = W

    #print( W )
    SCENRE = np.zeros( NR )
    for r in range( NR ):
        #SCENRE[r] = sum( FNRE[:,r] * W )
        #SCENRE[r] = sum( ( FNRE[:,r] ** n ) * W )
        #SCENRE[r] = sum( ( FNRE[:,r] ** n ) * W ) / ringSize[r]
        SCENRE[r] = sum( ( (FNRE[:,r]/ringSize[r]) ** n ) * W )

    return SCENRE
# enddef getSCENRE()


# Get Self-Consistent Effective Numbers of Ring Electrons (SCENRE) for a given 
# set of Kekule structures
#   RG: 2D array storing the rings of the given molecular graph
#   BD: List of numpy arrays storing the bonds for the given Kekule structures
#    n: power law for evaluation of weights of rings (n = 2 is recommended)
# Return:
#   SCENRE: List of Formal Numbers of Ring Electrons for the Kekule structures
def getSCENRE_penta( RG, BD, n ):
    MAX_ITER = 100 # Maximum number of iterations
    TOL = 1E-8 # Tolerence in weights

    NR = len( RG ) # Number of rings
    NK = len( BD ) # Number of Kekule structures

    ringSize = [] # ring sizes

    SCENRE = []
    FNRE = [] # Formal Numbers of Ring Electrons
    for bnd in BD: # Run over Kekule structures
        N = np.max( bnd ) # Number of atoms
        # Get linkage matrix:
        lm = np.zeros( ( N, N ), dtype='i' )
        for b in bnd:
            lm[ b[0]-1 ][ b[1]-1 ] = 1
            lm[ b[1]-1 ][ b[0]-1 ] = 1
        # Formal Numbers of Ring Electrons:
        fnre = np.zeros( NR, dtype='i' )
        for j in range( NR ):
            rg = RG[j]
            rgsize = len( rg )
            ringSize.append( rgsize )
            for k in range( rgsize ):
                iat1 = rg[k] - 1
                if k == rgsize-1:
                    iat2 = rg[0] - 1
                else:
                    iat2 = rg[k+1] - 1
                if lm[iat1][iat2] == 1:
                    fnre[j] += 2
        FNRE.append( fnre )

    FNRE = np.array( FNRE, dtype='f' )

    ringSign = np.zeros( NR )
    for r in range( NR ):
        if ringSize[r] == 6:
            ringSign[r] = 1
        elif ringSize[r] == 5:
            ringSign[r] = -1
        else:
            raise ValueError( 'Invalid ring size (%i) for ring %i'
                    % ( len( RG[r] ), r ) )
    #print( ringSign )

    # Self-consistently prediction of Effective Numbers of Ring Electrons:
    # Initial guess of weights of Kekule structures based on FNREs:
    p = FNRE ** n * ringSign
    #print( p )
    #input()
    # Total FNRE summed over all rings for each Kekule structure:
    p_kek = p.sum( axis=1 )
    #print( p_kek )
    #input()
    W = p_kek / p_kek.sum()

    p_ring = np.zeros( NR )
    W0 = W
    for loop in range( MAX_ITER ):
        # Get the ENRE of each ring:
        for r in range( NR ):
            #p_ring[r] = sum( ( FNRE[:,r] ** n ) * W )
            #p_ring[r] = sum( ( FNRE[:,r] ** n ) * W ) / ringSize[r]
            #p_ring[r] = sum( ( (FNRE[:,r]/ringSize[r]) ** n ) * W )
            p_ring[r] = sum( ( FNRE[:,r] ** n * ringSign[r] ) * W )
        #W_ring = p_ring / p_ring.sum()
        W_ring = np.abs( p_ring ) / p_ring.sum()
        # Evaluate the weight of each Kekule structure:
        for k in range( NK ):
            p_kek[k] = sum( ( FNRE[k,:] ** n * ringSign ) * W_ring )
        #print( p_kek )
        #input()
        W = p_kek / p_kek.sum()
        #print( W )

        # Check convergence:
        if max( abs( W - W0 ) ) < TOL:
            print( 'Self-consistent weights converged after %i iterations' %
                    ( loop+2 ) )
            break

        W0 = W

    #print( W )
    SCENRE = np.zeros( NR )
    for r in range( NR ):
        #SCENRE[r] = sum( FNRE[:,r] * W )
        #SCENRE[r] = sum( ( FNRE[:,r] ** n ) * W )
        #SCENRE[r] = sum( ( FNRE[:,r] ** n ) * W ) / ringSize[r]
        SCENRE[r] = sum( ( (FNRE[:,r]/ringSize[r]) ** n ) * W )

    return SCENRE
# enddef getSCENRE_penta()


# Get the squared Self-Consistent Effective Ring Electron Density (SCERED^2) 
# for a given set of Kekule structures
#   RG: 2D array storing the rings of the given molecular graph
#   BD: List of numpy arrays storing the bonds for the given Kekule structures
# Return:
#   ered2: List of the SCERED^2 values of rings
#       W: List of the predicted weights for the Kekule structures
def getSCERED2( RG, BD ):
    NR = len( RG ) # Number of rings
    NK = len( BD ) # Number of Kekule structures

    MAX_ITER = 100 # Maximum number of iterations
    TOL = 1E-6 / NK # Tolerence in weights

    ringSize = [] # ring sizes

    SCERED = []
    FNRE = [] # Formal Numbers of Ring Electrons
    for bnd in BD: # Run over Kekule structures
        N = np.max( bnd ) # Number of atoms
        # Get linkage matrix:
        lm = np.zeros( ( N, N ), dtype='i' )
        for b in bnd:
            lm[ b[0]-1 ][ b[1]-1 ] = 1
            lm[ b[1]-1 ][ b[0]-1 ] = 1
        # Formal Numbers of Ring Electrons:
        fnre = np.zeros( NR, dtype='i' )
        for j in range( NR ):
            rg = RG[j]
            rgsize = len( rg )
            ringSize.append( rgsize )
            for k in range( rgsize ):
                iat1 = rg[k] - 1
                if k == rgsize-1:
                    iat2 = rg[0] - 1
                else:
                    iat2 = rg[k+1] - 1
                if lm[iat1][iat2] == 1:
                    fnre[j] += 2
        FNRE.append( fnre )

    FNRE = np.array( FNRE, dtype='f' )

    ringSign = np.zeros( NR )
    for r in range( NR ):
        if ringSize[r] == 6:
            ringSign[r] = 1
        elif ringSize[r] == 5:
            ringSign[r] = -1
        else:
            raise ValueError( 'Invalid ring size (%i) for ring %i'
                    % ( len( RG[r] ), r ) )
    #print( ringSign )

    # Self-consistently prediction of Effective Ring Electron density:
    # Initial guess of weights of Kekule structures based on FNREs:
    W = np.ones( NK ) / NK
    W0 = W
    ered2 = np.zeros( NR ) # Initialize ERED^2
    lambda_k = np.zeros( NK ) # Initialize lambda_k
    for loop in range( MAX_ITER ):
        # Get the ERED of each ring:
        for r in range( NR ):
            ered2[r] = sum( (FNRE[:,r]/ringSize[r])**2 * W )
        W_r = ered2 / ered2.sum()
        # lambda_k for each Kekule structure:
        for k in range( NK ):
            lambda_k[k] = sum( (FNRE[k,:]/ringSize[r])**2 * ringSign * W_r )
        #print( lambda_k )
        #input()
        # Evaluate the weight of each Kekule structure:
        W = lambda_k / lambda_k.sum()
        #print( W )
        #input()

        # Check convergence:
        if max( abs( W - W0 ) ) < TOL:
            print( 'Self-consistent weights converged after %i iterations' %
                    ( loop+1 ) )
            break

        W0 = W

    #print( W )

    return ( ered2, W )
# enddef getSCERED2()


# Calculate SCERED^2 and predicted weights of Kekule structures from xyz
# --  This is an interface between runJob_SCERED() and getSCERED2():
# kekFileName: *.kek file name that stores the Kekule structures
#         xyz: Cartesian coordinates of carbon atoms
#          at: Indices of carbon atoms in the original molecules (which may 
#              contain hydrogen atoms)
def calcSCERED2( kekFileName, xyz, at ):
    # Adjacent matrix:
    A = adjMat( xyz )
    # List of neighbors:
    nblist = neighbors( A )
    # Ring list:
    RG = rings( nblist, at )         
    NR = len( RG )

    _, BD, _, _ = readKekule( kekFileName )
    NK = len( BD ) # Number of Kekule structures

    ( SCERED2, W_pred ) = getSCERED2( RG, BD )
    #print( SCERED2, W_pred )
    W_pred *= 100

    # Print out:
    print( 'SCERED for the rings:' )
    print( '-'*60 )
    print( '%5s  %12s    %s' % ('Ring', 'SCERED^2', 'Atoms') )
    for r in range( NR ):
        print( '%5i  ' % (r+1), end='' )
        print( '%12.8f ' % SCERED2[r], end='' )
        for a in RG[r]:
            print( ' %3i' % a, end='' )
        print()
    print( '-'*60 )

    ndig = max( int( np.log10(NK) ) + 1, 3 )
    fmt_str_No = '%%%is' % ndig
    fmt_num_No = '%%%ii' % ndig
    print( 'Predicted weights for the Kekule structures:' )
    print( '-'*80 )
    print( fmt_str_No % 'No.', end='' )
    print( '  %12s  %s' % ('Weight(%)', 'Kekule structure') )
    for k in range( NK ):
        print( fmt_num_No % (k+1), end='' )
        if W_pred[k] < 1E-4:
            print( '  %12.4E' % W_pred[k], end='' )
        else:
            print( '  %12.8f' % W_pred[k], end='' )
        print( ' ', end='' )
        for bnd in BD[k]:
            print( ' %i-%i' % ( bnd[0], bnd[1] ), end='' )
        print()
    print( '-'*80 )

    return
# enddef calcSCERED2()
