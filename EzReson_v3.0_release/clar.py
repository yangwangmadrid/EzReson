# Enumeration of all Clar resonators of a given molecular graph
#
#  Created on Dec 23, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    Aug 10, 2021
#    - Added function clarMax() to  enumerate Clar structures
#
#    Jul 28, 2021
#    - Added function getFNREClar() to get number of double bonds and number
#    of adjacent Clar sextets to the rings 
#
#    Jul 20, 2021
#    - Improved greatly (by more than four times) the efficiency of enumeration
#      of Clar resonators, and added function removeExactBonds() in topo.py
#
#    Jun 30, 2021
#    - Allowed enumeration of Clar resonators for addition patterns of PAHs
#      including the cases where there are bonds belonging to no rings
#
#    Dec 28, 2020
#    - Fixed a bug in dertermining sextet in the remaining Kekule bonds
#
#    Dec 27, 2020
#    - Functions clar() & enumClar() started working
#

global g_MAX_NSX
g_MAX_NSX = 0


from topo import *
from util import *
import copy as cp
import os
import re


# Recursive enumeration of all Kekule structures of a molecular fragment, with
# returned results as a list stored in memory:
def kekuleList( keklist, K, B, NAt, iat ):
    # Number of bonds in the list:
    NB = len( B )

    # Finished:
    if NB == 0:       
        Knew = K
        Bnew = B

        # If this is an invalid Kekulean path:
        NK = len( K )
        if NK != NAt // 2:
            return ( keklist, Knew, Bnew )

        # NOTE: Output atomic indices starting from 1 (not 0):
        keklist.append( K )

        return keklist, Knew, Bnew

    # In case that atom iat is not present in the bond list B:
    for i in range( NB ):
        if B[i][0] == iat or B[i][1] == iat:
            break
        keklist, Knew, Bnew = kekuleList( keklist, K, B, NAt, iat+1 )
        return keklist, Knew, Bnew

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
        Bnew = removeExactBonds( B, bnd )
        # Recursive calling:
        keklist, Knew, Bnew = kekuleList( keklist, Knew, Bnew, NAt, iat+1 );

    return keklist, Knew, Bnew
# enddef kekuleList()


# Recursive enumeration of all Kekule structures of a molecular fragment, with
# the condition that no sextets are present in any of the Kekule structures
# The returned results is a list stored in memory:
#       Rg6: A list to store indices of sextet rings;
#             1: the corresponding ring is a Clar's sextet
#             0: the corresponding ring is not a Clar's sextet
#    rglist: List of atomic indices for all rings
#   keklist: List of Kekule strutures (storing atomic indices of all bonds)
#         B: List of atomic indices of all bonds
#         A: List of indices of all atoms
#       iRg: Index of the starting atom to search paths of possible Kekule
#            structures. For initial call of this function, set iat as 0
def kekuleList_nonSextets( Rg6, rglist, keklist, K, B, A, atIx, iat ):
    #print( 'K = ', end='' )
    #for k in K:
    #    print( ' %i-%i' % ( k[0]+1, k[1]+1 ), end='' )
    #print()
    #print( 'B = ', end='' )
    #for b in B:
    #    print( ' %i-%i' % ( b[0]+1, b[1]+1 ), end='' )
    #print()
    #print( 'A = ', end='' )
    #for a in A:
    #    print( ' %i' % ( a+1 ), end='' )
    #print()
    #input()

    # Number of bonds in the list:
    NB = len( B )

    # Finished:
    if NB == 0:   
        # If this is an invalid Kekulean path:
        NK = len( K )
        NAt = len( atIx ) # Total number of atoms
        if NK != NAt // 2:
            return keklist, K, B, A

        # If this is not a sextet-free Kekulean path:
        NRg6 = len(Rg6)
        existSextet = False
        for i in range( NRg6 ):
            if Rg6[i] == 0 and len( rglist[i] ) == 6:
                nBond = 0
                for bnd in K:
                    #print( bnd[0], bnd[1], rglist[i] )
                    if bnd[0] in rglist[i] and bnd[1] in rglist[i]:
                        nBond += 1
                        if nBond == 3:
                           existSextet = True
                           #print( ' Found a sextet: ', rglist[i] )
                           return keklist, K, B, A

        # Now, accept the valid Kekulean path:
        # NOTE: Output atomic indices starting from 1 (not 0):
        keklist.append( K )
        return keklist, K, B, A


    # If iat has already been determined, then continue for the next atom:
    if atIx[ iat ] not in A:
        keklist, Knew, Bnew, Anew = kekuleList_nonSextets( Rg6, rglist, \
                keklist, K, B, A, atIx, iat+1 );
        return keklist, Knew, Bnew, Anew


    # Find all bonds involving atom iat:
    nbond_iat = 0
    B_iat = []
    for i in range( NB ):
        if B[i][0] == atIx[ iat ] or B[i][1] == atIx[ iat ]:
            B_iat.append( B[i] )
            nbond_iat += 1
            if nbond_iat >= 3: # At most tricoordinate:
                break

    # Run over all bonds involving atom iat:
    nbond_iat = 0
    Knew = cp.copy( K )
    Bnew = cp.copy( B )
    Anew = cp.copy( A )
    #print( A )
    #print( B )
    #print( K )
    #print( atIx[iat], B_iat )
    #input()
    for bnd in B_iat:
        # Update list of Kekule structures:
        Knew = cp.copy( K )
        Knew.append( bnd )
        # Remove atoms bnd[0] and bnd[1]:
        Anew = cp.copy( A )
        Anew.remove( bnd[0] )
        Anew.remove( bnd[1] )
        # Remove all bonds involving atom iat from bond list B:
        Bnew = cp.copy( B )
        for b in B_iat:
            Bnew = removeExactBonds( Bnew, b )
        # Remove all bonds involving atom iat2 (which is bonded to iat) from 
        # bond list B:
        if atIx[ iat ] == bnd[0]:
            iat2 = bnd[1]
        else:
            iat2 = bnd[0]
        # Find all bonds involving atom iat2:
        nbond_iat = 0
        B_iat2 = []
        for i in range( len( Bnew ) ):
            if Bnew[i][0] == iat2 or Bnew[i][1] == iat2:
                B_iat2.append( Bnew[i] )
                nbond_iat += 1
                if nbond_iat >= 3: # At most tricoordinate:
                    break
        for b in B_iat2:
            Bnew = removeExactBonds( Bnew, b )

        # Recursive calling:
        keklist, Knew, Bnew, Anew = kekuleList_nonSextets( Rg6, rglist, \
                keklist, Knew, Bnew, Anew, atIx, iat+1 );

    return keklist, Knew, Bnew, Anew
# enddef kekuleList_nonSextets()



# Recursive enumeration of all Clar resonators:
#
#    writer: For writing each Clar's resonantor to external file once detected
#       Rg6: A list to store indices of sextet rings;
#             1: the corresponding ring is a Clar's sextet
#             0: the corresponding ring is not a Clar's sextet
#            -1: the corresponding ring is undetermined yet
#    rglist: List of atomic indices for all rings
#  nbrglist: List of ring indices for all neighboring rings of each ring
#         B: Bond list
#   atNonRg: List of the atoms that belong to no rings:
#       iRg: Index of the starting ring to search paths of possible Clar
#            resonators. For initial call of clar(), set iRg as 0
def clar( writer, Rg6, rglist, nbrglist, B, atNonRg, iRg ):
    #print( 'New fast algorithm for clar() ...' )

    # Finished:
    if iRg == len( rglist ):
        #------------------------------------------------------------#
        #*** Find all Kekule structures from the remaining atoms: ***
        # Get all remaining atoms:
        AtLeft = atNonRg
        NRg6 = len(Rg6)
        for i in range( NRg6 ):
            if Rg6[i] == 0:
                for a in rglist[i]:
                    if a not in AtLeft:
                        AtLeft.append( a )
        for i in range( NRg6 ):
            if Rg6[i] == 1:
                for a in rglist[i]:
                    if a in AtLeft:
                        AtLeft.remove( a )
        NAtLeft = len( AtLeft )

        # Get nblist for all remaining atoms:
        BLeft = []
        for bnd in B:
            if bnd[0] in AtLeft and bnd[1] in AtLeft:
                BLeft.append( bnd )

        # Determine the remaining Kekule bonds:
        keklist, Knew, Bnew, Anew = kekuleList_nonSextets( \
                Rg6, rglist, [], [], BLeft, AtLeft, AtLeft, 0 )

        # Now, this is a valid Clar's structure and is to be output:
        for K in keklist:
            # Write rings of sextets first:
            for i in [ j for j, x in enumerate(Rg6) if x > 0 ]:
                writer.write( ' %i' % (i+1) )
            # A terminal zero to separate sextets & Kekule bonds
            writer.write( ' 0' )
            # Write Kekule bonds:
            for bnd in K:
                writer.write( ' %i %i' % (bnd[0]+1, bnd[1]+1) )
            writer.write( '\n' )
        #------------------------------------------------------------#
        return


    # If iRg has already been considered previously, then check & skip:
    if Rg6[ iRg ] == 1: # This is a Clar's sextet
        for nb in nbrglist[ iRg ]:
            if Rg6[ nb ] == 1: # Abort this path
                return
            if Rg6[ nb ] == -1: # Unexplored neighbor:
                Rg6[ nb ] = 0
        clar( writer, Rg6, rglist, nbrglist, B, atNonRg, iRg+1 )
        return
    elif Rg6[ iRg ] == 0: # This is not a Clar's sextet
        clar( writer, Rg6, rglist, nbrglist, B, atNonRg, iRg+1 )
        return

    # Only consider six-membered, because we are talking about sextets:
    if len( rglist[ iRg ] ) != 6:
        Rg6[ iRg ] = 0
        clar( writer, Rg6, rglist, nbrglist, B, atNonRg, iRg+1 )
        return


    # There are only two possibilities for ring iRg:
    # -- 1. It forms a sextet:
    Rg61 = Rg6.copy()
    Rg61[ iRg ] = 1
    # Set all neighboring rings:
    flag = True
    for nb in nbrglist[ iRg ]:
        if Rg61[ nb ] == 1: # Abort this path
            flag = False
            break
        if Rg61[ nb ] == -1: # Unexplored neighbor:
            Rg61[ nb ] = 0
    if flag:
        clar( writer, Rg61, rglist, nbrglist, B, atNonRg, iRg+1 )

    # -- 2. It is not a sextet:
    Rg6[ iRg ] = 0
    clar( writer, Rg6, rglist, nbrglist, B, atNonRg, iRg+1 )

    return
# enddef clar()


# Recursive enumeration of all Clar structures (i.e., Clar resonators with
# maxinum number of Clar sextets)
#
#    writer: For writing each Clar's resonantor to external file once detected
#       Rg6: A list to store indices of sextet rings;
#             1: the corresponding ring is a Clar's sextet
#             0: the corresponding ring is not a Clar's sextet
#            -1: the corresponding ring is undetermined yet
#    rglist: List of atomic indices for all rings
#  nbrglist: List of ring indices for all neighboring rings of each ring
#         B: Bond list
#   atNonRg: List of the atoms that belong to no rings:
#       iRg: Index of the starting ring to search paths of possible Clar
#            resonators. For initial call of clarMax(), set iRg as 0
def clarMax( writer, Rg6, rglist, nbrglist, B, atNonRg, iRg ):

    # Finished:
    if iRg == len( rglist ):
        # First, check if the number of sextets is no smaller than the current
        # maximum value (g_MAX_NSX)
        global g_MAX_NSX
        NSX = len( np.argwhere( np.array(Rg6) > 0 ) ) # Number of sextets
        if NSX < g_MAX_NSX:
            return

        #------------------------------------------------------------#
        #*** Find all Kekule structures from the remaining atoms: ***
        # Get all remaining atoms:
        AtLeft = atNonRg
        NRg6 = len(Rg6)
        for i in range( NRg6 ):
            if Rg6[i] == 0:
                for a in rglist[i]:
                    if a not in AtLeft:
                        AtLeft.append( a )
        for i in range( NRg6 ):
            if Rg6[i] == 1:
                for a in rglist[i]:
                    if a in AtLeft:
                        AtLeft.remove( a )
        NAtLeft = len( AtLeft )

        # Get nblist for all remaining atoms:
        BLeft = []
        for bnd in B:
            if bnd[0] in AtLeft and bnd[1] in AtLeft:
                BLeft.append( bnd )

        # Determine the remaining Kekule bonds:
        keklist, Knew, Bnew, Anew = kekuleList_nonSextets( \
                Rg6, rglist, [], [], BLeft, AtLeft, AtLeft, 0 )

        if len( keklist ) == 0:
            return

        # Now, this is a valid Clar's structure:
        ### 1. Check if the number of sextets is the new maximum:
        if NSX > g_MAX_NSX:
            g_MAX_NSX = NSX
        ### 2. Output:
        for K in keklist:
            # Write rings of sextets first:
            for i in [ j for j, x in enumerate(Rg6) if x > 0 ]:
                writer.write( ' %i' % (i+1) )
            # A terminal zero to separate sextets & Kekule bonds
            writer.write( ' 0' )
            # Write Kekule bonds:
            for bnd in K:
                writer.write( ' %i %i' % (bnd[0]+1, bnd[1]+1) )
            writer.write( '\n' )
        #------------------------------------------------------------#
        return


    # If iRg has already been considered previously, then check & skip:
    if Rg6[ iRg ] == 1: # This is a Clar's sextet
        for nb in nbrglist[ iRg ]:
            if Rg6[ nb ] == 1: # Abort this path
                return
            if Rg6[ nb ] == -1: # Unexplored neighbor:
                Rg6[ nb ] = 0
        clarMax( writer, Rg6, rglist, nbrglist, B, atNonRg, iRg+1 )
        return
    elif Rg6[ iRg ] == 0: # This is not a Clar's sextet
        clarMax( writer, Rg6, rglist, nbrglist, B, atNonRg, iRg+1 )
        return

    # Only consider six-membered, because we are talking about sextets:
    if len( rglist[ iRg ] ) != 6:
        Rg6[ iRg ] = 0
        clarMax( writer, Rg6, rglist, nbrglist, B, atNonRg, iRg+1 )
        return


    # There are only two possibilities for ring iRg:
    # -- 1. It forms a sextet:
    Rg61 = Rg6.copy()
    Rg61[ iRg ] = 1
    # Set all neighboring rings:
    flag = True
    for nb in nbrglist[ iRg ]:
        if Rg61[ nb ] == 1: # Abort this path
            flag = False
            break
        if Rg61[ nb ] == -1: # Unexplored neighbor:
            Rg61[ nb ] = 0
    if flag:
        clarMax( writer, Rg61, rglist, nbrglist, B, atNonRg, iRg+1 )

    # -- 2. It is not a sextet:
    Rg6[ iRg ] = 0
    clarMax( writer, Rg6, rglist, nbrglist, B, atNonRg, iRg+1 )

    return
# enddef clarMax()


# Enumeration of all Clar resonators from Cartesian coordinates:
#
#  ifClarMax: If enumerate Clar structures, otherwise Clar resonators are
#             enumerated
#     atExcl: List of indices (starting from 1) of the atoms that are 
#             excluded in the pi-system
#
def enumClar( writer, elem, xyz, ifClarMax, \
        atExcl=[], clarFileName_atExcl='' ):
    xyz, at = xyzNonH( elem, xyz ) # For non-hydrogen atoms

    # Exclude the atoms in atExcl[]:
    xyz0 = xyz.copy()
    at0 = at.copy()
    xyz = []
    at = []
    for i in range( len(at0) ):
        if i+1 not in atExcl:
            xyz.append( xyz0[i,:] )
            at.append( at0[i] )
    xyz = np.array( xyz )
    NAt = xyz.shape[0]
    # Write xyz to xyzFileName_atExcl:
    if len( atExcl ) > 0:
        xyzFileName_atExcl = re.sub( '\.clar$', '.xyz', clarFileName_atExcl )
        with open( xyzFileName_atExcl, 'w' ) as fw:
            fw.write( '%i\n' % NAt )
            for ix in atExcl:
                fw.write( ' %i' % ix )
                fw.write( '\n' )
            for i in range( NAt ):
                fw.write( 'C %12.6f %12.6f %12.6f\n' % 
                        ( xyz[i,0], xyz[i,1] ,xyz[i,2] ) )
        print( '---File %s written' % xyzFileName_atExcl )


    # Let at[] be the natural appearance indices of atoms:
    # This an ad hoc code, future generalization is needed for coordinates with
    # randomly ordered carbon atoms.
    at = range( NAt )

    # Adjacent matrix:
    A = adjMat( xyz )

    # List of neighbors:
    nblist = neighbors( A )
    #i = 0
    #for nb in nblist:
    #    i += 1
    #    print( '%-5i: ' % i, end='' )
    #    for ix in nb:
    #        print( ' %i' % (ix+1), end='' )
    #    print()


    # Bond list:
    B = bonds( A, at )
    #print( B )

    # Ring list:
    rglist = rings( nblist, at )
    # Write rings:
    NRg = len( rglist )
    writer.write( '%i\n' % NRg )
    for i in range( NRg ):
        for a in rglist[i]:
            writer.write( ' %i' % (a+1) )
        writer.write( '\n' )

    # Find the atoms that belong to no rings:
    atNonRg = at # Initialize
    for rg in rglist:
        for ix in rg:
            if ix in atNonRg:
                # Remove ix from atNonRg[]
                atNonRg = [ i for i in atNonRg if i != ix ]
    #print( atNonRg )
    #input()

    NRg = len( rglist ) # Number of rings
    Rg6 = [ -1 ] * NRg # Initialize sextet info list
    RgLeft = list( range( NRg ) )
    nbrglist = ringNeighbors( rglist )
    #clar( writer, Rg6, RgLeft, rglist, nbrglist, B, atNonRg, 0 )
    #print( Rg6, rglist, nbrglist, B, atNonRg )
    #input()
    if ifClarMax: # For Clar structures
        # First write to a temperary file:
        tmpFileName = ( 'tmp-%i.clar' %  os.getpid() )
        with open( tmpFileName, 'w' ) as writer_tmp:
            clarMax( writer_tmp, Rg6, rglist, nbrglist, B, atNonRg, 0 )
        print( 'Maximum possible number of Clar sextets is %s' % g_MAX_NSX )
        # Remove the Clar resonators that do not have max. num. of sextets:
        with open( tmpFileName ) as f:
            for line in f:
                NSX = 0
                for s in line.split():
                    # Read sextets:
                    if s == '0':
                        break
                    NSX += 1
                    if NSX == g_MAX_NSX:
                        writer.write( line )
        os.remove( tmpFileName )
    else: # For Clar resonantors
        clar( writer, Rg6, rglist, nbrglist, B, atNonRg, 0 )

    return
# enddef enumClar()


#  Generate *.clarmax file that includes all Clar structures (i.e., the Clar 
#  resonantors with max. num. of sextets) from the already enumerated Clar
#  resonators
#    clarmaxFileName:  the *.clarmax file where the Clar structures will be
#                      extracted, and the Clar resonantors have been already 
#                      pre-generated in  the *.clar file
#    atExcl:           list of indices (starting from 1) of the atoms that are 
#                      excluded in the pi-system
#
def clarMax_from_clar( clarmaxFileName, atExcl=[] ):
    clarFileName =  re.sub( '\.clarmax$', '.clar', clarmaxFileName )
    clarmaxtmpFileName = clarmaxFileName + '.tmp'
    fw = open( clarmaxFileName, 'w' ) # First part containing ring info.
    #print( clarmaxFileName )
    fw_tmp = open( clarmaxtmpFileName, 'w' )
    maxNSX = 0 # Maximum number of Clar's sextets
    NR = 0 # Number of rings
    iline = 0
    with open( clarFileName ) as f:
        for line in f:
            iline += 1
            if iline == 1:
                NR = int( line )
            if iline <= NR+1:
                fw.write( line )
                continue
            #print( line, end='' )
            NSX = 0
            for s in line.split():
                # Read sextets:
                if s == '0':
                    break
                NSX += 1
            if NSX > maxNSX:
                maxNSX = NSX
                fw_tmp.close()
                fw_tmp = open( clarmaxtmpFileName, 'w' )
                fw_tmp.write( line )
            elif NSX == maxNSX:
                fw_tmp.write( line )
    fw.close()
    fw_tmp.close()
    #print( maxNSX)

    # Append file fw_tmp to fw:
    with open( clarmaxFileName, 'a' ) as outp:
        with open( clarmaxtmpFileName ) as inp:
            for line in inp:
                outp.write( line )
    # Remove the temporary file:
    os.remove( clarmaxtmpFileName )

    return
# enddef clarMax_from_clar()



# Enumeration of Clar resonators from a given Kekule structure
#
#    writer: For writing each Clar's resonator to external file once detected
#    rglist: List of atomic indices for all rings
#  nbrglist: List of ring indices for all neighboring rings of each ring
#            NOTE: nbrglist is only among Fries sextets, i.e., ixRg_FS[]
#       bnd: List of all bonds that define the Kekule structure
#   ixRg_FS: Ring indices of Fries sextets
#       Rg6: A list to store indices of sextet rings;
#             1: the corresponding ring is a Clar's sextet
#             0: the corresponding ring is not a Clar's sextet
#            -1: the corresponding ring is undetermined yet
#            NOTE: Rg6 is only among Fries sextets, i.e., ixRg_FS[]
#       iRg: Index of the starting ring to search paths of possible Clar
#            resonators. For the initial call, set iRg as 0
def clar_from_kekule( writer, rglist, nbrglist, bnd, bnd_FS, ixRg_FS, Rg6, \
        iRg ):
    # Finished:
    if iRg == len( ixRg_FS ):
        # Get the remaining Kekule bonds by remove all bonds of Clar sextets
        # from bnd[]:
        for i in range( len(Rg6) ):
            if Rg6[i] == 1:
                for b in bnd_FS[i]:
                    bnd = removeExactBonds( bnd, b )
        # Find if there exist any sextet in the remaining Kekule bonds:
        for i in range( len(Rg6) ):
            if Rg6[i] == 0:
                nBond = 0
                for b in bnd:
                    if b[0] in rglist[ ixRg_FS[i] ] and \
                            b[1] in rglist[ ixRg_FS[i] ]:
                        nBond += 1
                if nBond == 3:
                    #print( ' Found a sextet: ring', ixRg_FS[i]+1 )
                    return

        # Write rings of sextets first:
        for i in [ j for j, x in enumerate(Rg6) if x > 0 ]:
            writer.write( ' %i' % ( ixRg_FS[i]+1 ) )
        # A terminal zero to separate sextets & Kekule bonds
        writer.write( ' 0' )
        # Write Kekule bonds:
        for b in bnd:
            writer.write( ' %i %i' % ( b[0]+1, b[1]+1 ) )
        writer.write( '\n' )
        #input()
        return


    # If iRg has already been considered previously, then check & skip:
    if Rg6[ iRg ] == 1: # This is a Clar's sextet
        for nb in nbrglist[ iRg ]:
            if Rg6[ nb ] == 1: # Abort this path
                return
            if Rg6[ nb ] == -1: # Unexplored neighbor:
                Rg6[ nb ] = 0
        clar_from_kekule( writer, rglist, nbrglist, bnd, bnd_FS, ixRg_FS, \
                Rg6, iRg+1 )
        return
    elif Rg6[ iRg ] == 0: # This is not a Clar's sextet
        clar_from_kekule( writer, rglist, nbrglist, bnd, bnd_FS, ixRg_FS, \
                Rg6, iRg+1 )
        return


    # There are only two possibilities for the sextet ring iRg:
    # -- 1. It forms a sextet:
    Rg61 = Rg6.copy()
    Rg61[ iRg ] = 1
    # Set all neighboring rings:
    flag = True
    for nb in nbrglist[ iRg ]:
        if Rg61[ nb ] == 1: # Abort this path
            flag = False
            break
        if Rg61[ nb ] == -1: # Unexplored neighbor:
            Rg61[ nb ] = 0
    if flag:
        clar_from_kekule( writer, rglist, nbrglist, bnd, bnd_FS, ixRg_FS, \
                Rg61, iRg+1 )

    # -- 2. It is not a sextet:
    Rg6[ iRg ] = 0
    clar_from_kekule( writer, rglist, nbrglist, bnd, bnd_FS, ixRg_FS, Rg6, \
            iRg+1 )

    return
# enddef clar_from_kekule()


# Enumeration of Clar resonators from the Kekule structures in a *.kek file:
#
#    writer: For writing each Clar's resonator to external file once detected
#    rglist: List of atomic indices for all rings
#  nbrglist: List of ring indices for all neighboring rings of each ring
#         N: Total number of atoms
def clar_by_reading_kekule( kekFileName, writer, rglist, nbrglist, N ):
    print( 'Enumeration of Clar resonators from the Kekule structures in '
            'file %s' % kekFileName )

    NR = len( rglist ) # Number of rings

    with open( kekFileName ) as f:
        for line in f:
            # For each kekule structure:
            bnd = []
            i = 0
            for s in line.split():
                if i == 0:
                    i = 1
                    b1 = int(s) - 1 # Index staring from 0, NOT 1
                else:
                    i = 0
                    b2 = int(s) - 1 # Index staring from 0, NOT 1
                    bnd.append( [ b1, b2 ] )
            #print( bnd )
            # Get linkage matrix:
            lm = np.zeros( ( N, N ), dtype='i' )
            for b in bnd:
                lm[ b[0] ][ b[1] ] = 1
                lm[ b[1] ][ b[0] ] = 1
            # Find the rings that are Fries sextets:
            ixRg_FS = [] # Ring indices of Fries sextets
            bnd_FS = [] # All three bonds in each of the Fries sextets
            for j in range( NR ):
                rgsize = len( rglist[j] )
                b = []
                for k in range( rgsize ):
                    iat1 = rglist[j][k]
                    if k == rgsize-1:
                        iat2 = rglist[j][0]
                    else:
                        iat2 = rglist[j][k+1]
                    if lm[ iat1 ][ iat2 ] == 1:
                        if iat1 < iat2:
                            b.append( [ iat1, iat2 ] )
                        else:
                            b.append( [ iat2, iat1 ] )
                if len( b ) == 3:
                    ixRg_FS.append( j )
                    bnd_FS.append( b )
            #print( bnd_FS )
            # Get neighboring ring list only among ixRg_FS[]:
            ixRg_FS = np.array( ixRg_FS )
            nbrglist_FS = []
            for i in ixRg_FS:
                nb_FS = []
                for nb in nbrglist[i]:
                    ix = np.argwhere( ixRg_FS == nb )
                    if len( ix ) > 0:
                        nb_FS.append( int( ix[0] ) )
                nbrglist_FS.append( nb_FS )
            #print( ixRg_FS )
            #print( nbrglist_FS )
            Rg6 = [ -1 ] * len( ixRg_FS )  # Initialize sextet info list
            clar_from_kekule( writer, rglist, nbrglist_FS, bnd, bnd_FS, \
                    ixRg_FS, Rg6, 0 )
            #input()

    return
# enddef clar_by_reading_kekule()


# Enumeration of all Clar resonators from Kekule structures:
# NOTE: This method is much less (by three times) efficient than enumClar().
# So, we should abandon it.
# For this reason, the code is NOT finished in the sense that the  generated
# *.clar file contains duplicated Clar resonantors and need to be further
# uniq-ed.
#
#    kekFileName:  the *.kek file where the Kekule structures have been already 
#                  pre-generated
#         atExcl:  list of indices (starting from 1) of the atoms that are 
#                  excluded in the pi-system
def enumClar_from_Kekule( kekFileName, writer, elem, xyz, atExcl=[], \
        clarFileName_atExcl=''  ):
    xyz, at = xyzNonH( elem, xyz ) # For non-hydrogen atoms

    # Exclude the atoms in atExcl[]:
    xyz0 = xyz.copy()
    at0 = at.copy()
    xyz = []
    at = []
    for i in range( len(at0) ):
        if i+1 not in atExcl:
            xyz.append( xyz0[i,:] )
            at.append( at0[i] )
    xyz = np.array( xyz )
    NAt = xyz.shape[0]
    # Write xyz to xyzFileName_atExcl:
    if len( atExcl ) > 0:
        xyzFileName_atExcl = re.sub( '\.clar$', '.xyz', clarFileName_atExcl )
        with open( xyzFileName_atExcl, 'w' ) as fw:
            fw.write( '%i\n' % NAt )
            for ix in atExcl:
                fw.write( ' %i' % ix )
                fw.write( '\n' )
            for i in range( NAt ):
                fw.write( 'C %12.6f %12.6f %12.6f\n' % 
                        ( xyz[i,0], xyz[i,1] ,xyz[i,2] ) )
        print( '---File %s written' % xyzFileName_atExcl )


    # Let at[] be the natural appearance indices of atoms:
    # This an ad hoc code, future generalization is needed for coordinates with
    # randomly ordered carbon atoms.
    at = range( NAt )

    # Adjacent matrix:
    A = adjMat( xyz )

    # List of neighbors:
    nblist = neighbors( A )

    # Bond list:
    B = bonds( A, at )
    #print( B )

    # Ring list:
    rglist = rings( nblist, at )
    # Write rings:
    NRg = len( rglist )
    writer.write( '%i\n' % NRg )
    for i in range( NRg ):
        for a in rglist[i]:
            writer.write( ' %i' % (a+1) )
        writer.write( '\n' )

    # Find the atoms that belong to no rings:
    atNonRg = at # Initialize
    for rg in rglist:
        for ix in rg:
            if ix in atNonRg:
                # Remove ix from atNonRg[]
                atNonRg = [ i for i in atNonRg if i != ix ]
    #print( atNonRg )
    #input()

    NRg = len( rglist ) # Number of rings
    Rg6 = [ -1 ] * NRg # Initialize sextet info list
    RgLeft = list( range( NRg ) )
    nbrglist = ringNeighbors( rglist )

    clar_by_reading_kekule( kekFileName, writer, rglist, nbrglist, NAt )

    return
# enddef enumClar_from_Kekule()


# Read Clar resonators from *.clar file and convert them to LP[] and BD[]
#   clarFileName: Name of *.clar file where Clar resonators are stored
#   div:         If div = (n,m), then divide clarFileName in n parts and only
#                consider the m-th part of the whole set of Clar resonators
def readClar( clarFileName, div=(1,1) ):
    # Get the rings:
    NR = 0 # Number of rings
    RG = []
    iline = 0
    with open( clarFileName ) as f:
        for line in f:
            iline += 1
            if iline == 1:
                NR = int( line )
                continue
            if iline > NR+1:
                break

            # Parse the atoms in each ring:
            RG.append( [ int(x) for x in line.split() ] )
    #print( RG )
    #input()

    # Division of the whole file:
    n, m = div  # Divide into n parts and take the m-th partition
    # Number of lines of Clar resonators in the file:
    nline = sum( 1 for line in open( clarFileName ) ) - (NR+1)
    wid = nline // n
    istart = (m-1) * wid + 1  # NOTE: istart == 1 (not 0) when m == 1
    iend = istart + wid - 1
    # Take care of the residual lines:
    if m == n and iend < nline:
        iend = nline

    LP = [] # Lone pairs (currently not supported)
    BD = [] # Bonds not foming a Clar's sextet
    SX = [] # Clar's sexets
    iline = 0
    iline0 = 0
    with open( clarFileName ) as f:
        for line in f:
            iline += 1
            iline0 += 1
            # Skip the lines for ring info:
            if iline0 <= NR + 1:
                iline = 0
                continue

            if iline < istart or iline > iend:
                continue
            flag = 0 # 0 for reading sextets; 1 for reading bonds
            bnd = []
            sxt = []
            i = 0
            for s in line.split():
                # Read sextets:
                if flag == 0:
                    # End of reading sextets:
                    if s == '0':
                        flag = 1
                        continue
                    sxt.append( int(s) )
                # Read bonds:
                else:
                    if i == 0:
                        i = 1
                        b1 = int(s)
                    else:
                        i = 0
                        b2 = int(s)
                        bnd.append( [ b1, b2 ] )
            SX.append( np.array( sxt, dtype='i' ) )
            BD.append( np.array( bnd, dtype='i' ) )
            LP.append( np.array( [], dtype='i' ) )
    return LP, BD, SX, RG
# enddef readClar()


# String of a given Clar resonator:
def clar_str( lp, bd, sx, RG, ifExplicitRings=False ):
    s = '' # Initialization

    if len( sx ) > 0:
        for k in sx:
            if not ifExplicitRings: # Only indicate ring index
                s += '%i ' % k
            else: # Give explicit atom indices of ring
                rgSz = len( RG[k-1] ) # NOTE: ring ix. in sx start from 1 not 0
                for i in range( rgSz ):
                    s += '%i' % RG[k-1][i]
                    if i < rgSz - 1:
                        s += '='
                    else:
                        s += ' '

    if len( bd ) > 0:
        i = 1
        bd1 = bd.copy()
        bd1 = bd1.reshape( bd.shape[0]*bd.shape[1] )
        nBD = len( bd1 )
        for k in bd1:
            s += '%i' % k
            if i % 2 == 1:
                s += '-'
            elif i < nBD:
                s += ' '
            i += 1

    if len( lp ) > 0:
        s = ': '.join( map( str, lp ) ) + ': '
 
 
    return s
# enddef clar_str()


# Get Formal Number of Ring Electrons considering Clar Sextets on the rings 
# for a given set of Clar resonantors
#         N: Total number of (carbon) atoms
#        RG: 2D array storing the rings of the given molecular graph
#  nbrglist: List of ring indices for all neighboring rings of each ring
#        BD: List of numpy arrays storing the bonds for the given Clar reson.
#        SX: List of numpy arrays storing the sextets for the given Clar reson.
# Return:
#       NDB: Number of Double Bonds on each ring (NDB)
#       NCS: Number of Clar Sextets adjacent to each ring (NCS)
def getFNREClar( N, RG, nbrglist, BD, SX ):
    NR = len( RG )
    NClar = len( SX ) # Number of Clar resonators

    NDB = [] # Number of Double Bonds on each ring (NDB)
    NCS = [] # Number of Clar Sextets adjacent to each ring (NCS)
    for iClar in range( NClar ): # Run over Clar resonators
        sxt = SX[iClar]
        bnd = BD[iClar]
        # Get linkage matrix:
        lm = np.zeros( ( N, N ), dtype='i' )
        for b in bnd:
            lm[ b[0]-1 ][ b[1]-1 ] = 1
            lm[ b[1]-1 ][ b[0]-1 ] = 1
        # Number of Double Bonds (NDB):
        ndb = np.zeros( NR, dtype='i' )
        # Number of Clar Sextets adjacent to each ring:
        ncs = np.zeros( NR, dtype='i' )
        for j in range( NR ):
            if (j+1) in sxt: # Skip Clar sextets
                ndb[j] = 6 # Use 6 to denote a Clar sextet
                ncs[j] = 0
                continue
            # Get number of double bonds:
            rg = RG[j]
            rgsize = len( rg )
            for k in range( rgsize ):
                iat1 = rg[k] - 1
                if k == rgsize-1:
                    iat2 = rg[0] - 1
                else:
                    iat2 = rg[k+1] - 1
                if lm[iat1][iat2] == 1:
                    ndb[j] += 1
            # Get number of adjacent Clar sextets:
            for nbrg in nbrglist[j]:
                if (nbrg+1) in sxt:
                    ncs[j] += 1
        NDB.append( ndb ) 
        NCS.append( ncs ) 
    return ( np.array( NDB ), np.array( NCS ) )
# enddef getFNREClar()


# Get the squared Self-Consistent Effective Ring Electron Density (SCERED^2) 
# for a given set of Clar resonators
#         N: Total number of (carbon) atoms
#        RG: 2D array storing the rings of the given molecular graph
#  nbrglist: List of ring indices for all neighboring rings of each ring
#        BD: List of numpy arrays storing the bonds for the given Clar reson.
#        SX: List of numpy arrays storing the sextets for the given Clar reson.
#       NDB: Number of Double Bonds on each ring (NDB)
#       NCS: Number of Clar Sextets adjacent to each ring (NCS)
# Return:
#   ered2: List of the SCERED^2 values of rings
#       W: List of the predicted weights for the Kekule structures
def getSCERED2_Clar( RG, NDB, NCS ):

    return
# enddef getSCERED2_Clar()
