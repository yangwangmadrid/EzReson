# Enumeration of all Kekule structures of a given molecular graph
#
#  Created on Apr 23, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#

from topo import *
from util import *
import copy as cp

# Iterative enumeration of all Kekule structures:
def kekule( writer, K, B, NAt, iat ):
    # Number of bonds in the list:
    NB = len( B )

    # Finished:
    if NB == 0:       
        Knew = K
        Bnew = B

        # If this is an invalid Kekulean path:
        NK = len( K )
        if NK != NAt // 2:
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
        Knew, Bnew = kekule( writer, K, B, NAt, iat+1 )
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
        Knew, Bnew = kekule( writer, Knew, Bnew, NAt, iat+1 );
        if iat == 0:
            print( 'Progress %.2f%% ...' % ( i/len(B_iat)*100 ) )

    return Knew, Bnew
# enddef kekule()


# Enumeration of all Kekule structures from Cartesian coordinates:
def enumKekule( writer, elem, xyz ):
    xyz, at = xyzNonH( elem, xyz ) # For non-hydrogen atoms
    NAt = xyz.shape[0]

    # Adjacent matrix:
    A = adjMat( xyz )

    # Bond list:
    B = bonds( A )

    # Assign back the original atomic labels in B[]:
    B0 = []
    for bnd in B:
        B0.append( np.array([ at[bnd[0]], at[bnd[1]] ]) )

    kekule( writer, [], B0, NAt, 0 )

    return
# enddef enumKekule()


# Read Kekule structures from *.kek file and convert them to LP[] and BD[]
#   kekFileName: Name of *.kek file where Kekule structures are stored
#   div:         If div = (n,m), then divide kekFileName in n parts and only
#                consider the m-th part of the whole set of Kekule structures
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
    iline = 0
    with open( kekFileName ) as f:
        for line in f:
            iline += 1
            if iline < istart or iline > iend:
                continue
            bnd = []
            i = 0
            for s in line.split():
                if i == 0:
                    i = 1
                    b1 = int(s)
                else:
                    i = 0
                    b2 = int(s)
                    bnd.append( [ b1, b2 ] )
            BD.append( np.array( bnd, dtype='i' ) )
            LP.append( np.array( [], dtype='i' ) )
    return LP, BD
# enddef readKekule()
