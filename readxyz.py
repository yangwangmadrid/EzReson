# Read *.xyz file for Cartesian coordinates of atoms
#
#  Created on Apr 23, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#

import re
import numpy as np

#==============================================================================
# Read an xyz file
#
# Return:
#   elem: list of element symbols
#   xyz:  Cartesian coordinates as an N*3 matrix
def readxyz( inputFileName ):
    iline = 0
    elem =[]
    xyz =[]
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            # Number of atoms:
            if iline == 1:
                NAt = int( (line.split())[0] )
                #print( 'NAt = %i' % NAt )
                continue

            # Title:
            if iline == 2:
                #print( 'Title = %s' % line )
                continue

            # Coordinates:
            fields = line.split()
            if len(fields) >= 4:
                elem.append( fields[0] )
                xyz.append( np.array([ float( fields[1] ), float( fields[2] ), \
                                       float (fields[3] ) ]) )
    # Check consistency of number of atoms:
    NAt1 = len( elem )
    if NAt1 != NAt:
        raise ValueError( 'Expected %i atoms, but %i are found' % (NAt, NAt1) )
    xyz = np.array( xyz )

    return ( elem, xyz )
# enddef readFchk()
