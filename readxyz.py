# Read *.xyz, *.gjf files for Cartesian coordinates of atoms
#
#  Updates:
#
#    Sep 5, 2020:
#    - Added function readgjf() and generalized readxyz() to read atomic 
#      coordinates from a Guassian's *.gjf file
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
    # Check if it is a Gaussian input file (*.gjf):
    pat_gjf = re.compile( '^\s*#\s*.*' )
    with open( inputFileName ) as f:
        for line in f:
            if pat_gjf.match( line ):
                return readgjf( inputFileName )

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
# enddef readgjf()


#==============================================================================
# Read a gjf file
#
# Return:
#   elem: list of element symbols
#   xyz:  Cartesian coordinates as an N*3 matrix
def readgjf( inputFileName ):
    iline = 0
    elem =[]
    xyz =[]
    flag = 0
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            # An empty line:
            if len( line.strip() ) == 0:
                flag += 1
                continue

            # Charge & Multiplicity:
            if flag == 2:
                if len( line.split() ) != 2:
                    print( 'ERROR: Invalid gjf format at line %i:' % iline )
                    print( line )
                    exit(1)
                flag += 1
                continue

            # Coordinates:
            if flag == 3:
                fields = line.split()
                if len(fields) >= 4:
                    elem.append( fields[0] )
                    xyz.append( np.array([ float( fields[1] ), \
                            float( fields[2] ), float (fields[3] ) ]) )

    xyz = np.array( xyz )

    return ( elem, xyz )
# enddef readgjf()
