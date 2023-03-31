# Read *.xyz file for Cartesian coordinates of atoms
#
#  Created on Apr 23, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    Jan 12, 2021:
#    - Added function readgjf() to read coordinates from Gaussian's input file
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
# enddef readxyz()


#==============================================================================
# Read a gjf file (Gaussian's input file)
#
# Return:
#   elem: list of element symbols
#   xyz:  Cartesian coordinates as an N*3 matrix
def readgjf( inputFileName ):
    pat_header = re.compile( '^\s*#.*' )
    pat_empty = re.compile( '^\s*$' )

    iline = 0
    elem =[]
    xyz =[]
    flag = 0
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            # Header:
            if flag == 0 and pat_header.match( line ):
                flag = 1
                continue

            # Empty lines:
            if flag > 0 and pat_empty.match( line ):
                flag += 1
                continue

            # Charge & multiplicity:
            if flag == 3:
                #print( line )
                flag = 4
                continue

            # Coordinates:
            if flag == 4:
                fields = line.split()
                if len(fields) >= 4:
                    elem.append( fields[0] )
                    xyz.append( np.array([ float( fields[1] ), \
                            float( fields[2] ), float (fields[3] ) ]) )
                else:
                    break
    xyz = np.array( xyz )

    return ( elem, xyz )
# enddef readgjf()


#==============================================================================
# Read an out file (Gaussian's output file)
#
# Return:
#   elem: list of element symbols
#   xyz:  Cartesian coordinates as an N*3 matrix
def readGout( inputFileName ):
    pat_header = re.compile( '^\s*Number\s+Number\s+Type\s+X\s+Y\s+Z\s*$' )

    iline = 0
    elem =[]
    xyz =[]
    flag = 0
    # Get the LAST header line for the coordinates section:
    iline_start = 0
    iline = 0
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            # Header:
            if pat_header.match( line ):
                iline_start = iline
                continue

    if iline_start == 0:
        raise RuntimeError( 'No atomic coordinates found in file %s' %
                inputFileName )
    iline_start += 2

    iline = 0
    elem =[]
    xyz =[]
    flag = 0
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            # Reading the atomic coordinates:
            if iline >= iline_start:
                fields = line.split()
                if len( fields ) != 6:
                    break
                elem.append( fields[1] )
                xyz.append( [ float( fields[3] ), float( fields[4] ), \
                        float( fields[5] ) ] )

    xyz = np.array( xyz )

    return ( elem, xyz )
# enddef readGout()
