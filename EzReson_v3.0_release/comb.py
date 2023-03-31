# Enumeration of all possible combinations of lone pairs and bonds
#
#  Created on Mar 5, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    Feb 25, 2022:
#    - Extended function lewis_str() to open-shell systems
#

import numpy as np
import os.path

class LewisStructure:
    def __init__( self, NAt, NE ):
        self.stat = 0  # 0 ==> normal;  1 ==> EOF
        self.NAt = NAt
        self.NE = NE
        self.nPair = NE // 2
        self.LP = np.array( [], dtype='i' )
        self.BD = np.array( [], dtype='i' )
        self.nLP = 0
        self.nBD = 0
# endclass LewisStructure

# Total number of combinations for choosing m single objects out of n 
# distinguishable objects
#
#   N = C_n^m = n*(n-1)*...*(n-m+1) / factorial(m)
#
def numSingleComb( n, m ):
    num = 1
    for k in range( 1, m+1 ):
        num = num * (n-k+1) // k;
    return num
# enddef numSingleComb()


# Total number of combinations for choosing p pairs out of n distinguishable
# objects
#
#  N = n*(n-1)*...*(n-2*p+1)/(2^p*factorial(p) )
#
#  Derivation: 
#    N = C_{n}^2 * C_{n-2}^2 * ... * C_{n-2*p+1}^2 / factorial(p)
#      = n*(n-1)/2 * (n-2)*(n-3)/2 * ... * (n-2*p+2)*(n-2*p+1)/2 / factorial(p)
#      = n*(n-1)*...*(n-2*p+1)/(2^p*factorial(p))
#
def numPairComb( n, p ):
    num = 1
    for k in range( 1, p+1 ):
        num = num * (n-2*k+2)*(n-2*k+1) // 2 // k
    return num
# enddef numPairComb()


# Total number of combinations for choosing m single objects and p pairs, 
# at the same time, out of n distinguishable objects
#
#  N = numSingleComb( n, m ) * numPairComb( n-m, p )
#
def numAllComb( n, m, p ):
    return numSingleComb( n, m ) * numPairComb( n-m, p )
# enddef numAllComb()


#
# Total number of possible Lewis structures for an n-center-m-electon system
# NOTE: Only valid for closed-shell systems at the moment
#
def numLewis( n, m ):
    if m % 2 != 0:
        raise ValueError( 'Number of electrons must be an even number!' )
    
    max_nLP = m // 2
    num = 0
    for nLP in range( max_nLP+1 ):
        num += numAllComb( n, nLP, max_nLP - nLP )

    return num
# enddef numLewis()

# Enumerate all combinations of choosing m objects and p pairs, 
# at the same time, from n distinguishable objects
# This is a wrapper of comb_all()
#
def allComb( n, m, p ):
    a = range(n)
    b1 = np.zeros( m, dtype='i' )
    b2 = np.zeros( (p,2), dtype='i' ) 
    # NOTE: a, b[] and p-1[] are all indices starting from 0 (not from 1)
    comb_all( a, b1, b2, n, m-1, m-1, n, p-1, p-1, a )
# enddef allComb()


# Enumerate all combinations of choosing p pairs from n distinguishable objects
# This is a wrapper of comb_pair()
#
def pairComb( n, p ):
    a = range(n)
    b = np.zeros( (p,2), dtype='i' )  
    # NOTE: a, b[] and p-1[] are all indices starting from 0 (not from 1)
    comb_pair( a, b, p-1, p-1, a )
# enddef pairComb()


# Same as allComb() but output to a file:
# Enumerate all combinations of choosing m objects and p pairs, 
# at the same time, from n distinguishable objects
# This is a wrapper of comb_all()
#
def allComb_write( writer, n, m, p ):
    a = range(n)
    b1 = np.zeros( m, dtype='i' )
    b2 = np.zeros( (p,2), dtype='i' ) 
    # NOTE: a, b[] and p-1[] are all indices starting from 0 (not from 1)
    comb_all_write( writer, a, b1, b2, n, m-1, m-1, n, p-1, p-1, a )
# enddef allComb_write()

# Same as pairComb() but output to a file:
# Enumerate all combinations of choosing p pairs from n distinguishable objects
# This is a wrapper of comb_pair()
#
def pairComb_write( writer, n, p ):
    a = range(n)
    b = np.zeros( (p,2), dtype='i' )  
    # NOTE: a, b[] and p-1[] are all indices starting from 0 (not from 1)
    comb_pair_write( writer, a, b, p-1, p-1, a )
# enddef pairComb()


# Auxiliary function for comb_pair() 
# a[]: 1*n 1D array
# b_ip[,]: 1*2 2D array 
def remaining_arr( a, b_ip ):
    n = len(a)
    m = 0
    a1 = [];
    for j in range(n):
        if a[j] < b_ip[0] and a[j] != b_ip[1]:
            a1 = np.append( a1, a[j] )
    return a1
# enddef remaining_arr()

# Combinations for choosing p pairs, out of n distinguishable objects:
def comb_pair(  a, b, p, p0, a0 ):
# a[] is a 1D array, storing indices of the balls (objects)
# b[,] is an N*2 2D array, recording all combinations of pair indices
# p is the the current number of pairs to select
# p0 is a fixed number, indicating the total number of pairs to select
# a0 is fixed == range(n), the full indices of n balls (i.e., 0, 1, 2, ..., n-1)

    # For special cases where there is no pair:
    if p0 == -1:
        print()
        return

    n = len(a) # Current total number of balls

    # For the last possible position of selected balls (as stored in b[]),
    # there are the following choices: n, or n-1, or n-2, ..., or m
    for i in range( n-1, 0, -1 ): # i = n-1, n-2, ..., 1
        b[p,0] = a[i]  # Assign the last position to the current b[]
        for j in range( i-1, -1, -1 ): # j = i-1, i-2, ..., 0
            b[p,1] = a[j]  # Assign the last position to the current b[]
            if p > 0: # If it is not the first selection
                # Update the last position as p-1;
                # There are (p-1) pairs left to choose, out of len(a1) balls,
                # where a1[] stores indices of the remaining balls
                a1 = remaining_arr( a, b[p,:] )
                comb_pair( a1, b, p-1, p0, a0 )
            else: # When m == 1, all m0 balls are picked out. 
                 # So, print the result.
                 for k in range( p0+1 ):
                     # If indices start from 1:
                     print( '(%i, %i) ' % ( a0[ b[k,1] ]+1, a0[ b[k,0] ]+1 ), \
                            end='' )
                     # If indices start from 0:
                     #print( '(%i, %i) ' % ( a0[ b[k,0] ], a0[ b[k,1] ] ), \
                     #       end='' )
                     if k == p0:
                         print()

# enddef comb_pair()


# Combinations for choosing m objects and p pairs, at the same time, 
# out of n distinguishable objects
#
def comb_all( a, b1, b2, n1, m, m0, n2, p, p0, a0 ):
# m0 records the total number of balls to select
# p0 records the total number of pairs to select

    # For the last possible position of selected balls (as stored in b1[]),
    # there are the following choices: n1, or n1-1, or n1-2, ..., or m
    for i1 in range( n1-1, m-1, -1 ):  # i1 = n1-1, n1-2, ..., m-1
        if m > -1:
            b1[m] = i1 # Assign the last position to the current b1[]
        if m > 0: # If it is not the first selection
            # Update the last position as i1-1;
            # There are (m-1) balls left to choose;
            # m0 records the total number of balls to select
            comb_all( a, b1, b2, i1, m-1, m0, n2, p, p0, a0 )
        else: # When m == 0, all m0 balls are picked out.
            # Write a 0 to indicate the beginning of a set of m single objects:
            # print( '0 ', end='' )

            for k in range( m0+1 ):
                # If indices start from 1:
                print( '%i ' % (a0[ b1[k] ]+1), end='' )
                # If indices start from 0:
                # print( '%i ' % a0[ b1[k] ], end='' )
            print()

            a_x = np.setdiff1d( a, b1 )
            comb_pair( range( len(a_x) ), b2, p, p0, a_x )

            if m == -1: # Force to end the loop
                break

# enddef comb_all()


# Same as comb_pair() but output to a file:
# Combinations for choosing p pairs, out of n distinguishable objects:
def comb_pair_write(  writer, a, b, p, p0, a0 ):
# a[] is a 1D array, storing indices of the balls (objects)
# b[,] is an N*2 2D array, recording all combinations of pair indices
# p is the the current number of pairs to select
# p0 is a fixed number, indicating the total number of pairs to select
# a0 is fixed == range(n), the full indices of n balls (i.e., 0, 1, 2, ..., n-1)

    # For special cases where there is no pair:
    if p0 == -1:
        return

    n = len(a) # Current total number of balls

    # For the last possible position of selected balls (as stored in b[]),
    # there are the following choices: n, or n-1, or n-2, ..., or m
    for i in range( n-1, 0, -1 ): # i = n-1, n-2, ..., 1
        b[p,0] = a[i]  # Assign the last position to the current b[]
        for j in range( i-1, -1, -1 ): # j = i-1, i-2, ..., 0
            b[p,1] = a[j]  # Assign the last position to the current b[]
            if p > 0: # If it is not the first selection
                # Update the last position as p-1;
                # There are (p-1) pairs left to choose, out of len(a1) balls,
                # where a1[] stores indices of the remaining balls
                a1 = remaining_arr( a, b[p,:] )
                comb_pair_write( writer, a1, b, p-1, p0, a0 )
            else: # When m == 1, all m0 balls are picked out. 
                 # So, print the result.
                 for k in range( p0+1 ):
                     # If indices start from 1:
                     writer.write( ' %i %i' % 
                             ( a0[ b[k,1] ]+1, a0[ b[k,0] ]+1 ) )
                     if k == p0:
                         writer.write( '\n' )

# enddef comb_pair_write()


# Same as comb_all() but output to a file:
# Combinations for choosing m objects and p pairs, at the same time, 
# out of n distinguishable objects
#
def comb_all_write( writer, a, b1, b2, n1, m, m0, n2, p, p0, a0 ):
# m0 records the total number of balls to select
# p0 records the total number of pairs to select

    # For the last possible position of selected balls (as stored in b1[]),
    # there are the following choices: n1, or n1-1, or n1-2, ..., or m
    for i1 in range( n1-1, m-1, -1 ):  # i = n1-1, n1-2, ..., m-1
        if m > -1:
            b1[m] = i1 # Assign the last position to the current b1[]
        if m > 0: # If it is not the first selection
            # Update the last position as i1-1;
            # There are (m-1) balls left to choose;
            # m0 records the total number of balls to select
            comb_all_write( writer, a, b1, b2, i1, m-1, m0, n2, p, p0, a0 )
        else: # When m == 0, all m0 balls are picked out.
            # Write a 0 to indicate the beginning of a set of m single objects:
            writer.write( '0' )

            for k in range( m0+1 ):
                # If indices start from 1:
                writer.write( ' %i' % (a0[ b1[k] ]+1) )
                # If indices start from 0:
                # print( '%i ' % a0[ b1[k] ], end='' )
            writer.write( '\n' )

            a_x = np.setdiff1d( a, b1 )
            comb_pair_write( writer, range( len(a_x) ), b2, p, p0, a_x )

            if m == -1: # Force to end the loop
                break

# enddef comb_all_write()


# Generate all possible lewis structures and output to a file:
def lewis_write( NAt, NE ):
    outp = 'LEWIS_%ic_%ie.dat' % ( NAt, NE )

    # Verify even number of electrons (closed-shell system):
    if NE %2 != 0:
        raise ValueError( 'Number of electrons must be an even number!' )

    # Number of electron pairs that may be either LPs or BDs:
    nPair = NE // 2

    with open( outp, 'w' ) as writer:
        for nLP in range( nPair+1 ):
            allComb_write( writer, NAt, nLP, nPair - nLP )

# enddef lewis_write()


# Read one Lewis structure from a file that contains Lewis structures, 
# generated by lewis_write()
def lewis_readone( reader, LS, maxNLP ):
    if LS.stat == 1: # EOF
        return

    for line in reader:
        # print( line )
        nf = []  # List for storing numbers of each text field
        for sf in line.split():
            nf.append( int( sf ) )

        # Skip blank lines:
        if len( nf ) == 0:
            continue

        # If a new combination of LPs begins:
        if nf[0] == 0:
            LS.nLP = len( nf ) - 1
            LS.nBD = LS.nPair - LS.nLP
            # Stop if nLP excceeds maxNLP:
            if LS.nLP > maxNLP:
                LS.stat = 1
                return LS

            LS.LP = np.array( nf[ 1: ], dtype='i' )
            # print( LS.LP )

            # In case that there are only LPs and no BDs:
            if LS.nLP == LS.nPair:
                LS.BD = []
                return LS
        else:
            LS.BD = np.array( nf, dtype='i' ).reshape( (LS.nBD,2) )
            # print( LS.BD )
            return LS

    LS.stat = 1 # Reached EOF

    return LS
# enddef lewis_readone()


# Read a bunch of Lewis strcuture with a maximum number of LPs

def lewis_read( NAt, NE, maxNLP=np.Inf ):
    inp = 'LEWIS_%ic_%ie.dat' % (NAt, NE)
    # Check file does not exit, then make the enumeration:
    if not os.path.exists( inp ):
        print( 'File', inp, 'not found' )
        print( 'Generating', inp, 'for Lewis structures ...' )
        lewis_write( NAt, NE )

    LS = LewisStructure( NAt, NE  )
    print( 'Reading', inp, 'for Lewis structures ...' )
    LP = [];
    BD = [];
    with open( inp, 'r' ) as f:
        while True:
            LS = lewis_readone( f, LS, maxNLP )

            if LS.stat == 1:
                break

            #print( 'LS ', iLS+1, LS1.LP, ' | ', LS1.BD )
            #print( 'LS %i ' % (iLS+1), *LS.LP, ' | ', *LS.BD )
            #print( lewis_str( LS.LP, LS.BD ) )
            LP.append( LS.LP )
            BD.append( LS.BD )

    return ( LP, BD )

# enddef lewis_read()


# String of a given Lewis structure:
def lewis_str( LP, BD, RA=[], RB=[] ):
    # Initialize:
    s = ''

    # Alpha radical sites:
    if len( RA ) > 0:
        s += '^ '.join( map( str, RA ) ) + '^ '

    # Beta radical sites:
    if len( RB ) > 0:
        s += 'v '.join( map( str, RB ) ) + 'v '

    # Lone pairs:
    if len( LP ) > 0:
        s += ': '.join( map( str, LP ) ) + ': '

    if len( BD ) == 0:
        return s

    i = 1
    BD1 = BD.copy()
    BD1 = BD1.reshape( BD.shape[0]*BD.shape[1] )
    nBD = len( BD1 )
    for k in BD1:
        s += '%i' % k
        if i % 2 == 1:
            s += '-'
        elif i < nBD:
            s += ' '
        i += 1

    return s
# enddef lewis_str()


# Parse the strings of a given set of Lewis structures:
def str2Lewis( strList ):
    LP = []
    BD = []
    for s in strList:
        lps = []
        bds = []
        for f in s.split():
            # Lone pair:
            if ':' in f:
                n = f.replace( ':', '' )
                #if len(n) != 1 or int(n) < 0:
                if int(n) < 0:
                    raise ValueError( 'Invalid Lewis structure:', s )
                lps.append( np.array( n, dtype='i' ) )
            # Bond:
            elif '-' in f:
                b = f.split( '-' )
                if len(b) != 2 or int(b[0]) < 0 or int(b[1]) < 0:
                    raise ValueError( 'Invalid Lewis structure:', s )
                bds.append( np.array( b, dtype='i' ) )
            # Invalid:
            else:
                raise ValueError( 'Invalid Lewis structure:', s )
        LP.append( np.array( lps, dtype='i' ) )
        BD.append( np.array( bds, dtype='i' ) )

    return LP, BD
# enddef str2Lewis()

