# Simple Hueckel Molecular Orbital (HMO) calculations
#
# Created on Apr 27, 2020
# Written by Yang Wang (yangwang@yzu.edu.cn)
#
#
#  Updates:
#
#    May 10, 2022:
#    - Extended Huckel calculations to open-shell systems
#
#    Oct 2, 2021:
#    - Supported *.gjf file as the input

import numpy as np
from topo import *
from readxyz import *

class huckel:
    def __init__( self, xyz, NAt, NE, NEa, NEb, mult, C, E, Etot, gap, chg, \
            spin, A, D, Da, Db, occ, occ_a, occ_b, BOs, BO, F, WBI ):
        self.xyz = xyz      # Cartesian coordinates
        self.NAt = NAt      # Number of atoms/MOs
        self.NE = NE        # Number of electrons
        self.NEa = NEa      # Number alpha of electrons
        self.NEb = NEb      # Number beta of electrons
        self.mult = mult    # Spin multiplicity
        self.C = C          # Coefficients of MO, C(iat,iMO)
        self.E = E          # Eigen energies in units of -beta
        self.Etot = Etot    # Total pi energy in units of -beta
        self.gap = gap      # HOMO-LUMO gap in units of -beta
        self.chg = chg      # Atomic charges
        self.spin = spin    # Atomic spin densities
        self.A = A          # Adjacency matrix
        self.D = D          # Density matrix of all electrons
        self.Da = Da        # Density matrix of alpha electrons
        self.Db = Db        # Density matrix of beta electrons
        self.occ = occ      # Occupanies of MOs
        self.occ_a = occ_a  # Occupanies of alpha MOs
        self.occ_b = occ_b  # Occupanies of beta MOs
        self.BOs = BOs      # Coulson BOs for ALL pairs of atoms no matter 
                            # they are bonded or not
        self.BO = BO        # Coulson bond orders only for bonded pairs:
                            #   BO(:,1): index of 1st atom
                            #   BO(:,2): index of 2nd atom  ( BO(:,1)<BO(:,2) )
                            #   BO(:,3): bond order in units of -beta
        self.F = F          # Coulson free valences
        self.WBI = WBI      # Wiberg bond index matrix with nondiagonal entries 
                            # being WBIs and diagonal entries being valence
#==============================================================================
# endclass huckel


# Perform simple Hueckel Molecular Orbital (HMO) calculations:
#
# inp: There are many options for this argument:
#      1) <string>: Input filename containing Cartesian coordinates of atoms
#      2) <NAt*3 ndarray, dtype='f'>: Cartesian coordinates of atoms
#      3) <NAt*NAt ndarray, dtype='i'>: Adjacency matrix
# q:   Net charge of the system (an integer)
def hmo( inp, q=0 ):

    if type(inp) is str:
        if inp[-4:] == '.gjf':
            elem, xyz = readgjf( inp )
        else:
            elem, xyz = readxyz( inp )
        # Remove hydrogens:
        xyz, _ = xyzNonH( elem, xyz )
        A = adjMat( xyz )
    elif type(inp) is tuple:
        elem, xyz = inp
        # Remove hydrogens:
        xyz, _ = xyzNonH( elem, xyz )
        A = adjMat( xyz )
    elif type(inp) is np.ndarray:
        if inp.dtype is np.dtype( 'float' ):
            xyz = inp
            A = adjMat( xyz )
        elif inp.dtype is np.dtype( 'int' ):
            A = inp
            xyz = []
        else:
            raise ValueError( 'Invalid data type', type(inp), 
                              'of argument inp for function hmo()' )
    else:
        raise ValueError( 'Invalid data type', type(inp), 
                          'of argument inp for function hmo()' )


    NAt = len( A )       # Number of atoms or MOs
    NE = NAt - q         # Total number of electrons

    E, C = np.linalg.eigh( A )
    E = -E # in units of -beta
    # Sort E from negative to positive values
    ix_sort = np.argsort( E )
    E = E[ ix_sort ]
    C = C[ :, ix_sort ]

    # Applying Hund's rule and Aufbau principle to fill the levels with
    # alpha and beta electrons:
    # (1) Degeneracies of of levels:
    #E1 = E
    #Deg = np.ones( NAt, dtype='i' )
    E1 = []
    Deg = []
    j = 0
    for k in range( NAt ):
        if k > 0 and abs( E1[j-1] - E[k] ) < 1E-13:
            Deg[j-1] += 1
        else:
            E1.append(  E[k] )
            Deg.append( 1 )
            j += 1

    # (2) Fill the electrons:
    occ_a = np.zeros( NAt, dtype='i' ) # alpha MO occupancies
    occ_b = np.zeros( NAt, dtype='i' ) # beta MO occupancies
    counter = 0
    ia = 0
    ib = 0
    for k in range( len(Deg) ):
        deg_k = Deg[k]
        # Fill with alpha electrons:
        for j in range( deg_k ):
            if counter == NE: # No eletron left
                break
            counter += 1
            occ_a[ ia ] = 1
            ia += 1
        # Fill with beta electrons:
        for j in range( deg_k ):
            if counter == NE: # No eletron left
                break
            counter += 1
            occ_b[ ib ] = 1
            ib += 1

    NEa = sum( occ_a )   # Number of alpha electrons
    NEb = sum( occ_b )   # Number of beta electrons
    mult = NEa - NEb + 1 # Spin multiplicity

    # MO occupancies:
    occ_a = np.r_[ np.ones( NEa ), np.zeros(NAt-NEa) ]
    occ_b = np.r_[ np.ones( NEb ), np.zeros(NAt-NEb) ]
    occ = occ_a + occ_b
    occ_spin = occ_a - occ_b

    # Total energy:
    Etot = ( occ * E ).sum()

    # HOMO-LUMO gap:
    if NEa != NEb:  # open-shell
        gap = 0
    else:
        gap = E[ NEa ] - E[ NEa-1 ]

    # Density matrix:
    Da = C[ :, range(NEa) ] @ C[ :, range(NEa) ].T
    Db = C[ :, range(NEb) ] @ C[ :, range(NEb) ].T
    D = Da + Db

    # Atomic charges:
    chg = 1 - np.diag( D )

    # Atomic spin charges:
    spin = np.diag( Da - Db )

    # Coulson bond order matrix == density matrix of HMOs:
    BOs = D

    # Coulson bond orders and free valences only for connected bonds:
    BO = []
    F = np.zeros( NAt ) # initialize
    for i1 in range( NAt ):
        for i2 in range( i1+1, NAt ):
            if A[ i1, i2 ] == 1:
                BO.append( [ i1+1, i2+1, BOs[ i1, i2 ] ] )
                F[ i1 ] += BOs[ i1, i2 ]
                F[ i2 ] += BOs[ i1, i2 ]
    BO = np.array( BO )
    F = 1.73205080756887719318 - F # (3+sqrt(3)) - (F+3), with three sigma-bonds

    # Wiberg bond indices:
    WBI = D * D
    for i in range( NAt ):
        WBI[ i, i ] = 0.
        WBI[ i, i ] = WBI[ i, : ].sum()

    return huckel( xyz, NAt, NE, NEa, NEb, mult, C, E, Etot, gap, chg, spin, \
            A, D, Da, Db, occ, occ_a, occ_b, BOs, BO, F, WBI )
# enddef hmo()
