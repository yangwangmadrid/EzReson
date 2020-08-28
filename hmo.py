#Simple Hueckel Molecular Orbital (HMO) calculations
#
# Created on Apr 27, 2020
# Written by Yang Wang (yangwang@yzu.edu.cn)

import numpy as np
from topo import *
from readxyz import *

class huckel:
    def __init__( self, NAt, NE, NEa, NEb, C, E, Etot, gap, chg, spin, \
            A, D, Da, Db, occ, occ_a, occ_b, BOs, BO, F, WBI ):
        self.NAt = NAt      #  Number of atoms/MOs
        self.NE = NE        # Number of electrons
        self.NEa = NEa      # Number alpha of electrons
        self.NEb = NEb      # Number beta of electrons
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
        #print( inp )
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
        else:
            raise ValueError( 'Invalid data type', type(inp), 
                              'of argument inp for function hmo()' )
    else:
        raise ValueError( 'Invalid data type', type(inp), 
                          'of argument inp for function hmo()' )


    NAt = len( A ) # Number of atoms or MOs
    NE = NAt - q      # Number of electrons
    NEa = NE // 2   # Number of alpha electrons
    NEb = NE - NEa  # Number of beta electrons

    E, C = np.linalg.eigh( A )
    E = -E # in units of -beta
    # Sort E from negative to positive values
    ix_sort = np.argsort( E )
    E = E[ ix_sort ]
    C = C[ :, ix_sort ]

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

    return huckel( NAt, NE, NEa, NEb, C, E, Etot, gap, chg, spin, \
            A, D, Da, Db, occ, occ_a, occ_b, BOs, BO, F, WBI )
# enddef hmo()
