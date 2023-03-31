# Read NBO's disk files
#
#  Created on Feb, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#

import re
import numpy as np

class NAOInformation:
  def __init__( self, NNAO, NNAODep, elem, atIx, naoIx, nl, lm, typ,
          spin, NAt, NE, NEa, NEb, NBas, NMO 
          ):
     self.NNAO = NNAO      # Number of NAOs
     self.NNAODep = NNAODep   # Number of linear dependent NAOs
     self.elem = elem      # Element string of atom
     self.atIx = atIx      # Atomic index of NAO, starting from 1 (not 0)
     self.naoIx = naoIx    # NAO indices for each atom, starting from 1 (not 0)
     self.nl = nl          # Principal and angular numbers of NAO
     self.lm = lm          # Angular and magnetic quantum numbers of NAO
     self.typ = typ        # Type of NAO
     self.spin = spin      # Spin --> 0: restricted; 1: unrestricted
     self.NAt = NAt        # Number of atoms
     self.NE = NE          # Number of electrons
     self.NEa = NEa        # Number of alpha electrons
     self.NEb = NEb        # Number of beta electrons
     self.NBas = NBas  # Number of basis functions
     self.NMO = NMO        # Number of molecular orbitals
#==============================================================================
# endclass NAOInformation

#==============================================================================
# Read NBO's disk file for symmetric matrix by keyword 'MATKEY=W',
# such as density matrix, Fock matrix and overlap matrix.
#
# Return values:
#   For spin-restricted case, return a matrix (ndarray): M
#   For spin-unrestricted case, return a tuple of matrices: ( Ma, Mb )
#
# Note: For overlap matrix, there is no spin distinguished
#
def readNBOSymmMat( inputFileName, numOrb ):
    # Initialize the matrix:
    Ma = np.empty( (numOrb, numOrb) )

    pat_ruler = re.compile( r'^\s*-+\n$' )
    pat_alphaspin = re.compile( r'^\s*ALPHA\s+SPIN' )
    pat_betaspin = re.compile( r'^\s*BETA\s+SPIN' )

    flag = 0
    iRow = 0
    iCol = 0
    iline = 0
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            # Locate the ruler line:
            if pat_ruler.match( line ):
                flag = 1
                continue

            # Check if it is a spin-unrestricted calculation:
            if flag == 1:
                if pat_alphaspin.match( line ):
                    spin = 1
                elif pat_betaspin.match( line ):
                    spin = 2
                    Mb = empty( (numOrb, numOrb) )
                else:
                    spin = 0

                flag += 1
                # print( 'spin =', spin )
                if spin > 0:
                    continue

            # Parsing the data:
            if flag != 2:
                continue

            nf = 0
            for x in line.split():
                # print( '%i, %i' % (iRow, iCol) )
                nf += 1
                if spin <= 1: # alpha or all
                    Ma[iRow][iCol] = float(x)
                    Ma[iCol][iRow] = Ma[iRow][iCol] # symmetric matrix
                else: # beta
                    Mb[iRow][iCol] = float(x)
                    Mb[iCol][iRow] = Mb[iRow][iCol] # symmetric matrix
                iRow += 1

                # Check if a column is complete:
                if iRow == iCol+1:
                    iRow = 0
                    iCol += 1
                elif iRow > numOrb:
                    print( 'More data than expected ({0}) at line {1}:'.format( 
                            numOrb, iline ) )
                    print( line )
                    exit()

            if nf == 0:
                print( 'Empty data at line {0}:'.format( iline ) )
                print( line )
                exit()

            # Check if all alpha spin is done:
            if spin <= 1 and iCol == numOrb:
                if spin == 0:
                    break
                else: # beta
                    flag = 1
                    iCol = 0

            # Check if all beta spin is done:
            if spin == 2 and iCol > numOrb:
                print( 'More data than expected ({0}) at line {1}:'.format( \
                        numOrb, iline ) )
                print( line )
                exit()

    if spin <= 1: # spin-restricted
        M = Ma
    else:         # spin-unrestricted
        M = ( Ma, Mb )
    
    return M
#==============================================================================
# enddef readNBOSymmMat()


#==============================================================================
# Read NBO's disk file for general matrix by keyword 'MATKEY=W',
# such as AONAO and AOMO coefficient matrices.
#
# Return values:
#   For spin-independent case, return a matrix (ndarray): M
#   For spin-dependent case, return a tuple of matrices: Ma, Mb
#
# Note: For overlap matrix, there is no spin distinguished
#
def readNBOMat( inputFileName, numBas, numOrb ):
    # Initialize the matrix:
    Ma = np.empty( (numBas, numOrb) )

    pat_ruler = re.compile( r'^\s*-+\n$' )
    pat_alphaspin = re.compile( r'^\s*ALPHA\s+SPIN' )
    pat_betaspin = re.compile( r'^\s*BETA\s+SPIN' )

    flag = 0
    iBas = 0
    iOrb = 0
    iline = 0
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            # Locate the ruler line:
            if pat_ruler.match( line ):
                flag = 1
                continue

            # Check if it is a spin-dependent calculation:
            if flag == 1:
                if pat_alphaspin.match( line ):
                    spin = 1
                elif pat_betaspin.match( line ):
                    spin = 2
                    Mb = np.empty( (numBas, numOrb) )
                else:
                    spin = 0

                flag += 1
                # print( 'spin =', spin )
                if spin > 0:
                    continue

            # Parsing the data:
            if flag != 2:
                continue

            nf = 0
            for x in line.split():
                # print( '%i, %i' % (iBas, iOrb) )
                nf += 1
                if spin <= 1: # alpha or all
                    Ma[iBas][iOrb] = float(x)
                else: # beta
                    Mb[iBas][iOrb] = float(x)
                iBas += 1

                # Check if an MO is complete:
                if iBas == numBas:
                    iBas = 0
                    iOrb += 1
                elif iBas > numBas:
                    print( 'More data than expected ({0}) at line {1}:'.format( 
                            numBas, iline ) )
                    print( line )
                    exit()

            if nf == 0:
                print( 'Empty data at line {0}:'.format( iline ) )
                print( line )
                exit()

            # Check if all alpha spin is done:
            if spin <= 1 and iOrb == numOrb:
                if spin == 0:
                    break
                else: # beta
                    flag = 1
                    iOrb = 0

            # Check if all beta spin is done:
            if spin == 2 and iOrb > numOrb:
                print( 'More data than expected ({0}) at line {1}:'.format( \
                        numOrb, iline ) )
                print( line )
                exit()

    if spin <= 1: # spin-restricted
        M = Ma
    else:         # spin-unrestricted
        M = ( Ma, Mb )
    
    return M
#==============================================================================
# enddef readNBOMat()


#==============================================================================
# Read NAO information from Gaussian's output file
#
# Return:
#   NAOInformation NAOInfo
#
def readNAOInfo( inputFileName ):
    pat_basis = re.compile( 
            '^\s*\d+\s*basis\s+functions,*\s*\d+\s*primitive\s+'
            'gaussians,*\s*\d+\s*cartesian basis functions\s*'
            )

    pat_electron = re.compile( 
            '^\s*\d+\s*alpha\s+electrons\s*\d+\s*beta\s+electrons\s*'
            )

    pat_atom = re.compile( 
            '^\s*NAtoms=\s+\d+\s+.*=\s+\d+.*'
            )

    pat_indpBasis = re.compile( 
            '^\s*NBsUse=\s+\d+.*EigRej=\s*[\d|-]+.*'
            )

    pat_nao_header = re.compile( 
            '^\s*NAO\s+Atom\s+No\s+lang\s+Type.*\s+Occupancy\s*'
            )

    pat_ruler = re.compile( '^\s*'+'-'*10+'\s*' )

    pat_nao = re.compile(
             '^\s*\d+\s*\w+\s*\d+\s*[\w|\(|\)]+\s*\w+\(\s*\w+\)\s*'
             '\d*\.\d*.*'
            )

    # Initialize lists of NAO info:
    atIx = [] # Atomic index of NAO
    nl = []   # Principal and angular numbers of NAO
    lm = []   # Angular and magnetic quantum numbers of NAO
    typ = []  # Type of NAO
    elem = [] # Element string of atom

    iNAO_last = 0 # for determing the end of parsing NAO info.
    flag = 0
    iline = 0
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            # Basis functions line:
            if pat_basis.match( line ):
                NBas = int( (line.split())[0] )
                #print( line )
                #print( NBas )
                continue

            # Numbers of electrons line:
            if pat_electron.match( line ):
                NEa = int( (line.split())[0] )
                NEb = int( (line.split())[3] )
                NE = NEa + NEb
                #print( line )
                #print( "%i %i %i" % (NEa, NEb, NE) )
                continue

            # Numbers of atoms line:
            if pat_atom.match( line ):
                NAt = int( (line.split())[1] )
                #print( line )
                #print( "%i" % (NAt) )
                continue

            # Numbers of independent basis functions line:
            if pat_indpBasis.match( line ):
                NMO = int( (line.split())[1] )
                #print( line )
                #print( "%i" % (NMO) )
                continue

            # Header of NAO info line:
            if flag == 0 and pat_nao_header.match( line ):
                flag = 1
                pat_spin = re.compile( '^.*Spin\s*' )
                if pat_spin.match( line ):
                    spin = 1
                else:
                    spin = 0
                #print( line )
                #print( 'Spin = %i' % spin )
                continue

            # Ruler preceeding NAO info line:
            if flag == 1 and pat_ruler.match( line ):
                flag = 2
                #print( line )
                continue

            # NAO info line:
            if flag == 2 and pat_nao.match( line ):
                #print( line )
                # Note that there may appear two '(**)' fields, as the first
                # may be the angular/magnetic orbital like 'f(c1)', and the
                # second the orbital type like 'Ryd( 5g)'
                field = line.rsplit( '(', maxsplit=1 )
                field1 = field[0];
                field = field[1].rsplit( ')', maxsplit=1 )
                field2 = field[0];
                
                s1 = field1.split()
                #-- Index of NAO: --
                iNAO = int( s1[0] ) - 1
                #print( 'iNAO %i' % iNAO )
                # Check if reading of total NAO info is finished for
                # spin-unrestricted cases:
                if iNAO < iNAO_last:
                    break
                iNAO_last = iNAO

                #-- Angular and magnetic quantum numbers of NAO: --
                lm.append( s1[-2] )
                #print( 'lm %s' % lm[iNAO] )
                #-- Type of NAO: --
                typ.append( s1[-1] )
                #print( 'type %s' % typ[iNAO] )
                #-- Element string of atom: --
                elem.append( re.sub( '\d', '', s1[1] ) )
                #print( 'elem %s' % elem[iNAO] )
                #-- Atomic index of the current NAO: --
                atIx.append( int( re.sub( '[A-z]', '', s1[-3] ) ) )
                #print( 'atIx %i' % atIx[iNAO] )

                #-- Principal and angular numbers of NAO: --
                nl.append( field2 )
                continue

    # Number of NAOs:
    NNAO = len( elem )
    #print( 'NNAO %i' % NNAO )
    # Number of discarded basis functions due to lindear dependency:
    NNAODep = NBas - NNAO
    if NNAODep < 0:
        raise ValueError( 'Number of NAOs (%i) is greater than'
                'number of basis functions (%i)' % (NNAO, NBas) )

    # Check if number of atoms is consistent:
    # print( atIx )
    if NAt != max( atIx ):
        raise ValueError( 'Number of atoms (%i) is not consistet with'
                'NAt in NAO info (%i)' % (NAt, max( atIx )) )

    # NAO indices for each atom:
    naoIx = []
    iat = 1
    istart = 0
    ix = np.arange( 1, NNAO+1 )
    for i in ix:
        if atIx[i-1] > iat:
            naoIx.append( ix[ istart : i-1 ] )
            iat += 1
            istart = i-1
    naoIx.append( ix[ istart : NNAO+1 ] )
    #print( naoIx )

    #---- Collection all NAO info: ----
    return NAOInformation( NNAO, NNAODep, elem, atIx, naoIx, nl, lm, typ,
          spin, NAt, NE, NEa, NEb, NBas, NMO
          )
#==============================================================================
# enddef readNAOInfo()
