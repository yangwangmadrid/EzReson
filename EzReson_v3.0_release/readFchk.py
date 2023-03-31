# Read Gaussian's fchk file for infos such as orbital energies and coefficients
#
#  Created on Feb, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#

global DEBUG_MODE
#DEBUG_MODE = True
DEBUG_MODE = False

# Bohr --> Angstrom, see https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
global BOHR_TO_ANG
BOHR_TO_ANG = 0.529177210903

import re
import numpy as np

################################################################################
if DEBUG_MODE:
    print( '*'*80 )
    print( '    DEBUGGING MODE' )
    print( '        (Set DEBUG_MODE = False to swith off)' )
    print( '*'*80 )

def debug( obj, **options ):
    if DEBUG_MODE:
        print( obj, end=options.get('end')  )
################################################################################


class FchkInformation:
  # Note: both E and C are tuples
  #   For spin-unrestricted calculations, both E and C are tuples:
  #     E = (Ea, Eb)
  #     C = (Ca, Cb)
  #   For spin-restricted calculations, both E and C are np arrays (matrices):
  #     E = E
  #     C = C
  #
  def __init__( self, chg, mult, NE, NEa, NEb, NAt, Z, Zeff, M, xyz, spin, 
          NBas, NBasIndp, atIx, aoIx, E_tot, E_SCF, E, C, QMull
          ):
     self.chg = chg               # Total charge of the moelcule
     self.mult = mult             # Spin multiplicity of the moelcule
     self.NE = NE                 # Number of electrons
     self.NEa = NEa               # Number of alpha electrons
     self.NEb = NEb               # Number of beta electrons
     self.NAt = NAt               # Number of atoms
     self.Z = Z                   # Atomic numbers of elements
     self.Zeff = Zeff             # Effective atomic numbers (when ECP is used)
     self.M = M                   # Atomic masses (of the most abundant isotope)
     self.xyz = xyz * BOHR_TO_ANG # Cartesian coordinates of atoms
     self.spin = spin             # Spin --> 0: restricted; 1: unrestricted
     self.NBas = NBas             # Number of basis functions
     self.NBasIndp = NBasIndp     # Number of independent basis functions
     self.atIx = atIx             # Atomic index of AO
     self.aoIx = aoIx             # AO indices for each atom
     self.E_tot = E_tot           # Total energy
     self.E_SCF = E_SCF           # Total SCF energy
     self.E = E                   # Energies of CMOs
     self.C = C                   # Coefficients of CMOs
     self.QMull = QMull           # Mulliken charges
#==============================================================================
# endclass FchkInformation


#==============================================================================
# Read Gaussian's fchk file
#
# Return:
#   FchkInformation FchkInfo
#
def readFchk( inputFileName ):
    pat_NAt = re.compile( 
            '^\s*Number\s+of\s+atoms\s+[A-Z]\s+\d+.*'
            )

    pat_chg = re.compile( 
            '^\s*Charge\s+[A-Z]\s+-*\d+.*'  # Possibly having '-' sign
            )

    pat_mult = re.compile( 
            '^\s*Multiplicity\s+[A-Z]\s+\d+.*'
            )

    pat_NE = re.compile( 
            '^\s*Number\s+of\s+electrons\s+[A-Z]\s+\d+.*'
            )

    pat_NEa = re.compile( 
            '^\s*Number\s+of\s+alpha\s+electrons\s+[A-Z]\s+\d+.*'
            )

    pat_NEb = re.compile( 
            '^\s*Number\s+of\s+beta\s+electrons\s+[A-Z]\s+\d+.*'
            )

    pat_NBas = re.compile( 
            '^\s*Number\s+of\s+basis\s+functions\s+[A-Z]\s+\d+.*'
            )

    pat_NBasIndp = re.compile( 
            '^\s*Number\s+of\s+independent\s+functions\s+[A-Z]\s+\d+.*'
            )

    pat_Z = re.compile( 
            '^\s*Atomic\s+numbers\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_Z = 0
    i_Z = 0

    pat_Zeff = re.compile( 
            '^\s*Nuclear\s+charges\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_Zeff = 0
    i_Zeff = 0

    pat_xyz = re.compile( 
            '^\s*Current\s+cartesian\s+coordinates\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_xyz = 0
    i_xyz = 0

    pat_M = re.compile( 
            '^\s*Real\s+atomic\s+weights\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_M = 0
    i_M = 0

    pat_shTyp = re.compile( 
            '^\s*Shell\s+types\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_shTyp = 0
    i_shTyp = 0

    pat_shTyp = re.compile( 
            '^\s*Shell\s+types\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_shTyp = 0
    i_shTyp = 0

    pat_shAtMap = re.compile( 
            '^\s*Shell\s+to\s+atom\s+map\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_shAtMap = 0
    i_shAtMap = 0

    pat_E_SCF = re.compile( 
            '^\s*SCF\s+Energy\s+[A-Z]\s+-*\d+\.*\d*[E|D]*\d*.*'
            )

    pat_E_tot = re.compile( 
            '^\s*Total\s+Energy\s+[A-Z]\s+-*\d+\.*\d*[E|D]*\d*.*'
            )

    pat_Ea = re.compile( 
            '^\s*Alpha\s+Orbital\s+Energies\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_Ea = 0
    i_Ea = 0

    pat_Eb = re.compile( 
            '^\s*Beta\s+Orbital\s+Energies\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_Eb = 0
    i_Eb = 0

    pat_Ca = re.compile( 
            '^\s*Alpha\s+MO\s+coefficients\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_Ca = 0
    i_Ca = 0

    pat_Cb = re.compile( 
            '^\s*Beta\s+MO\s+coefficients\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_Cb = 0
    i_Cb = 0

    pat_QMull = re.compile( 
            '^\s*Mulliken\s+Charges\s+[A-Z]\s+[A-Z]+=\s+\d+.*'
            )
    n_QMull = 0
    i_QMull = 0

    spin = 0; # Initialize: spin-restricted
    iline = 0
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            # Number of atoms:
            if pat_NAt.match( line ):
                NAt = int( (line.split())[-1] )
                debug( 'NAt = %i' % NAt )
                continue

            # Charge:
            if pat_chg.match( line ):
                chg = int( (line.split())[-1] )
                debug( 'chg = %i' % chg )
                continue

            # Multiplicity:
            if pat_mult.match( line ):
                mult = int( (line.split())[-1] )
                debug( 'mult = %i' % mult )
                continue

            # Number of electrons:
            if pat_NE.match( line ):
                NE = int( (line.split())[-1] )
                debug( 'NE = %i' % NE )
                continue

            # Number of alpha electrons:
            if pat_NEa.match( line ):
                NEa = int( (line.split())[-1] )
                debug( 'NEa = %i' % NEa )
                continue

            # Number of beta electrons:
            if pat_NEb.match( line ):
                NEb = int( (line.split())[-1] )
                debug( 'NEb = %i' % NEb )
                continue

            # Number of basis functions:
            if pat_NBas.match( line ):
                NBas = int( (line.split())[-1] )
                debug( 'NBas = %i' % NBas )
                continue

            # Number of independent functions:
            if pat_NBasIndp.match( line ):
                NBasIndp = int( (line.split())[-1] )
                debug( 'NBasIndp = %i' % NBasIndp )
                continue

            # Atomic numbers
            if n_Z == 0 and pat_Z.match( line ):
                n_Z = int( (line.split())[-1] )
                Z = np.zeros( n_Z, dtype='i' )
                continue
            elif n_Z > 0:
                for x in line.split():
                    Z[i_Z] = int( x )
                    i_Z += 1
                # Check if completed:
                if i_Z == n_Z:
                    n_Z = 0
                    debug( 'Z = ', end='' )
                    debug( Z )
                continue
                
            # Nuclear charges
            if n_Zeff == 0 and pat_Zeff.match( line ):
                n_Zeff = int( (line.split())[-1] )
                Zeff = np.zeros( n_Zeff, dtype='d' )
                continue
            elif n_Zeff > 0:
                for x in line.split():
                    Zeff[i_Zeff] = float( x )
                    i_Zeff += 1
                # Check if completed:
                if i_Zeff == n_Zeff:
                    n_Zeff = 0
                    debug( 'Zeff = ', end='' )
                    debug( Zeff )
                continue
                
            # Cartesian coordinates
            if n_xyz == 0 and pat_xyz.match( line ):
                n_xyz = int( (line.split())[-1] )
                # Check consistentcy with NAt:
                if n_xyz != NAt*3:
                    raise ValueError( 'Number of atoms (%i) is not consistent'
                            'with the number of Cartesian coordinates (%i)' %
                            ( NAt, n_xyz ) )
                xyz = np.zeros( n_xyz, dtype='d' )
                continue
            elif n_xyz > 0:
                for x in line.split():
                    xyz[i_xyz] = float( x )
                    i_xyz += 1
                # Check if completed:
                if i_xyz == n_xyz:
                    xyz = xyz.reshape( NAt, 3 )
                    n_xyz = 0
                    debug( 'xyz = ', end='' )
                    debug( xyz )
                continue
                
            # Atomic masses
            if n_M == 0 and pat_M.match( line ):
                n_M = int( (line.split())[-1] )
                M = np.zeros( n_M, dtype='d' )
                continue
            elif n_M > 0:
                for x in line.split():
                    M[i_M] = float( x )
                    i_M += 1
                # Check if completed:
                if i_M == n_M:
                    n_M = 0
                    debug( 'M = ', end='' )
                    debug( M )
                continue

            # Shell types
            if n_shTyp == 0 and pat_shTyp.match( line ):
                n_shTyp = int( (line.split())[-1] )
                shTyp = np.zeros( n_shTyp, dtype='i' )
                continue
            elif n_shTyp > 0:
                for x in line.split():
                    shTyp[i_shTyp] = float( x )
                    i_shTyp += 1
                # Check if completed:
                if i_shTyp == n_shTyp:
                    n_shTyp = 0
                    debug( 'shTyp = ', end='' )
                    debug( shTyp )
                continue

            # Shell to atom map
            if n_shAtMap == 0 and pat_shAtMap.match( line ):
                n_shAtMap = int( (line.split())[-1] )
                shAtMap = np.zeros( n_shAtMap, dtype='i' )
                continue
            elif n_shAtMap > 0:
                for x in line.split():
                    shAtMap[i_shAtMap] = float( x )
                    i_shAtMap += 1
                # Check if completed:
                if i_shAtMap == n_shAtMap:
                    n_shAtMap = 0
                    debug( 'shAtMap = ', end='' )
                    debug( shAtMap )
                continue

            # SCF energy:
            if pat_E_SCF.match( line ):
                E_SCF = float( (line.split())[-1] )
                debug( 'E_SCF = %.16f' % E_SCF )
                continue

            # Total energy:
            if pat_E_tot.match( line ):
                E_tot = float( (line.split())[-1] )
                debug( 'E_tot = %.16f' % E_tot )
                continue

            # Alpha Orbital Energies
            if n_Ea == 0 and pat_Ea.match( line ):
                n_Ea = int( (line.split())[-1] )
                Ea = np.zeros( n_Ea, dtype='d' )
                continue
            elif n_Ea > 0:
                for x in line.split():
                    Ea[i_Ea] = float( x )
                    i_Ea += 1
                # Check if completed:
                if i_Ea == n_Ea:
                    n_Ea = 0
                    debug( 'Ea = ', end='' )
                    debug( Ea )
                continue

            # Beta Orbital Energies
            if n_Eb == 0 and pat_Eb.match( line ):
                spin = 1; # spin-unrestricted
                debug( 'spin = ', end='' )
                debug( spin )
                n_Eb = int( (line.split())[-1] )
                Eb = np.zeros( n_Eb, dtype='d' )
                continue
            elif n_Eb > 0:
                for x in line.split():
                    Eb[i_Eb] = float( x )
                    i_Eb += 1
                # Check if completed:
                if i_Eb == n_Eb:
                    n_Eb = 0
                    debug( 'Eb = ', end='' )
                    debug( Eb )
                continue

            # Alpha MO coefficients
            if n_Ca == 0 and pat_Ca.match( line ):
                n_Ca = int( (line.split())[-1] )
                # Check consistentcy with NAt:
                if n_Ca != NBas*NBasIndp:
                    raise ValueError( 'Number of basis functions (%i) and ' 
                            'number of MOs are not consistent (%i)'
                            'with the number of alpha MO coefficients (%i)' %
                            ( NBas, NBasIndp, n_Ca ) )
                Ca = np.zeros( n_Ca, dtype='d' )
                continue
            elif n_Ca > 0:
                for x in line.split():
                    Ca[i_Ca] = float( x )
                    i_Ca += 1
                # Check if completed:
                if i_Ca == n_Ca:
                    # Note: use order='F' to reshape column-wisely 
                    Ca = Ca.reshape( NBas, NBasIndp, order='F' )
                    n_Ca = 0
                    debug( 'Ca = ', end='' )
                    debug( Ca )
                continue

            # Beta MO coefficients
            if n_Cb == 0 and pat_Cb.match( line ):
                n_Cb = int( (line.split())[-1] )
                # Check consistentcy with NAt:
                if n_Cb != NBas*NBasIndp:
                    raise ValueError( 'Number of basis functions (%i) and ' 
                            'number of MOs (%i) are not consistent '
                            'with the number of beta MO coefficients (%i)' %
                            ( NBas, NBasIndp, n_Cb ) )
                Cb = np.zeros( n_Cb, dtype='d' )
                continue
            elif n_Cb > 0:
                for x in line.split():
                    Cb[i_Cb] = float( x )
                    i_Cb += 1
                # Check if completed:
                if i_Cb == n_Cb:
                    # Note: use order='F' to reshape column-wisely 
                    Cb = Cb.reshape( NBas, NBasIndp, order='F' )
                    n_Cb = 0
                    debug( 'Cb = ', end='' )
                    debug( Cb )
                continue

            # Mulliken charges
            if n_QMull == 0 and pat_QMull.match( line ):
                n_QMull = int( (line.split())[-1] )
                QMull = np.zeros( n_QMull, dtype='d' )
                continue
            elif n_QMull > 0:
                for x in line.split():
                    QMull[i_QMull] = float( x )
                    i_QMull += 1
                # Check if completed:
                if i_QMull == n_QMull:
                    n_QMull = 0
                    debug( 'QMull = ', end='' )
                    debug( QMull )
                continue
    # End of with open( inputFileName ) as f:

    # Determine atomic index of AO:
    atIx = np.zeros( NBas, dtype='i' ) # Atomic index of AO
    k = 0
    for i in range( 0, len(shTyp) ):
        # Number of AOs in each shell:
        el = shTyp[i];
        if el == 0 or el == 1 or el <= -2: # Pure AOs
            n_AOInSh = 2*abs(el)+1
        elif el == -1: # sp
            n_AOInSh = 4 # 1(s) + 3(p)
        else:
            el_ab = abs(el)
            n_AOInSh = (el_ab+1)*(el_ab+2)//2 # Combination number: C_{n+2}^{2}

        for j in range( 0, n_AOInSh ):
            atIx[k] = shAtMap[i]
            k += 1
    # Check consistency:
    if k != NBas:
        raise ValueError( 'Number of basis functions (%i) is not consistent '
                'with that (%i) given by shell information' % (NBas,k) )
    # Check if number of atoms is consistent:
    if NAt != atIx[-1]:
        raise ValueError( 'Number of atoms (%i) is not consistet with'
                'atomic indices of AOs (%i)' % (NAt,  atIx[-1]) )
    debug( 'atIx = ', end='' )
    debug( atIx )

    # AO indices for each atom:
    aoIx = []
    iat = 1
    istart = 0
    ix = np.arange( 1, NBas+1 )
    for i in ix:
        if atIx[i-1] > iat:
            aoIx.append( ix[ istart : i-1 ] )
            iat += 1
            istart = i-1
    aoIx.append( ix[ istart : NBas+1 ] )
    debug( 'aoIx = ', end='' )
    debug( aoIx )

    # Merge obital energies and coefficients:
    if spin == 0:
        E = Ea
        C = Ca
    else:
        E = ( Ea, Eb )
        C = ( Ca, Cb )

    #---- Collection all NAO info: ----
    return FchkInformation( chg, mult, NE, NEa, NEb, NAt, Z, Zeff, M, xyz, spin,
          NBas, NBasIndp, atIx, aoIx, E_tot, E_SCF, E, C, QMull
          )
#==============================================================================
# enddef readFchk()
