# EzReson:
#  A general efficient package for performing chemical resonance analysis of a 
#  DFT wave function
#
#  Created on Aug 26, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#
#  Updates:
#
#    May 16, 2022:
#    -  Added the `SciForm option for Kekule-WFRT jobs for both closed- and
#    open-shell cases
#
#    May 8, 2022:
#    -  Added the ProjDiv keyword to allow one to divide the complete set of
#    Kekule structures into many parts, for each of which the WF and DM 
#    projections are calculations separately, thus requiring much less memory.
#    This is useful for a huge number of Kekule structures, which would
#    otherwise encounter out of memory problem.
#
#    May 2, 2022:
#    -  Started to perform open-shell PROJ jobs for pi-conjugated
#    monoradicals and biradicals
#    -  Showing elapsed time at the end of output of the program
#
#    Feb 11, 2022:
#    -  Started to perform open-shell WFRT analysis for pi-conjugated
#    monoradicals and biradicals
#
#    Dec, 2021:
#    - Allowed to output overlap, Hamiltonian and interaction-energy matrices
#    between Clar resonators or specified Lewis structures
#
#    Oct 13, 2021:
#    - Fixed a bug for parsing the input file: 
#        The colon symbol used for LPs in Lewis structures conflicts with that
#        used in the MATLAB-style specification of atoms or LMOs
#
#    Oct 12, 2021:
#    - Extended function parseLewis() for open-shell Lewis structures
#
#    Oct 2, 2021:
#    - Allowed to write overlap and Hamiltonian matrices to external files
#
#    Jun 29, 2021:
#    - Allowed to specify LMOs and atoms using MATLAB range style (iBeg:iEnd)
#
#    Jun 28, 2021:
#    - Able to perform PROJ and WFRT jobs for Clar resonators, with options to
#    customize projCut and ClarStyle
#
#    Jun 23, 2021:
#    - Added a new job type, ENUM, for sole enumeration of Kekule or Clar 
#      structures
#
#    Jun 21, 2021:
#    - Enabled enumeration of Clar resonators
#
#    Dec 2, 2020:
#    - Fixed the bugs in runJob_WFRT() for Huckel+Kekule calculations when a
#    *.kek needs to be generated
#

VERSION = '3.0 (May 2022)'

import re
import sys
sys.stdout.flush()
import os.path
import time
from datetime import timedelta
import numpy as np
import globals
from readFchk import *
from readNBOMat import *
from transform import *
from hmo import *
from wfrt import *
from wfrt_clar import *

# To measure elapsed time:
time_start = time.time()
def printElapsedTime():
    time_end = time.time()
    seconds = time_end - time_start
    d = seconds // 86400
    seconds -= d*86400
    h = seconds// 3600
    seconds -= h*3600
    m = seconds// 60
    s = seconds - m*60
    print( 'Elapsed time:  %4i days %2i hours %2i minutes %6.3f seconds' % 
            ( d, h, m, s ) )

# Elapsed time


# Global variables:
if __name__ == "__main__": 
    globals.initialize() 


# Parameters in the input file:
class param:
    def __init__( self, basename, job, lmos, atoms, maxnlp, projcut, projdiv, \
            writeraos, flipraos, lewis, kekule, clar, huckel, precdmrt, \
            degcridmrt, clarstyle, clarmax, atomsexcl, sciform, \
            overlap, hamilt, eint ):
        self.basename = basename
        self.job = job
        self.lmos = lmos
        self.atoms = atoms
        self.maxnlp = maxnlp
        self.projcut = projcut
        self.projdiv = projdiv
        self.writeraos = writeraos
        self.flipraos = flipraos
        self.lewis = lewis
        self.kekule = kekule
        self.clar = clar
        self.huckel = huckel
        self.precdmrt = precdmrt
        self.degcridmrt = degcridmrt
        self.clarstyle = clarstyle
        self.clarmax = clarmax
        self.atomsexcl = atomsexcl
        self.sciform = sciform
        self.overlap = overlap
        self.hamilt = hamilt
        self.eint = eint
#==============================================================================
# endclass para


#==============================================================================
# Read the input file
#
# Return:
#   An instance of class papa storing all input values
#
def readControlFile( inputFileName ):

    pat_validline = re.compile( 
            '^\s*[^#]\s*[A-z]+\s*=\s*.*', re.IGNORECASE
            )

    pat_basename = re.compile( 
            '^\s*File\s*=\s*.*', re.IGNORECASE
            )

    pat_job = re.compile( 
            '^\s*Job\s*=\s*.*', re.IGNORECASE
            )

    pat_lmos = re.compile( 
            '^\s*LMOs\s*=\s*.*', re.IGNORECASE
            )

    pat_atoms = re.compile( 
            '^\s*Atoms\s*=\s*.*', re.IGNORECASE
            )

    pat_maxnlp = re.compile( 
            '^\s*MaxNLP\s*=\s*.*', re.IGNORECASE
            )

    pat_projcut = re.compile( 
            '^\s*ProjCut\s*=\s*.*', re.IGNORECASE
            )

    pat_projdiv = re.compile( 
            '^\s*ProjDiv\s*=\s*.*', re.IGNORECASE
            )

    pat_writeraos = re.compile( 
            '^\s*WriteRAOs\s*=\s*.*', re.IGNORECASE
            )

    pat_flipraos = re.compile( 
            '^\s*FlipRAOs\s*=\s*.*', re.IGNORECASE
            )

    pat_lewis = re.compile( 
            '^\s*Lewis\s*=\s*.*', re.IGNORECASE
            )

    pat_kekule = re.compile( 
            '^\s*Kekule\s*=\s*.*', re.IGNORECASE
            )

    pat_clar = re.compile( 
            '^\s*Clar\s*=\s*.*', re.IGNORECASE
            )

    pat_huckel = re.compile( 
            '^\s*Huckel\s*=\s*.*', re.IGNORECASE
            )

    pat_precdmrt = re.compile( 
            '^\s*PrecDMRT\s*=\s*.*', re.IGNORECASE
            )

    pat_degcridmrt = re.compile( 
            '^\s*DegCriDMRT\s*=\s*.*', re.IGNORECASE
            )

    pat_clarstyle = re.compile( 
            '^\s*ClarStyle\s*=\s*.*', re.IGNORECASE
            )

    pat_clarmax = re.compile( 
            '^\s*ClarMax\s*=\s*.*', re.IGNORECASE
            )

    pat_sciform = re.compile( 
            '^\s*SciForm\s*=\s*.*', re.IGNORECASE
            )

    pat_atomsexcl = re.compile( 
            '^\s*AtomsExcl\s*=\s*.*', re.IGNORECASE
            )

    pat_scered = re.compile( 
            '^\s*SCERED\s*=\s*.*', re.IGNORECASE
            )

    pat_overlap = re.compile( 
            '^\s*Overlap\s*=\s*.*', re.IGNORECASE
            )

    pat_hamilt = re.compile( 
            '^\s*Hamiltonian\s*=\s*.*', re.IGNORECASE
            )

    pat_eint = re.compile( 
            '^\s*Eint\s*=\s*.*', re.IGNORECASE
            )

    # Default values:
    basename = ''
    job = ''
    lmos = []
    atoms = []
    maxnlp = inf
    projcut = 0.
    projdiv = 1
    writeraos = False
    flipraos = None
    lewis = []
    kekule = False
    clar = False
    huckel = False
    precdmrt = 0.99
    degcridmrt = 1E-3
    clarstyle = 'short'
    clarmax = False
    atomsexcl = []
    sciform = False
    overlap = False
    hamilt = False
    eint = False

    iline = 0
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            if pat_validline.match( line ):
                fields = re.sub( '^.*=', '', line )
                if 'lewis' not in line.lower():
                    fields = re.sub( '\s*:\s*', ':', fields )
                if len( fields.strip() ) == 0:
                    print( 'ERROR: Value missing at line %i:' % iline )
                    print( '    -->%s' % line, end='' )
                    print( 'Aborted' )
                    exit(1)

            # Basename of DFT output files:
            if pat_basename.match( line ):
                basename = fields.split()[0]
                # print( 'basename = %s' % basename )
                continue

            # Job type:
            if pat_job.match( line ):
                job = fields.split()[0].upper()
                #print( 'job = %s' % job )
                continue

            # Indices of LMOs:
            if pat_lmos.match( line ):
                if fields.strip().lower() == 'pi':
                    lmos = -1
                    continue
                for s in fields.split():
                    s_BegEnd = s.split( ':' )
                    if len( s_BegEnd ) == 1:
                        lmos.append( int(s) )
                    elif len( s_BegEnd ) == 2:
                        for i in range( int(s_BegEnd[0]), int(s_BegEnd[1])+1 ):
                            lmos.append( i )
               #print( 'lmos = ', lmos )
                #input()
                continue

            # Indices of atoms:
            if pat_atoms.match( line ):
                for s in fields.split():
                    s_BegEnd = s.split( ':' )
                    if len( s_BegEnd ) == 1:
                        atoms.append( int(s) )
                    elif len( s_BegEnd ) == 2:
                        for i in range( int(s_BegEnd[0]), int(s_BegEnd[1])+1 ):
                            atoms.append( i )
                #print( 'atoms = ', atoms )
                #input()
                continue

            # Maximum number of lone pairs:
            if pat_maxnlp.match( line ):
                s = fields.split()[0]
                if s.lower() == 'inf':
                    maxnlp = inf
                else:
                    maxnlp = int( s )
                continue

            # Cutoff value of projections relative to major contributor:
            if pat_projcut.match( line ):
                s = fields.split()[0]
                projcut = float( s )
                continue

            # Number of division of resonators for a projection job:
            if pat_projdiv.match( line ):
                s = fields.split()[0]
                projdiv = int( s )
                continue

            # Whether write RAOs to an external *.fchk file or not:
            if pat_writeraos.match( line ):
                s = fields.split()[0]
                if s.upper() == 'TRUE' or s == '1':
                    writeraos = True
                elif s.upper() == 'FALSE' or s == '0':
                    writeraos = False
                else:
                    print( 'ERROR: Invalid value of parameter WriteRAOs:', s )
                    print( '  Valid options are TRUE, FALSE, 1, 0' )
                    print( 'Aborted' )
                    exit(1)
                continue

            # Indices of atoms whose RAO phase is to be flipped:
            if pat_flipraos.match( line ):
                flipraos = []
                s = fields.split()
                if len(s) == 1:
                    n = int( s[0] )
                    if n == -1:
                        flipraos = None
                    elif n >= 0:
                        flipraos = [ n ]
                    else:
                        print( 'ERROR: Invalid value of parameter FlipRAOs:', \
                                s )
                        print( '  Valid values are -1, 0 and positive numbers' )
                        print( 'Aborted' )
                        exit(1)
                else:
                    for si in s:
                        flipraos.append( int(si) )
                continue

            # Indices of the specified Lewis structures:
            if pat_lewis.match( line ):
                LP = []
                BD = []
                RA = []
                RB = []
                ifOpenShell = False
                for s in fields.split():
                    # If s is an index of Lewis structures:
                    if ':' not in s and '/' not in s:
                        lewis.append( int(s) )
                    else:
                        lp, bd, ra, rb = parseLewis( s )
                        LP.append( lp )
                        BD.append( bd )
                        RA.append( ra )
                        RB.append( rb )
                        if len(ra) > 0 or len(rb) > 0:
                            ifOpenShell = True
                # Check validity of the input Lewis structures:
                # Num. of alpha electrons:
                NEa0 = len(LP[0]) + len(BD[0]) + len( RA[0] ) 
                # Num. of beta electrons:
                NEb0 = len(LP[0]) + len(BD[0]) + len( RB[0] ) 
                for i in range( 1, len(LP) ):
                    NEai = len(LP[i]) + len(BD[i]) + len( RA[i] ) 
                    if NEai  != NEa0:
                        print( 'ERROR: Lewis structure #%i has a different '
                                'number of alpha electrons' % (i+1) )
                        exit(1)
                    NEbi = len(LP[i]) + len(BD[i]) + len( RB[i] ) 
                    if NEbi  != NEb0:
                        print( 'ERROR: Lewis structure #%i has a different '
                                'number of beta electrons' % (i+1) )
                        exit(1)

                # Closed-shell case:
                if not ifOpenShell:
                    print( 'Closed-shell Lewis structures read in' )
                    lewis = ( LP, BD )
                else:
                    lewis = ( LP, BD, RA, RB )
                    print( 'Open-shell Lewis structures read in' )
                #print( 'lewis = ', lewis )
                #input()
                continue

            # Whether using Kekule structures only:
            if pat_kekule.match( line ):
                s = fields.split()[0]
                if s.upper() == 'TRUE' or s == '1':
                    kekule = True
                elif s.upper() == 'FALSE' or s == '0':
                    kekule = False
                else:
                    print( 'ERROR: Invalid value of parameter Kekule:', s )
                    print( '  Valid options are TRUE, FALSE, 1, 0' )
                    print( 'Aborted' )
                    exit(1)
                continue

            # Whether using Clar resonators only:
            if pat_clar.match( line ):
                s = fields.split()[0]
                if s.upper() == 'TRUE' or s == '1':
                    clar = True
                elif s.upper() == 'FALSE' or s == '0':
                    clar = False
                else:
                    print( 'ERROR: Invalid value of parameter Clar:', s )
                    print( '  Valid options are TRUE, FALSE, 1, 0' )
                    print( 'Aborted' )
                    exit(1)
                continue

            # Whether in the framework of Huckel molecular orbitals theory:
            if pat_huckel.match( line ):
                s = fields.split()[0]
                if s.upper() == 'TRUE' or s == '1':
                    huckel = True
                elif s.upper() == 'FALSE' or s == '0':
                    huckel = False
                else:
                    print( 'ERROR: Invalid value of parameter Huckel:', s )
                    print( '  Valid options are TRUE, FALSE, 1, 0' )
                    print( 'Aborted' )
                    exit(1)
                continue

            # Minimum reproducibility that must be achieved:
            if pat_precdmrt.match( line ):
                s = fields.split()[0]
                precdmrt = float( s )
                continue

            # Criterium to determine degenerate Lewis structures:
            if pat_degcridmrt.match( line ):
                s = fields.split()[0]
                degcridmrt = float( s )
                continue

            # Style of printing Clar resonators:
            if pat_clarstyle.match( line ):
                clarstyle = fields.split()[0].lower()
                if clarstyle != 'short' and clarstyle != 'explicit':
                    print( 'ERROR: Invalid value of parameter ClarStyle:',
                            clarstyle )
                    print( '  Valid options are Short and Explicit' )
                    print( 'Aborted' )
                    exit(1)
                continue

            # Whether writing the Clar structures:
            if pat_clarmax.match( line ):
                s = fields.split()[0]
                if s.upper() == 'TRUE' or s == '1':
                    clarmax = True
                elif s.upper() == 'FALSE' or s == '0':
                    clarmax = False
                else:
                    print( 'ERROR: Invalid value of parameter ClarMax:', s )
                    print( '  Valid options are TRUE, FALSE, 1, 0' )
                    print( 'Aborted' )
                    exit(1)
                continue

            # Indices of excluded atoms:
            if pat_atomsexcl.match( line ):
                for s in fields.split():
                    s_BegEnd = s.split( ':' )
                    if len( s_BegEnd ) == 1:
                        atomsexcl.append( int(s) )
                    elif len( s_BegEnd ) == 2:
                        for i in range( int(s_BegEnd[0]), int(s_BegEnd[1])+1 ):
                            atomsexcl.append( i )
                # Sort in ascending order:
                atomsexcl.sort()
                #print( 'atomsexcl = ', atomsexcl )
                #input()
                continue

            # Whether showing the numerical results in scientific format:
            if pat_sciform.match( line ):
                s = fields.split()[0]
                if s.upper() == 'TRUE' or s == '1':
                    sciform = True
                elif s.upper() == 'FALSE' or s == '0':
                    sciform = False
                else:
                    print( 'ERROR: Invalid value of parameter SciForm:', s )
                    print( '  Valid options are TRUE, FALSE, 1, 0' )
                    print( 'Aborted' )
                    exit(1)
                continue

            # Whether writing the overlap matrix of resonance structures to an
            # external file:
            if pat_overlap.match( line ):
                s = fields.split()[0]
                if s.upper() == 'TRUE' or s == '1':
                    overlap = True
                elif s.upper() == 'FALSE' or s == '0':
                    overlap = False
                else:
                    print( 'ERROR: Invalid value of parameter Overlap:', s )
                    print( '  Valid options are TRUE, FALSE, 1, 0' )
                    print( 'Aborted' )
                    exit(1)
                continue

            # Whether writing the Hamiltonian matrix of resonance structures 
            # to an external file:
            if pat_hamilt.match( line ):
                s = fields.split()[0]
                if s.upper() == 'TRUE' or s == '1':
                    hamilt = True
                elif s.upper() == 'FALSE' or s == '0':
                    hamilt = False
                else:
                    print( 'ERROR: Invalid value of parameter Hamilt:', s )
                    print( '  Valid options are TRUE, FALSE, 1, 0' )
                    print( 'Aborted' )
                    exit(1)
                continue

            # Whether writing the matrix of interaction energy betwenen the 
            # resonance structures to an external file:
            if pat_eint.match( line ):
                s = fields.split()[0]
                if s.upper() == 'TRUE' or s == '1':
                    eint = True
                elif s.upper() == 'FALSE' or s == '0':
                    eint = False
                else:
                    print( 'ERROR: Invalid value of parameter Eint:', s )
                    print( '  Valid options are TRUE, FALSE, 1, 0' )
                    print( 'Aborted' )
                    exit(1)
                continue

    # End of with open( inputFileName ) as f

    # Check validity of input values:
    # -- No basename specified:
    if not basename:
        print( 'ERROR: ', end='' )
        print( 'Missing \"File = ...\" to specify the basename of the DFT '
               'output files' )
        print( 'Aborted' )
        exit(1)

    # -- No job specified:
    if not job:
        print( 'ERROR: ', end='' )
        print( 'Missing \"Job = ...\" to specify the job type' )
        #print( '  Valid options are ENUM, LMO, WFRT, DMRT, '
        #        'EMRT, PROJ, SCERED' )
        print( '  Valid options are ENUM, LMO, WFRT, DMRT, '
                'PROJ, SCERED' )
        print( 'Aborted' )
        exit(1)

    # -- For Kekule/Clar resonators, # of electrons must equal # of atoms:
    if kekule or clar:
        if isinstance(lmos, list) and ( len( lmos )*2 != len( atoms ) ) \
                and not huckel:
            print( 'ERROR: ', end='' )
            print( 'For kekule/clar structures, the number of LMOs (%i) '
                    'is not half of the number of atoms (%i)' % 
                    ( len( lmos ), len( atoms ) ) )
            print( 'Aborted' )
            exit(1)

    # -- Replusive parameters:
    if len( lewis ) > 0 and kekule:
        print( 'ERROR: ', end='' )
        print( 'Lewis structures cannot be specified while Kekule option is ' \
                'turned on' )
        print( 'Aborted' )
        exit(1)
    if len( lewis ) > 0 and clar:
        print( 'ERROR: ', end='' )
        print( 'Lewis structures cannot be specified while Clar option is ' \
                'turned on' )
        print( 'Aborted' )
        exit(1)
    if kekule and clar:
        print( 'ERROR: ', end='' )
        print( 'The Kekule and the Clar options cannot be turned on at the ' \
                'same time' )
        print( 'Aborted' )
        exit(1)
    if kekule and clarmax:
        print( 'ERROR: ', end='' )
        print( 'The Kekule and the ClarMax options cannot be turned on ' \
                'at the same time' )
        print( 'Aborted' )
        exit(1)

    # -- The atomsexcl option is only for Clar resonators
    if (not clar) and (not clarmax) and len( atomsexcl ) > 0:
        print( 'ERROR: ', end='' )
        print( 'The AtomsExcl option is only valid '
                'when the Clar option is on' )
        print( 'Aborted' )
        exit(1)

    # -- Currently only support closed-shell Clar resonators:
    if len( atomsexcl ) > 0 and len( atomsexcl ) % 2 == 1:
        print( 'ERROR: ', end='' )
        print( 'The number of excluded atoms must be even' )
        print( 'Aborted' )
        exit(1)

    # -- The Clar and ClarMax options are mutually exclusive except for the
    # ENUM job
    if clar and clarmax and job != 'ENUM':
        print( 'ERROR: ', end='' )
        print( 'The Clar and ClarMax options are mutually exclusive '
                'except for the ENUM job' )
        print( 'Aborted' )
        exit(1)

    return param( basename, job, lmos, atoms, maxnlp, projcut, projdiv, \
            writeraos, flipraos, lewis, kekule, clar, huckel, precdmrt, \
            degcridmrt, clarstyle, clarmax, atomsexcl, sciform, overlap, \
            hamilt, eint )
#==============================================================================
# enddef readControlFile()


#==============================================================================
# Parse a string of Lewis structure
#
def parseLewis( s ):
    lp = [] # Lone pair list
    bd = [] # bond list
    ra = [] # alpha-radical site list
    rb = [] # beta-radical site list
    f = s.split( '/' )
    for ep in f: # for each electron pair:
        if ':' in ep: # A lone pair:
            n = ep.replace( ':', '' )
            if len(n) != 1:
                raise ValueError( 'Invalid Lewis structure:', s )
            if int(n) < 1:
                raise ValueError( 'Invalid Lewis structure:', s )

            lp.append( ep.replace( ':', '' ) )
        elif '-' in ep: # A covalent bond
            b = ep.split( '-' )
            if len(b) != 2:
                raise ValueError( 'Invalid Lewis structure:', s )
            b1 = int( b[0] )
            b2 = int( b[1] )
            if b1 < 1 or b2 < 1:
                raise ValueError( 'Invalid Lewis structure:', s )
            bd.append( b1 )
            bd.append( b2 )
        elif '\'' in ep or '.' in ep: # A biradical
            # Alpha sites:
            if '\'' in ep and '.' in ep:
                print( 'Error: Alpha and beta sites must be specified '
                        'separately' )
                raise ValueError( 'Invalid Lewis structure:', s )
            if '\'' in ep:
                r_a = ep.split( '\'' )
                for at in r_a:
                    if len(at) == 0:
                        continue
                    atIx = int( at )
                    if atIx >= 1:
                        ra.append( atIx )
            # Beta sites:
            if '.' in ep:
                r_b = ep.split( '.' )
                for at in r_b:
                    if len(at) == 0:
                        continue
                    atIx = int( at )
                    if atIx >= 1:
                        rb.append( atIx )
        else:
            raise ValueError( 'Invalid Lewis structure:', s )

    lp = np.array( lp, dtype='i' )
    nbd = len(bd) // 2
    bd = np.array( bd, dtype='i' ).reshape( nbd, 2 )
    ra = np.array( ra, dtype='i' )
    rb = np.array( rb, dtype='i' )

    return ( lp, bd, ra, rb )
#==============================================================================
# enddef parseLewis()

#==============================================================================
# Perform an ENUM job
#
def runJob_ENUM( basename, kekule, clar, clarmax, atomsexcl ):
    elem, xyz = readxyz( basename + '.xyz' )
    # Remove hydrogens:
    xyz, _ = xyzNonH( elem, xyz )
    if kekule:
        print( 'Performing combinatorial enumeration of '
                'Kekule structures ...' )
        # Generate all possible Kekule structures:
        kekFileName = basename + '.kek'
        if os.path.exists( kekFileName ):
            print( 'Kekule structures have been previously enumerated in '
                    'file', kekFileName )
        else:
            with open( kekFileName, 'w' ) as writer:
                print( 'Enumerating Kekule structures for %s ...' % basename )
                enumKekule( writer, elem, xyz )
                print( 'File', kekFileName, 'written' )

    if clar:
        print( 'Performing combinatorial enumeration of '
                'Clar resonators ...' )
        # Generate all possible Clar resonators:
        clarFileName = basename + '.clar'
        # In the case of excluding some atoms:
        if len( atomsexcl ) > 0:
            str_atomsexcl = '_pos'
            for ix in atomsexcl:
                str_atomsexcl += '_%i' % ix
            clarFileName = basename + str_atomsexcl + '.clar'
        if os.path.exists( clarFileName ):
            print( 'Clar resonators have been previously enumerated in '
                    'file', clarFileName )
        else:
            with open( clarFileName, 'w' ) as writer:
                print( 'Enumerating Clar resonators for %s ...' % basename )
                if len( atomsexcl ) > 0:
                    print( '--> Excluding atoms:', end='' )
                    for ix in atomsexcl:
                        print( ' %i' % ix, end='' )
                    print()
                enumClar( writer, elem, xyz, False, atomsexcl, clarFileName )
                #enumClar_from_Kekule(re.sub('\.clar$', '.kek', clarFileName),\
                #        writer, elem, xyz, atomsexcl, clarFileName )
                print( 'File', clarFileName, 'written' )

    if clarmax:
        clarFileName = basename + '.clar'
        # In the case of excluding some atoms:
        if len( atomsexcl ) > 0:
            str_atomsexcl = '_pos'
            for ix in atomsexcl:
                str_atomsexcl += '_%i' % ix
            clarFileName = basename + str_atomsexcl + '.clar'

        # 1. If Clar resonators have already been generated
        if os.path.exists( clarFileName ):
            print( 'Extracting Clar structures from file %s ...' % 
                    clarFileName )
            # Generate *.clarmax file that includes all Clar structures
            # (i.e., the Clar resonantors with max. num. of sextets)
            clarmaxFileName = clarFileName + 'max'
            if os.path.exists( clarmaxFileName ):
                print( 'Clar structures have been previously enumerated in '
                        'file', clarmaxFileName )
            else:
                clarMax_from_clar( clarmaxFileName, atomsexcl )
                print( 'File', clarmaxFileName, 'written' )

        # 2. Request for direct generation of Clar str. without Clar resonators
        else:
            print( 'Performing combinatorial enumeration of '
                    'Clar structures ...' )
            # Generate all possible Clar structures:
            clarmaxFileName = basename + '.clarmax'
            # In the case of excluding some atoms:
            if len( atomsexcl ) > 0:
                str_atomsexcl = '_pos'
                for ix in atomsexcl:
                    str_atomsexcl += '_%i' % ix
                clarmaxFileName = basename + str_atomsexcl + '.clarmax'
            if os.path.exists( clarmaxFileName ):
                print( 'Clar structures have been previously enumerated in '
                        'file', clarmaxFileName )
            else:
                with open( clarmaxFileName, 'w' ) as writer:
                    print( 'Enumerating Clar structures for %s ...' % 
                            basename )
                    if len( atomsexcl ) > 0:
                        print( '--> Excluding atoms:', end='' )
                        for ix in atomsexcl:
                            print( ' %i' % ix, end='' )
                        print()
                    enumClar( writer, elem, xyz, True, atomsexcl, \
                            clarmaxFileName )
                    print( 'File', clarmaxFileName, 'written' )

#==============================================================================
# enddef runJob_ENUM()

#==============================================================================
# Perform an LMO job
#
def runJob_LMO( basename ):
    print( 'Performing Pipek-Mezey localization of molecular orbitals ...' )
    FchkInfo = readFchk( basename + '.fchk' )
    # Check if AOs has no linear dependence (so that C is invertible):
    if FchkInfo.NBasIndp < FchkInfo.NBas:
        raise ValueError( 'Linear dependence found in the AOs\nAborted!' )
    naoInfo = readNAOInfo( basename +'.out' )
    CAONAO = readNBOMat( basename + '.33', naoInfo.NBas, naoInfo.NNAO  )
    CNAOLMO, ELMO = genLMOs( FchkInfo, naoInfo, CAONAO, basename )
#==============================================================================
# enddef runJob_LMO()


#==============================================================================
# Perform a WFRT job
#
def runJob_WFRT( basename, lmos, atoms, maxnlp, projcut, writeraos, \
        flipraos, lewis, kekule, clar, clarmax, huckel, clarstyle, atomsexcl, \
        sciform ):
    print( 'Performing Wave Function based Resonance Theory (WFRT) '
           'calculations ...' )
    if not huckel:
        FchkInfo = readFchk( basename + '.fchk' )
        # Check if AOs has no linear dependence (so that C is invertible):
        if FchkInfo.NBasIndp < FchkInfo.NBas:
            raise ValueError( 'Linear dependence found in the AOs\nAborted!' )
        naoInfo = readNAOInfo( basename +'.out' )
        CAONAO = readNBOMat( basename + '.33', naoInfo.NBas, naoInfo.NNAO  )
        CNAOLMO, ELMO = genLMOs( FchkInfo, naoInfo, CAONAO, basename )

        # In the case of pi LMOs:
        if isinstance( lmos, int ):
            # For closed-shell system:
            if FchkInfo.spin == 0:
                nHOMO = FchkInfo.NEa
                NLMO = len( atoms ) // 2
                # Indices of pi-LMOs. NOTE: starting from 1
                lmos = range( nHOMO - NLMO + 1, nHOMO + 1 )   
                print( 'The pi resonance system is automatically defined '
                        'by LMOs:' )
                for ix in lmos:
                    print( ' %i' % ix, end='' )
            else:
                NEa_b = FchkInfo.NEa - FchkInfo.NEb
                NLMOa = ( len( atoms ) + NEa_b ) // 2
                NLMOb = ( len( atoms ) - NEa_b ) // 2
                if ( NLMOa + NLMOb ) != len( atoms ) or \
                        ( NLMOa - NLMOb ) != NEa_b:
                    raise ValueError( 'Number of pi electrons is not '
                            'consistent with spin multiplicity\n'
                            'Check the list of Atoms and make sure the '
                            'nnumber of atoms is correct.' )
                nHOMOa = FchkInfo.NEa
                nHOMOb = FchkInfo.NEb
                # Indices of pi-LMOs. NOTE: starting from 1
                lmos_a = range( nHOMOa - NLMOa + 1, nHOMOa + 1 )   
                lmos_b = range( nHOMOb - NLMOb + 1, nHOMOb + 1 )   
                lmos = ( lmos_a, lmos_b )
                print( 'The pi resonance system is automatically defined '
                        'by the alpha and beta LMOs:' )
                print( 'Alpha LMOs:' )
                for ix in lmos_a:
                    print( ' %i' % ix, end='' )
                print()
                print( 'Beta LMOs:' )
                for ix in lmos_b:
                    print( ' %i' % ix, end='' )
            print()
    else:
        coordFile = basename + '.gjf'
        if not os.path.exists( coordFile ):
            coordFile = basename + '.xyz' 
        hmoSol = hmo( coordFile )

    # Wave-function based resonance theory (WFRT) analysis:
    if projcut != 0:
        print( 'Using Lewis structures with a projection cutoff of %.6E' % \
                projcut )
    if writeraos:
        rao_tag = basename
        print( 'RAOs are to be written in the file %s_RAO.fchk' % rao_tag )
    else:
        rao_tag = ''
    if kekule:
        print( 'Using only Kekule structures to expand the wave function' )
    if clar:
        print( 'Using only Clar resonators to expand the wave function' )
    if clarmax:
        print( 'Using only Clar structures to expand the wave function' )

    if kekule:
        # Generate all possible Kekule structures:
        kekFileName = basename + '.kek'
        if os.path.exists( kekFileName ):
            print( 'Kekule structures have been previously enumerated in '
                    'file', kekFileName )
        else:
            elem = []
            if huckel:
                elem = 'C' * hmoSol.NAt
            else:
                for Z in FchkInfo.Z:
                    elem.append( str(Z) )
            with open( kekFileName, 'w' ) as writer:
                print( 'Enumerating Kekule structures for %s ...' % basename )
                if huckel:
                    enumKekule( writer, elem, hmoSol.xyz )
                else:
                    enumKekule( writer, elem, FchkInfo.xyz, FchkInfo.spin, \
                            FchkInfo.mult )
                print( 'File', kekFileName, 'written' )

    if clar:
        # Generate all possible Clar structures:
        clarFileName = basename + '.clar'
        if os.path.exists( clarFileName ):
            print( 'Clar resonators have been previously enumerated in '
                    'file', clarFileName )
        else:
            elem = []
            if huckel:
                elem = 'C' * hmoSol.NAt
            else:
                for Z in FchkInfo.Z:
                    elem.append( str(Z) )
            with open( clarFileName, 'w' ) as writer:
                print( 'Enumerating Clar resonators for %s ...' % basename )
                if huckel:
                    enumClar( writer, elem, hmoSol.xyz, False, atomsexcl )
                else:
                    enumClar( writer, elem, FchkInfo.xyz, False, atomsexcl )
                print( 'File', clarFileName, 'written' )
    if clarmax:
        # Generate all possible Clar structures:
        clarmaxFileName = basename + '.clarmax'
        # In the case of excluding some atoms:
        if len( atomsexcl ) > 0:
            str_atomsexcl = '_pos'
            for ix in atomsexcl:
                str_atomsexcl += '_%i' % ix
            clarmaxFileName = basename + str_atomsexcl + '.clarmax'
        clarFileName = re.sub( '\.clarmax$', '.clar', clarmaxFileName )
        if os.path.exists( clarmaxFileName ):
            print( 'Clar structures have been previously enumerated in '
                    'file', clarmaxFileName )
        else:
            elem = []
            if huckel:
                elem = 'C' * hmoSol.NAt
            else:
                for Z in FchkInfo.Z:
                    elem.append( str(Z) )
            with open( clarmaxFileName, 'w' ) as writer:
                # 1. If *.clar file exists:
                if os.path.exists( clarFileName ):
                    print( 'Extracting Clar structures from file %s ...' % 
                            clarFileName )
                    clarMax_from_clar( clarmaxFileName, atomsexcl )
                    print( 'File', clarmaxFileName, 'written' )
                # 2. Request for direct generation of Clar structures:
                else:
                    print( 'Performing combinatorial enumeration of '
                            'Clar structures ...' )
                    with open( clarmaxFileName, 'w' ) as writer:
                        print( 'Enumerating Clar structures for %s ...' % 
                                basename )
                        if len( atomsexcl ) > 0:
                            print( '--> Excluding atoms:', end='' )
                            for ix in atomsexcl:
                                print( ' %i' % ix, end='' )
                            print()
                        if huckel:
                            enumClar( writer, elem, hmoSol.xyz, True, \
                                    atomsexcl, clarmaxFileName )
                        else:
                            enumClar( writer, elem, FchkInfo.xyz, True, \
                                    atomsexcl, clarmaxFileName )
                        print( 'File', clarmaxFileName, 'written' )

    print()

    # wfrt_hmo:
    if huckel:
        print( 'In the framework of simple Huckel molecular orbitals theory' )
        if len( lewis ) == 0 and not kekule and not clar and not clarmax:
            wfrt_hmo( hmoSol, atoms, maxnlp, projcut )
        elif kekule:
            wfrt_hmo_kekule( hmoSol, kekFileName, sciform )
        elif clar:
            wfrt_clar_hmo( hmoSol, clarFileName, 1, projcut, clarstyle )
        elif clarmax:
            wfrt_clar_hmo( hmoSol, clarmaxFileName, 1, projcut, clarstyle )
        else:
            wfrt_hmo_spec( hmoSol, lewis )
        exit(0)

    # wfrt:
    if len( lewis ) == 0 and not kekule and not clar and not clarmax:
        wfrt( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, maxnlp, \
            projcut, 'uni', rao_tag, flipraos )
    # wfrt_kekule:
    elif kekule:
        wfrt_kekule( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
                 kekFileName, sciform, 'uni', rao_tag )
    elif clar:
        wfrt_clar( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
                 clarFileName, 1, projcut, clarstyle, sciform, 'uni' )
    elif clarmax:
        wfrt_clar( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
                 clarmaxFileName, 1, projcut, clarstyle, sciform, 'uni' )
    # wfrt_spec:
    else:
        wfrt_spec( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
            lewis, 'uni', rao_tag, flipraos )
#==============================================================================
# enddef runJob_WFRT()



#==============================================================================
# Perform a EMRT job (ONLY in the Huckel framework)
#
def runJob_EMRT( basename, lmos, atoms, maxnlp, projcut, writeraos, \
        flipraos, lewis, kekule, clar, clarmax, huckel, clarstyle, atomsexcl, \
        sciform ):
    print( 'Performing Energy Minimization Resonance Theory (EMRT) '
           'calculations with the Huckel approximation...' )

    if not huckel:
        raise ValueError( 'EMRT calculation is only supported within the '
                'simple Huckel framework' )

    coordFile = basename + '.gjf'
    if not os.path.exists( coordFile ):
        coordFile = basename + '.xyz' 
    hmoSol = hmo( coordFile )

    if len( lewis ) == 0 and not kekule and not clar and not clarmax:
        pass
    elif kekule:
        pass
    elif clar:
        pass
    elif clarmax:
        pass
    else:
        emrt_spec( hmoSol, lewis )
    exit(0)


#==============================================================================
# Perform a DMRT job
def runJob_DMRT( basename, lmos, atoms, maxnlp, precdmrt, degcridmrt, \
        lewis, kekule, huckel ):
    print( 'Performing Density Matrix based Resonance Theory (DMRT) '
           'calculations ...' )
    if not huckel:
        FchkInfo = readFchk( basename + '.fchk' )
        # Check if AOs has no linear dependence (so that C is invertible):
        if FchkInfo.NBasIndp < FchkInfo.NBas:
            raise ValueError( 'Linear dependence found in the AOs\nAborted!' )
        naoInfo = readNAOInfo( basename +'.out' )
        CAONAO = readNBOMat( basename + '.33', naoInfo.NBas, naoInfo.NNAO  )
        CNAOLMO, ELMO = genLMOs( FchkInfo, naoInfo, CAONAO, basename )

        # In the case of pi LMOs:
        if isinstance( lmos, int ):
            nHOMO = FchkInfo.NEa
            NLMO = len( atoms ) // 2
            # Indices of pi-LMOs. NOTE: starting from 1
            lmos = range( nHOMO - NLMO + 1, nHOMO + 1 )   
            print( 'The pi resonance system is automatically defined '
                    'by LMOs:' )
            for ix in lmos:
                print( ' %i' % ix, end='' )
            print()
    else:
        hmoSol = hmo( basename + '.xyz' )
        #hmoSol = hmo( basename + '.gjf' )

    # Density-matrix based resonance theory (DMRT) analysis:
    print( 'Minimum reproducibility of %.6f%% is required' % \
            (precdmrt*100) )
    print( 'Using a criterium of %.6E to determine degenerate Lewis '
            'structures' % degcridmrt )
    if kekule:
        print( 'Using only Kekule structures to expand the wave function' )

    print()

    # dmrt_hmo:
    if huckel:
        print( 'In the framework of simple Huckel molecular orbitals theory' )
        # if len( lewis ) == 0 and not kekule:  # Not supported yet
        if len( lewis ) == 0:
            dmrt_hmo( hmoSol, maxnlp, 1., [], [], precdmrt, degcridmrt )
        else:
            dmrt_hmo_spec( hmoSol, lewis )
        exit(0)

    # dmrt:
    if len( lewis ) == 0 and not kekule:
        dmrt( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, maxnlp, \
            1., [], [], precdmrt, degcridmrt )
    # dmrt_kekule:
    elif kekule:
        # Generate all possible Kekule structures:
        kekFileName = basename + '.kek'
        if os.path.exists( kekFileName ):
            print( 'Kekule structures have been previously enumerated in '
                    'file', kekFileName )
        else:
            elem = []
            for Z in FchkInfo.Z:
                elem.append( str(Z) )
            with open( kekFileName, 'w' ) as writer:
                print( 'Enumerating Kekule structures for %s ...' % basename )
                enumKekule( writer, elem, FchkInfo.xyz )
                print( 'File', kekFileName, 'written' )

        dmrt_kekule( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
                 kekFileName, 'uni' )
    # dmrt_spec:
    else:
        dmrt_spec( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
            lewis, 'uni' )
#==============================================================================
# enddef runJob_DMRT()



#==============================================================================
# Perform a PROJ job
#
def runJob_PROJ( basename, lmos, atoms, maxnlp, writeraos, flipraos, \
        lewis, kekule, clar, clarmax, projdiv, huckel, clarstyle, atomsexcl, \
        sciform ):
    print( 'Performing wave function and density matrix projection '
           'calculations ...' )
    if not huckel:
        FchkInfo = readFchk( basename + '.fchk' )
        # Check if AOs has no linear dependence (so that C is invertible):
        if FchkInfo.NBasIndp < FchkInfo.NBas:
            raise ValueError( 'Linear dependence found in the AOs\nAborted!' )
        naoInfo = readNAOInfo( basename +'.out' )
        CAONAO = readNBOMat( basename + '.33', naoInfo.NBas, naoInfo.NNAO  )
        CNAOLMO, ELMO = genLMOs( FchkInfo, naoInfo, CAONAO, basename )

        # In the case of pi LMOs:
        if isinstance( lmos, int ):
            # For closed-shell system:
            if FchkInfo.spin == 0:
                nHOMO = FchkInfo.NEa
                NLMO = len( atoms ) // 2
                # Indices of pi-LMOs. NOTE: starting from 1
                lmos = range( nHOMO - NLMO + 1, nHOMO + 1 )
                print( 'The pi resonance system is automatically defined '
                        'by LMOs:' )
                for ix in lmos:
                    print( ' %i' % ix, end='' )
            else:
                NEa_b = FchkInfo.NEa - FchkInfo.NEb
                NLMOa = ( len( atoms ) + NEa_b ) // 2
                NLMOb = ( len( atoms ) - NEa_b ) // 2
                if ( NLMOa + NLMOb ) != len( atoms ) or \
                        ( NLMOa - NLMOb ) != NEa_b:
                    raise ValueError( 'Number of pi electrons is not '
                            'consistent with spin multiplicity\n'
                            'Check the list of Atoms and make sure the '
                            'nnumber of atoms is correct.' )
                nHOMOa = FchkInfo.NEa
                nHOMOb = FchkInfo.NEb
                # Indices of pi-LMOs. NOTE: starting from 1
                lmos_a = range( nHOMOa - NLMOa + 1, nHOMOa + 1 )
                lmos_b = range( nHOMOb - NLMOb + 1, nHOMOb + 1 )
                lmos = ( lmos_a, lmos_b )
                print( 'The pi resonance system is automatically defined '
                        'by the alpha and beta LMOs:' )
                print( 'Alpha LMOs:' )
                for ix in lmos_a:
                    print( ' %i' % ix, end='' )
                print()
                print( 'Beta LMOs:' )
                for ix in lmos_b:
                    print( ' %i' % ix, end='' )
            print()
    else:
        coordFile = basename + '.gjf'
        if not os.path.exists( coordFile ):
            coordFile = basename + '.xyz'
        hmoSol = hmo( coordFile )

    # Calculating WF and DM projections:
    if writeraos:
        rao_tag = basename
        print( 'RAOs are to be written in the file %s_RAO.fchk' % rao_tag )
    else:
        rao_tag = ''
    if kekule:
        print( 'Computing projections only for Kekule structures' )
    if clar:
        print( 'Computing projections only for Clar resonators' )
    if clarmax:
        print( 'Computing projections only for Clar structures' )

    if kekule:
        # Generate all possible Kekule structures:
        kekFileName = basename + '.kek'
        if os.path.exists( kekFileName ):
            print( 'Kekule structures have been previously enumerated in '
                    'file', kekFileName )
        else:
            elem = []
            if huckel:
                elem = 'C' * hmoSol.NAt
            else:
                for Z in FchkInfo.Z:
                    elem.append( str(Z) )
            with open( kekFileName, 'w' ) as writer:
                print( 'Enumerating Kekule structures for %s ...' % basename )
                if huckel:
                    enumKekule( writer, elem, hmoSol.xyz )
                else:
                    enumKekule( writer, elem, FchkInfo.xyz, FchkInfo.spin, \
                            FchkInfo.mult )
                print( 'File', kekFileName, 'written' )

    if clar:
        # Generate all possible Clar structures:
        clarFileName = basename + '.clar'
        if os.path.exists( clarFileName ):
            print( 'Clar resonators have been previously enumerated in '
                    'file', clarFileName )
        else:
            elem = []
            if huckel:
                elem = 'C' * hmoSol.NAt
            else:
                for Z in FchkInfo.Z:
                    elem.append( str(Z) )
            with open( clarFileName, 'w' ) as writer:
                print( 'Enumerating Clar resonators for %s ...' % basename )
                if huckel:
                    enumClar( writer, elem, hmoSol.xyz, False, atomsexcl )
                else:
                    enumClar( writer, elem, FchkInfo.xyz, False, atomsexcl )
                print( 'File', clarFileName, 'written' )

    if clarmax:
        # Generate all possible Clar structures:
        clarmaxFileName = basename + '.clarmax'
        # In the case of excluding some atoms:
        if len( atomsexcl ) > 0:
            str_atomsexcl = '_pos'
            for ix in atomsexcl:
                str_atomsexcl += '_%i' % ix
            clarmaxFileName = basename + str_atomsexcl + '.clarmax'
        clarFileName = re.sub( '\.clarmax$', '.clar', clarmaxFileName )
        if os.path.exists( clarmaxFileName ):
            print( 'Clar structures have been previously enumerated in '
                    'file', clarmaxFileName )
        else:
            elem = []
            if huckel:
                elem = 'C' * hmoSol.NAt
            else:
                for Z in FchkInfo.Z:
                    elem.append( str(Z) )
            with open( clarmaxFileName, 'w' ) as writer:
                # 1. If *.clar file exists:
                if os.path.exists( clarFileName ):
                    print( 'Extracting Clar structures from file %s ...' % 
                            clarFileName )
                    clarMax_from_clar( clarmaxFileName, atomsexcl )
                    print( 'File', clarmaxFileName, 'written' )
                # 2. Request for direct generation of Clar structures:
                else:
                    print( 'Performing combinatorial enumeration of '
                            'Clar structures ...' )
                    with open( clarmaxFileName, 'w' ) as writer:
                        print( 'Enumerating Clar structures for %s ...' % 
                                basename )
                        if len( atomsexcl ) > 0:
                            print( '--> Excluding atoms:', end='' )
                            for ix in atomsexcl:
                                print( ' %i' % ix, end='' )
                            print()
                        if huckel:
                            enumClar( writer, elem, hmoSol.xyz, True, \
                                    atomsexcl, clarmaxFileName )
                        else:
                            enumClar( writer, elem, FchkInfo.xyz, True, \
                                    atomsexcl, clarmaxFileName )
                        print( 'File', clarmaxFileName, 'written' )

    print()

    # proj_DM_WF_hmo:
    if huckel:
        print( 'In the framework of simple Huckel molecular orbitals theory' )
        if len( lewis ) == 0 and not kekule and not clar and not clarmax:
            if len( lewis ) == 0:
                lewis = -1
            proj_DM_WF_hmo( hmoSol, lewis )
        elif kekule:
            proj_DM_WF_kekule_hmo( hmoSol, kekFileName, projdiv )
        elif clar:
            proj_WF_clar_hmo( hmoSol, clarFileName, 1, clarstyle, sciform )
        elif clarmax:
            proj_WF_clar_hmo( hmoSol, clarmaxFileName, 1, clarstyle, sciform )
        exit(0)

    # proj_DM_WF:
    if len( lewis ) == 0 and not kekule and not clar and not clarmax:
        proj_DM_WF( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
                -1, 'uni', rao_tag, flipraos )
    # proj_DM_WF_kekule:
    elif kekule:
        proj_DM_WF_kekule( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, \
                lmos, kekFileName, projdiv, 'uni', rao_tag, flipraos )
    # proj_WF_clar:
    elif clar:
        proj_WF_clar( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
                 clarFileName, 1, clarstyle, sciform, 'uni', rao_tag, \
                 flipraos )
    elif clarmax:
        proj_WF_clar( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
                 clarmaxFileName, 1, clarstyle, sciform, 'uni', rao_tag, \
                 flipraos )
    # proj_DM_WF_spec:
    else:
        proj_DM_WF( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
            lewis, 'uni', rao_tag, flipraos )
#==============================================================================
# enddef runJob_PROJ()




#==============================================================================
# Perform a SCERED job
#
def runJob_SCERED( basename, kekule, clar, clarmax, atomsexcl ):
    if kekule:
        print( 'Performing SCERED calculations using Kekule structures ...' )
    elif clar:
        print( 'Performing SCERED calculations using Clar resonators ...' )
    else:
        print( 'ERROR: The SCERED job can only be performed '
                'when either the Kekule or the Clar option is turned on' )
        print( 'Aborted' )
        exit(1)


    elem, xyz = readxyz( basename + '.xyz' )
    NAt = len( elem )
    # Only considering C atoms:
    xyz0 = xyz.copy()
    at = [] # C atoms
    xyz = [] 
    for i in range( NAt ):
        if elem[i] == '6' or elem[i] == 'C':
            at.append( i+1 )
            xyz.append( xyz0[i] )
    xyz = np.array( xyz ) # Convert list to numpy array

    # Enumeration of Kekule/Clar resonators:
    if kekule:
        # Generate all possible Kekule structures:
        kekFileName = basename + '.kek'
        if os.path.exists( kekFileName ):
            print( 'Kekule structures have been previously enumerated in '
                    'file', kekFileName )
        else:
            with open( kekFileName, 'w' ) as writer:
                print( 'Enumerating Kekule structures for %s ...' % basename )
                enumKekule( writer, elem, xyz )
                print( 'File', kekFileName, 'written' )
    if clar:
        # Generate all possible Clar structures:
        clarFileName = basename + '.clar'
        if os.path.exists( clarFileName ):
            print( 'Clar resonators have been previously enumerated in '
                    'file', clarFileName )
        else:
            with open( clarFileName, 'w' ) as writer:
                print( 'Enumerating Clar resonators for %s ...' % basename )
                enumClar( writer, elem, xyz, True, atomsexcl )
                print( 'File', clarFileName, 'written' )
    print()

    # Calculate SCERED^2 for Kekule/Clar resonators:
    if kekule:
        calcSCERED2( kekFileName, xyz, at )
#==============================================================================
# enddef runJob_SCERED()




#==============================================================================
#  THE MAIN PROGRAM:
#==============================================================================

# Welcome message:
print( 'EzReson version', VERSION )
print( '  -- A program for resonance analysis of a DFT wave function' )
print( 'Written by Yang WANG (yangwang@yzu.edu.cn)' )
print( 'Copyright 2022 Yang Wang' )
print( '' )

if len(sys.argv) < 2:
    print( 'Usage: python ezreson.py <input-file>' )
    exit(1)

#-------- Parse the input file: --------
inputFile = sys.argv[1] # Get the input file name from command-line
# print( inputFile )
# Check if the input file exists:
if not os.path.isfile( inputFile ):
    print( 'Control file', inputFile, 'not found\nAborted' )
    exit(1)
p = readControlFile( inputFile )


#======== GLOBAL VARIABLES ========
globals.basename = p.basename
globals.ifWrite_Overlap = p.overlap
globals.ifWrite_Hamilt = p.hamilt
globals.ifWrite_Eint = p.eint


#-------- ENUM job --------
if p.job == 'ENUM':
    runJob_ENUM( p.basename, p.kekule, p.clar, p.clarmax, p.atomsexcl )
    printElapsedTime()
    exit(0)

#-------- LMO job --------
if p.job == 'LMO':
    runJob_LMO( p.basename )
    printElapsedTime()
    exit(0)

#-------- WFRT job --------
elif p.job == 'WFRT':
    runJob_WFRT( p.basename, p.lmos, p.atoms, p.maxnlp, p.projcut, \
            p.writeraos, p.flipraos, p.lewis, p.kekule, p.clar, p.clarmax, \
            p.huckel, p.clarstyle, p.atomsexcl, p.sciform )
    printElapsedTime()
    exit(0)

#-------- DMRT job --------
elif p.job == 'DMRT':
    runJob_DMRT( p.basename, p.lmos, p.atoms, p.maxnlp, p.precdmrt, \
            p.degcridmrt, p.lewis, p.kekule, p.huckel )
    printElapsedTime()
    exit(0)

#-------- PROJ job --------
elif p.job == 'PROJ':
    runJob_PROJ( p.basename, p.lmos, p.atoms, p.maxnlp, \
            p.writeraos, p.flipraos, p.lewis, p.kekule, p.clar, p.clarmax, \
            p.projdiv, p.huckel, p.clarstyle, p.atomsexcl, p.sciform )
    printElapsedTime()
    exit(0)

#-------- SCERED job --------
elif p.job == 'SCERED':
    runJob_SCERED( p.basename, p.kekule, p.clar, p.clarmax, p.atomsexcl )
    printElapsedTime()
    exit(0)

##-------- EMRT job --------
#elif p.job == 'EMRT':
#    runJob_EMRT( p.basename, p.lmos, p.atoms, p.maxnlp, p.projcut, \
#            p.writeraos, p.flipraos, p.lewis, p.kekule, p.clar, p.clarmax, \
#            p.huckel, p.clarstyle, p.atomsexcl, p.sciform )
#    exit(0)

#------ Invalid job type --------
else:
    print( 'ERROR: unrecognized job type', p.job )
    print( '  Valid options are ENUM, LMO, WFRT, DMRT, PROJ, SCERED' )
    print( 'Aborted' )
    printElapsedTime()
    exit(1)


#==============================================================================
