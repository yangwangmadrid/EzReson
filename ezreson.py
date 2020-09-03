# EzReson:
#  A general efficient package for performing chemical resonance analysis of a 
#  DFT wave function
#
#  Created on Aug 26, 2020
#  Written by Yang Wang (yangwang@yzu.edu.cn)
#

VERSION = '1.0 (Aug 2020)'

import re
import sys
sys.stdout.flush()
import os.path
import numpy as np
from readFchk import *
from readNBOMat import *
from transform import *
from hmo import *
from wfrt import *


# Parameters in the input file:
class param:
    def __init__( self, basename, job, lmos, atoms, maxnlp, projcut, \
            writeraos, flipraos, lewis, kekule, huckel, precdmrt, degcridmrt ):
        self.basename = basename
        self.job = job
        self.lmos = lmos
        self.atoms = atoms
        self.maxnlp = maxnlp
        self.projcut = projcut
        self.writeraos = writeraos
        self.flipraos = flipraos
        self.lewis = lewis
        self.kekule = kekule
        self.huckel = huckel
        self.precdmrt = precdmrt
        self.degcridmrt = degcridmrt
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

    pat_huckel = re.compile( 
            '^\s*Huckel\s*=\s*.*', re.IGNORECASE
            )

    pat_precdmrt = re.compile( 
            '^\s*PrecDMRT\s*=\s*.*', re.IGNORECASE
            )

    pat_degcridmrt = re.compile( 
            '^\s*DegCriDMRT\s*=\s*.*', re.IGNORECASE
            )

    # Default values:
    basename = ''
    job = ''
    lmos = []
    atoms = []
    maxnlp = inf
    projcut = 0.
    writeraos = False
    flipraos = None
    lewis = []
    kekule = False
    huckel = False
    precdmrt = 0.99
    degcridmrt = 1E-3

    iline = 0
    with open( inputFileName ) as f:
        for line in f:
            iline += 1

            if pat_validline.match( line ):
                fields = re.sub( '^.*=', '', line )
                if len( fields.strip() ) == 0:
                    print( 'ERROR: Value missing at line %i:' % iline )
                    print( '    -->%s' % line, end='' )
                    print( 'Abort' )
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
                for s in fields.split():
                    lmos.append( int(s) )
                # print( 'lmos = ', lmos )
                continue

            # Indices of atoms:
            if pat_atoms.match( line ):
                for s in fields.split():
                    atoms.append( int(s) )
                # print( 'atoms = ', atoms )
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
                    print( 'Abort' )
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
                        print( 'Abort' )
                        exit(1)
                else:
                    for si in s:
                        flipraos.append( int(si) )
                continue

            # Indices of the specified Lewis structures:
            if pat_lewis.match( line ):
                LP = []
                BD = []
                for s in fields.split():
                    # If s is an index of Lewis structures:
                    if ':' not in s and '/' not in s:
                        lewis.append( int(s) )
                    else:
                        lp, bd = parseLewis( s )
                        LP.append( lp )
                        BD.append( bd )
                if len(LP) > 0 or len(BD) > 0:
                    lewis = ( LP, BD )
                # print( 'lewis = ', lewis )
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
                    print( 'Abort' )
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
                    print( 'Abort' )
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

    # End of with open( inputFileName ) as f:

    # Check validity of input values:
    # -- No basename specified:
    if not basename:
        print( 'ERROR: ', end='' )
        print( 'Missing \"File = ...\" to specify the basename of the DFT '
               'output files' )
        print( 'Abort' )
        exit(1)

    # -- No job specified:
    if not job:
        print( 'ERROR: ', end='' )
        print( 'Missing \"Job = ...\" to specify the job type' )
        print( '  Valid options are LMO, WFRT, DMRT, PROJ' )
        print( 'Abort' )
        exit(1)

    # -- Replusive parameters:
    if len( lewis ) > 0 and kekule:
        print( 'ERROR: ', end='' )
        print( 'Lewis structures cannot be specified while Kekule option is ' \
                'turned on' )
        print( 'Abort' )
        exit(1)

    return param( basename, job, lmos, atoms, maxnlp, projcut, writeraos, \
            flipraos, lewis, kekule, huckel, precdmrt, degcridmrt )
#==============================================================================
# enddef readControlFile()


#==============================================================================
# Parse a string of Lewis structure
#
def parseLewis( s ):
    lp = []
    bd = []
    f = s.split( '/' )
    for ep in f: # for each electron pair:
        # A lone pair:
        if ':' in ep:
            n = ep.replace( ':', '' )
            if len(n) != 1:
                raise ValueError( 'Invalid Lewis structure:', s )
            if int(n) < 1:
                raise ValueError( 'Invalid Lewis structure:', s )

            lp.append( ep.replace( ':', '' ) )
        elif '-' in ep:
            b = ep.split( '-' )
            if len(b) != 2:
                raise ValueError( 'Invalid Lewis structure:', s )
            b1 = int( b[0] )
            b2 = int( b[1] )
            if b1 < 1 or b2 < 1:
                raise ValueError( 'Invalid Lewis structure:', s )
            bd.append( b1 )
            bd.append( b2 )

    lp = np.array( lp, dtype='i' )
    nbd = len(bd) // 2
    bd = np.array( bd, dtype='i' ).reshape( nbd, 2 )
    return ( lp, bd )
#==============================================================================
# enddef parseLewis()

#==============================================================================
# Perform an LMO job
#
def runJob_LMO( basename ):
    print( 'Performing Pipek-Mezey localization of molecular orbitals ...' )
    FchkInfo = readFchk( basename + '.fchk' )
    # Check if AOs has no linear dependence (so that C is invertible):
    if FchkInfo.NBasIndp < FchkInfo.NBas:
        raise ValueError( 'Linear dependence found in the AOs\nAbort!' )
    naoInfo = readNAOInfo( basename +'.out' )
    CAONAO = readNBOMat( basename + '.33', naoInfo.NBas, naoInfo.NNAO  )
    CNAOLMO, ELMO = genLMOs( FchkInfo, naoInfo, CAONAO, basename )
#==============================================================================
# enddef runJob_LMO()


#==============================================================================
# Perform a WFRT job
#
def runJob_WFRT( basename, lmos, atoms, maxnlp, projcut, writeraos, \
        flipraos, lewis, kekule, huckel ):
    print( 'Performing Wave Function based Resonance Theory (WFRT) '
           'calculations ...' )
    if not huckel:
        FchkInfo = readFchk( basename + '.fchk' )
        # Check if AOs has no linear dependence (so that C is invertible):
        if FchkInfo.NBasIndp < FchkInfo.NBas:
            raise ValueError( 'Linear dependence found in the AOs\nAbort!' )
        naoInfo = readNAOInfo( basename +'.out' )
        CAONAO = readNBOMat( basename + '.33', naoInfo.NBas, naoInfo.NNAO  )
        CNAOLMO, ELMO = genLMOs( FchkInfo, naoInfo, CAONAO, basename )
    else:
        hmoSol = hmo( basename + '.xyz' )

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

    if kekule:
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

    print()

    # wfrt_hmo:
    if huckel:
        print( 'In the framework of simple Huckel molecular orbitals theory' )
        if len( lewis ) == 0 and not kekule:  # Not supported yet
            wfrt_hmo( hmoSol, atoms, maxnlp, projcut )
        elif kekule:
            wfrt_hmo_kekule( hmoSol, kekFileName )
        else:
            wfrt_hmo_spec( hmoSol, lewis )
        exit(0)

    # wfrt:
    if len( lewis ) == 0 and not kekule:
        wfrt( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, maxnlp, \
            projcut, 'uni', rao_tag, flipraos )
    # wfrt_kekule:
    elif kekule:
        wfrt_kekule( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
                 kekFileName, 'uni' )
    # wfrt_spec:
    else:
        wfrt_spec( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
            lewis, 'uni', rao_tag, flipraos )
#==============================================================================
# enddef runJob_WFRT()



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
            raise ValueError( 'Linear dependence found in the AOs\nAbort!' )
        naoInfo = readNAOInfo( basename +'.out' )
        CAONAO = readNBOMat( basename + '.33', naoInfo.NBas, naoInfo.NNAO  )
        CNAOLMO, ELMO = genLMOs( FchkInfo, naoInfo, CAONAO, basename )
    else:
        hmoSol = hmo( basename + '.xyz' )

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
        lewis, kekule, huckel ):
    print( 'Performing wave function and density matrix projection '
           'calculations ...' )
    if not huckel:
        FchkInfo = readFchk( basename + '.fchk' )
        # Check if AOs has no linear dependence (so that C is invertible):
        if FchkInfo.NBasIndp < FchkInfo.NBas:
            raise ValueError( 'Linear dependence found in the AOs\nAbort!' )
        naoInfo = readNAOInfo( basename +'.out' )
        CAONAO = readNBOMat( basename + '.33', naoInfo.NBas, naoInfo.NNAO  )
        CNAOLMO, ELMO = genLMOs( FchkInfo, naoInfo, CAONAO, basename )
    else:
        hmoSol = hmo( basename + '.xyz' )

    # Calculating WF and DM projections:
    if writeraos:
        rao_tag = basename
        print( 'RAOs are to be written in the file %s_RAO.fchk' % rao_tag )
    else:
        rao_tag = ''
    if kekule:
        print( 'Using only Kekule structures to expand the wave function' )
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


    print()

    # proj_DM_WF_hmo:
    if huckel:
        print( 'In the framework of simple Huckel molecular orbitals theory' )
        if len( lewis ) == 0 and not kekule:
            if len( lewis ) == 0:
                lewis = -1
            proj_DM_WF_hmo( hmoSol, lewis )
        elif kekule:
            proj_DM_WF_kekule_hmo( hmoSol, kekFileName )
        exit(0)

    # proj_DM_WF:
    if len( lewis ) == 0 and not kekule:
        proj_DM_WF( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
                -1, 'uni', rao_tag, flipraos )
    # proj_DM_WF_kekule:
    elif kekule:
        proj_DM_WF_kekule( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, \
                lmos, kekFileName, 1, 'uni', rao_tag, flipraos )
    # proj_DM_WF_spec:
    else:
        proj_DM_WF( naoInfo, FchkInfo, CAONAO, CNAOLMO, ELMO, atoms, lmos, \
            lewis, 'uni', rao_tag, flipraos )
#==============================================================================
# enddef runJob_PROJ()




#==============================================================================
#  THE MAIN PROGRAM:
#==============================================================================

# Welcome message:
print( 'EzReson version', VERSION )
print( '  -- A program for resonance analysis of a DFT wave function' )
print( 'Written by Yang WANG (yangwang@yzu.edu.cn)' )
print( 'Copyright 2020 Yang Wang' )
print( '' )

if len(sys.argv) < 2:
    print( 'Usage: python ezreson.py <input-file>' )
    exit(1)

#-------- Parse the input file: --------
inputFile = sys.argv[1] # Get the input file name from command-line
# print( inputFile )
# Check if the input file exists:
if not os.path.isfile( inputFile ):
    print( 'Control file', inputFile, 'not found\nAbort' )
    exit(1)
p = readControlFile( inputFile )

#-------- LMO job --------
if p.job == 'LMO':
    runJob_LMO( p.basename )
    exit(0)

#-------- WFRT job --------
elif p.job == 'WFRT':
    runJob_WFRT( p.basename, p.lmos, p.atoms, p.maxnlp, p.projcut, \
            p.writeraos, p.flipraos, p.lewis, p.kekule, p.huckel )
    exit(0)

#-------- DMRT job --------
elif p.job == 'DMRT':
    runJob_DMRT( p.basename, p.lmos, p.atoms, p.maxnlp, p.precdmrt, \
            p.degcridmrt, p.lewis, p.kekule, p.huckel )
    exit(0)

#-------- PROJ job --------
elif p.job == 'PROJ':
    runJob_PROJ( p.basename, p.lmos, p.atoms, p.maxnlp, \
            p.writeraos, p.flipraos, p.lewis, p.kekule, p.huckel )
    exit(0)

#------ Invalid job type --------
else:
    print( 'ERROR: unrecognized job type', p.job )
    exit(1)

#==============================================================================
