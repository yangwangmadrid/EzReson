# EzReson: An efficient program for chemical resonance analysis

The program provides a general and efficient tool to perform a resonance analysis of molecules. It translates the wave function obtained from a standard DFT or Hartree-Fock calculation to the chemically more intepretable Lewis structures.

## How to cite
If you are using EzReson in your research papers or presentations, it is obligatory to cite the following works:
1. Yang Wang. Toward a Reliable and Widely Applicable Resonance Theory from DFT Calculations:
  Should it be Density Based or Wave Function Based? Submitted.
2. Yang Wang. An Efficient Wave-Function-Based Resonance Theory in the DFT Framework. Submitted.

## Copyright and license
The Author of the EzReson software is Yang Wang (yangwang@yzu.edu.cn; orcid.org/0000-0003-2540-2199). The EzReson program is released under GNU General Public License v3 (GPLv3).

## Disclaimer
The EzReson software is provided as it is, with no warranties. The Author shall not be liable for any use derived from it.
Feedbacks and bug reports are always welcome (yangwang@yzu.edu.cn). However, it is kindly reminded that the Author does not take on the responsibility of providing technical support.  


## How to install

### Requirements
- python >= 3.6
- numpy >= 1.18.0
- scipy >= 1.5.1

### Installation
1. Put the folder the EzReson package to any location as you like, which is referred to as the source directory hereafter. 
2. Open with a text editor the script file ezreson.sh in the source directory and set the EZREON_DIR variable as the path of the source directory.
3. Add the source directory to the global environment variable PATH in e.g., .bash_profile or .bashrc under your HOME.

Then, you are ready to go.


## How to use

### Gaussian calculations
1. In the input file (e.g., abc.gjf), add in the route section the keywords "fchk=All Pop=NBO6Read" while at the end of the file add "$NBO NOBOND AONAO=W $END". In this way, the checkpoint file "Test.FChk" and the NBO matrix file "abc.33" will be generated. Then, rename "Test.FChk" to "abc.fchk". The Gaussian output file should have the extension name of ".out" (If necessary, "abc.log" ought to be renamed as "abc.out").

Notes: 
i) It is strongly recommended to use the NBO program later than version 5.0. The free version of NBO 3.1 implemented in Gaussian package would be problematic and give unreliable results.
ii) Do not use "fchk=All" to generate the checkpoint file if there are two such jobs running at the same working directory at the same time. Otherwise, the two jobs will write the same "Test.FChk" file. Instead, add the "%chk=abc.chk" to obtain the checkpoint file "abc.chk". Then, use the formchk utility to convert "abc.chk" to "abc.fchk".

2. Make sure that the following four files, as the inputs for EzReson, are in the same working directory:
    - abc.gjf
    - abc.out
    - abc.fchk
    - abc.33 

3. Change to the working directory and prepare the input file for EzReson (vide infra), say, "abc_wfrt.in"

4. Simply execute the following command:
    ezreson abc_wfrt.in > abc_wfrt.out

As you see, you will find the result of resonance analysis by EzReson in file "abc_wfrt.out".


### A typical wave-function-based resonance theory (WFRT) analysis

1. Perform LMO

### WFRT analysis using Lewis structures with maximum number of lone pairs

### WFRT analysis using projection cutoff to truncate the set of Lewis structures

### WFRT analysis using specified Lewis structures

### WFRT analysis using all Kekule structures

### WFRT analysis in the framework of simple Hueckel molecular orbital (HMO) theory

### A typical density-matrix-based resonance theory (DMRT) analysis

### DMRT analysis in the HMO framework


### Projection calculations

resemble



## Cautions

RAOs
ozone


## Limitations


