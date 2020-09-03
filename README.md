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
1. In the input file (e.g., abc.gjf), add in the route section the keywords `fchk=All Pop=NBO6Read` while at the end of the file add `$NBO NOBOND AONAO=W $END`. In this way, the checkpoint file "Test.FChk" and the NBO matrix file "abc.33" will be generated. Then, rename "Test.FChk" to "abc.fchk". The Gaussian output file should have the extension name of ".out" (If necessary, "abc.log" ought to be renamed as "abc.out").

**NOTES**: 
- It is strongly recommended to use the NBO program later than version 5.0. The free version of NBO 3.1 implemented in Gaussian package would be problematic and give unreliable results.
- Do not use `fchk=All` to generate the checkpoint file if there are two such jobs running at the same working directory at the same time. Otherwise, the two jobs will write the same "Test.FChk" file. Instead, add the `%chk=abc.chk` to obtain the checkpoint file "abc.chk". Then, use the formchk utility to convert "abc.chk" to "abc.fchk".

2. Make sure that the following four files, as the inputs for EzReson, are in the same working directory:
- abc.gjf
- abc.out
- abc.fchk
- abc.33 

3. Change to the working directory and prepare the input file for EzReson (vide infra), say, "abc_wfrt.in"

4. Simply execute the following command:
    `ezreson abc_wfrt.in > abc_wfrt.out`

As you see, you will find the result of resonance analysis by EzReson in file "abc_wfrt.out".


### A typical wave-function-based resonance theory (WFRT) analysis

Let us take benzene as an example. Suppose that after the DFT calculations of benzene you have obtained the four output files, Ph.gjf, Ph.out, Ph.fchk and Ph.33.

**1. Perform an LMO calculation to obtain Pipek--Mezey localized MOs**

Prepare an input file, named "benzene_lmo.in" as:
```
File = Ph
Job = LMO
```

**NOTE**: The letters in EzReon's input file are *case-insensitive*.

Then, use the following command to run the LMO job:

`ezreson benzene_lmo.in > benzene_lmo.out`

After finishing the LMO calculation correctly, the following output files will
be generated:
- Ph_CNAOLMO.dat
- Ph_ELMO.dat
- Ph_LMO.fchk


**2. Identify the LMOs associated with the resonance subsystem**

In this particular case of benezene, we are to find the occupied LMOs corresponding to the pi-conjugate system.
To this end, open file "Ph_LMO.fchk" with visualization software like JMol or Gabedit. For JMol, after opening Ph_LMO.fchk, type in the script console:
`isosurface mo 21`
and the LMO-#21 (which is the HOMO) will be displayed. You will see that this LMO belongs to the pi-resonance system and is thus selected.
Then, keep on inspecting lower LMOs, #20, #19, ..., until all LMOs belonging to the resonance subsystem have been chosen. Since in this case there are 6 electrons in the resonance subsystem, you only need to identify 3 LMOs.

As for benzene, the resonating LMOs are identified as 19, 20 and 21.

**NOTE**: Other visualization softwares may not support reading *.fchk for visualization of orbitals. But you can always convert *.fchk file to *.cube file by the cubegen utility in the Gaussian suite. Then, a great variety of visualization tools are able to read the *.cube file to visualize molecular orbitals.


**3. Perform the WRFT analysis**

We have determined that the LMOs for resonance subsystem of benzene are #19, 20 and 21. We also see that the involved atoms are the carbon atoms, whose indices are 1, 2, 3, 4, 5 and 6 as indicated in Ph.gjf.

Accordingly, we prepare an input file, named "benzene_wfrt.in" as:
```
File = Ph
Job = WFRT
LMOs = 19 20 21
Atoms = 1 2 3 4 5 6
```

Then, use the following command to run the WFRT job:
 
`ezreson benzene_wfrt.in > benzene_wfrt.out`

**NOTE**: For the indices of atoms, *the order matters* in order to apply Rumer's rule for determination of linearly independent set of Lewis structures. For monocyclic systems, the ordered atoms should form a circle. For other systems, the choice is somewhat arbitrary, but it is recommended that the atoms be disposed to form a circle as much as possible.


**Fragment of the output**:
```
Performing Wave Function based Resonance Theory (WFRT) calculations ...
 
Reading LEWIS_6c_6e.dat for Lewis structures ...
215 Lewis structures in total to be considered
... ...
... ...
Solving linear equation for WFRT coefficients ...
40 Lewis structures violating Rumer's rule excluded
0 Lewis structures of small projections ( < 0.00E+00 ) excluded
Using totally 175 Lewis structures to expand the actual wave function
Normalizing the wave function by a factor of 0.9951914286 ...
---------------------------------------------------------------------------------------------------------
  No. Projection Coefficient      RE  REeff  Mulli. Bickel. Ros-Sc.  Lowdin    PWSO  Lewis structure
---------------------------------------------------------------------------------------------------------
    1   0.559795   0.3333335  175.96  -31.51  18.75%  22.85%   9.60%  14.59%  19.80%  1-2 3-4 5-6
    2   0.559795   0.3333335  175.96  -31.51  18.75%  22.85%   9.60%  14.59%  19.80%  2-3 4-5 1-6
   15   0.221154   0.2222227  347.79   -8.30   4.94%   2.70%   4.27%   2.77%   2.44%  4: 2-3 5-6
   16   0.221154   0.2222227  347.79   -8.30   4.94%   2.70%   4.27%   2.77%   2.44%  1: 2-3 5-6
   17   0.221154   0.2222219  347.80   -8.30   4.94%   2.70%   4.27%   2.77%   2.44%  6: 1-2 4-5
   18   0.221154   0.2222219  347.80   -8.30   4.94%   2.70%   4.27%   2.77%   2.44%  3: 1-2 4-5
   19   0.221154   0.2222219  347.80   -8.30   4.94%   2.70%   4.27%   2.77%   2.44%  2: 3-4 1-6
   20   0.221154   0.2222219  347.80   -8.30   4.94%   2.70%   4.27%   2.77%   2.44%  5: 3-4 1-6
... ...
... ...
---------------------------------------------------------------------------------------------------------
Reproducibility =  99.519143%
```

### WFRT analysis using Lewis structures with maximum number of lone pairs

Use the option `MAXNLP = <num>` in the input file to set the maximum number of lone pairs allowed in the Lewis structures

**NOTE**: `<num>` can be 0, 1, 2, ..., or inf 


### WFRT analysis using projection cutoff to truncate the set of Lewis structures

For instance, if you want to employ only the Lewis structures whose projection to the actual wave function is less than 0.1, just add the following line in the input file:
```
ProjCut = 0.1
```


### WFRT analysis using specified Lewis structures

You can perform the WFRT analysis using only a user-defined Lewis structures, which can be specified in the input file by using the `Lewis` keyword. There are two ways of specifying Lewis structures.

  - Specification with indices of Lewis structures:

    Here is an example:
    ```
    Lewis = 1 2 9 10
    ```
    It means that only the 1st, 2nd, 9th and 10th Lewis structures are to be used to expand the wave function. Note that the indices follow the descending order of the projection of the Lewis structures onto the actual wave function.

  - Specification by writing explicitly the Lewis structures

    An example:
    ```
    Lewis = 1-2/3-4/5-6  1:/2-3/5: 2:/4:/6:
    ```
    In this example, we have defined three Lewis structures. The first one is a Kekule structure, where `-` denotes the bond connect two atoms. In the second Lewis structure, each of atoms 1 and 5 has a lone pair and atoms 2 and 3 forms a bond. The last Lewis structure contains three lone pairs located, respectively at atoms 2, 4 and 6.


### WFRT analysis using all Kekule structures

To perform a WFRT analysis using only all possible Kekule structures, set the `Kekule` parameter to be TRUE:
```
Kekule = TRUE
```

**NOTE**: The two options `Kekule` and `Lewis` are repulsive to each other, i.e., they cannot be present at the same time in the input file.


### WFRT analysis in the framework of simple Hueckel molecular orbital (HMO) theory

### A typical density-matrix-based resonance theory (DMRT) analysis

### DMRT analysis in the HMO framework


### Projection calculations

resemble



## Cautions

RAOs
ozone


## Limitations


