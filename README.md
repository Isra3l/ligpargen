# **LigParGen v2.1**

**Author:** &nbsp;&nbsp;Israel Cabeza de Vaca Lopez</br>
**Email:**  &nbsp;&nbsp;&nbsp;israel.cabezadevaca@yale.edu // israel.cabezadevaca@icm.uu.se </br>
**Place:** &nbsp;&nbsp;&nbsp; William L. Jorgensen Lab at Yale University // Jens Carlsson Lab at Uppsala university</br>
**Date:** &nbsp;&nbsp;  2020-2021

**Description:** An automatic OPLS-AA parameter generator for small organic molecules using CM1A, 1.14CM1A and CM1A-LBCC charge models. LigParGen accepts any Open Babel molecular format including SMILES, PDB, MOL, MOL2, among others. Final OPLSAA parameter outputs will be written in topology/coordinate input files for BOSS, Q, Tinker, PQR, openMM, CHARMM/NAMD, Gromacs, LAMMPS, Desmond, and xplor softwares.

**Note:** This new version has been written from scratch but it is based on the **Leela Dodda** initial ligpargen python code.

New LigParGen features:

- The order and the name of the atoms will remain the same in the output files.
- This new version of the ligpargen includes a robust version of the alchemical transformation method to generate single and dual topologies for four different molecular mechanics softwares (BOSS, CHARMM/NAMD, Gromacs, and Tinker).
- Sanity checks to detect incorrect inputs have been implemented (incompatible charge model, net molecule charge, input format, ...).
- A log file to check the inputs, outputs, and intermediate processes has been created. This allows tracking the input information (charge model used, input molecule,...) and also provides useful warnings in case of error.
- Automatic net charge detection in the input molecule. If the user specifies a different charge than one automatically estimated, the log file will include a warnning.
- Atom XYZ positions in the molecule input remain unchanged in the output files.
- Hydrogen atoms will be added automatically just for SMILES inputs. Molecules in any other input format require to have all hydrogens to provide more flexibility of the protonation states. In this way, the user can avoid parameterization problems with tautomers or stereoisomers.
- Different molecule input format, net charges, charge models, and optimization steps can be used for each molecule in alchemical transformations.
- Bugs fixed (Q wrong torsion parameters, alchemical molecule overlap failure,...)
- Additional default input parameters have been included such as residue name, charge model,...

## **INSTALATION**

LigParGen requires the free BOSS software to generate the OPLSAA parameters.

1 - Download and install BOSS software from the official William L. Jorgensen lab website: http://zarbi.chem.yale.edu/software.html

BOSS is compiled for linux using 32 bits libraries so it can not run in windows using WSL. Alternatively, you can use a virtual machine in windows such as virtualBox to install a linux distro (ubuntu, centos, ...) and run LigParGen.

1.1 - Set the BOSSdir enviromental variable:

  - bashrc
  
            export BOSSdir=PATH_TO_BOSS_DIRECTORY
  - cshrc
        
            setenv BOSSdir PATH_TO_BOSS_DIRECTORY

    **TIP:** add this command line in your ~/.bashrc or ~/.cshrc file.

2 - Download and install conda (anaconda or minicoda):  

    wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh

or

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

3 - Create and activate an enviroment for python3.7 (Everything was tested with 3.7):

    conda create --name py37 python=3.7

    conda activate py37  
**Optional TIP:** add this line to your .bashrc/.cshrc,....

4- Install rdkit in py37 enviroment

    conda install -c rdkit rdkit

5 - Download and install LigParGen

    wget https://github.com/Isra3l/ligpargen.git

    pip install -e ligpargen

**Optional:** Check your installation by runing the tests included in the ligpargen folder.

    cd ligpargen;python -m unittest

**TIP:** Do not forget to activate your py37 enviroment before using LigParGen.

    conda activate py37

## **USAGE**

Input arguments list:

- ' -s ': SMILES input  (Ex. -s 'CCCC' )
- ' -i ': PDB, MOL and MOL2 files + any input format supported by Open Babel (Ex. -i benzene.pdb )
- ' -n ': Molecule file name (Ex. -n benzene - it will produce benzene.gmx.gro, benzene.charmm.rtf,...)
- ' -p ': Folder for the output files (Ex. -p bnz )
- ' -r ': Residue name for the output files (Ex. -r BNZ)
- ' -c ': Molecule net charge (Ex. -c +1)
- ' -o ': Number of optimizations (Ex. -o 3)
- ' -cgen ': Charge model to be used - CM1A or CM1A-LBCC (Ex. -cgen CM1A)

For alchemical transformations:

- ' -sb ': SMILES input for molecule B (Ex. -sb 'CCCCO' )
- ' -ib ': PDB, MOL and MOL2 files + any input format supported by Open Babel (Ex. -ib phenol.pdb )
- ' -cb ': Molecule net charge (Ex. -cb +1)
- ' -ob ': Number of optimizations (Ex. -ob 3)
- ' -cgenb ': Charge model to be used - CM1A or CM1A-LBCC (Ex. -cgenb CM1A)

**Notes:**

- CM1A is automatically scaled by 1.14 in neutral molecules.
- CM1A-LBCC is just for neutral molecules and it is also scaled by 1.14.
  
## **Examples**

Molecule template generation:

    ligpargen -i phenol.pdb -n phenol -p phenol -r MOL -c 0 -o 0 -cgen CM1A
    ligpargen -s 'c1ccc(cc1)O' -n phenol -p phenol -r MOL -c 0 -o 0 -cgen CM1A
    ligpargen -s 'c1ccc(cc1)O'
    ligpargen -s 'c1ccc(cc1)O' -n phenol -cgen CM1A-LBCC

Alchemical transformations:

    ligpargen -i phenol.pdb -ib benzene.pdb -n phenolToBenzene -p phenol2bnz -r A2B -c 0 -o 0 -cgen CM1A -cb 0 -ob 1 -cgenb CM1A-LBCC
    ligpargen -i phenol.pdb -ib benzene.pdb -n phenolToBenzene
    ligpargen -s 'c1ccc(cc1)O' -sb 'c1ccccc1' -n phenol_benzene
    ligpargen -s 'c1ccc(cc1)O' -sb 'c1ccccc1' -n phenol_benzene -o 0 -cgen CM1A -cb 0 -ob 1 -cgenb CM1A-LBCC
    ligpargen -i phenol.pdb -sb 'c1ccccc1' -n phenol_benzene -o 0 -cgen CM1A -cb 0 -ob 1 -cgenb CM1A-LBCC

Default values:

- -n : molecule
- -p : Current working directory
- -r : MOL
- -c,-cb : Automatically determined from input molecule
- -o,-ob : 0
- -cgen,-cgenb: CM1A-LBCC for neutral molecules and CM1A for charged molecules.

For help use the -h flag:

    ligpargen -h

Please do not forget to cite the following references:

1. LigParGen web server: an automatic OPLS-AA parameter generator for organic ligands  
        Leela S. Dodda  Israel Cabeza de Vaca  Julian Tirado-Rives William L. Jorgensen 
        Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W331–W336

2. 1.14*CM1A-LBCC: Localized Bond-Charge Corrected CM1A Charges for Condensed-Phase Simulations
