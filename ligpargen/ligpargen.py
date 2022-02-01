'''
Title:  LigParGen 2.1
Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@yale.edu // israel.cabezadevaca@icm.uu.se  
Place:  William L. Jorgensen Lab (Yale University) // Jens Carlsson Lab (Uppsala University)
Date:   2020-2021

Description: An automatic OPLS-AA parameter generator for small organic using CM1A, 1.14CM1A and CM1A-LBCC 

Notes: Based on the Leela Dodda ligpargen python code.

'''


import logging
import os
import shutil
import importlib.util

from .topology.Molecule import Molecule

from .tools import commandline
from .tools import utilities
from .tools import boss
from .tools import alchemify

from .inout import zmat
from .inout import Q
from .inout import gromacs
from .inout import charmm
from .inout import tinker
from .inout import openmm
from .inout import pqr
from .inout import desmond
from .inout import lammps
from .inout import xplor

import dataclasses

import logging
logger = logging.getLogger(__name__)


__author__ = "Israel Cabeza de Vaca Lopez"
__copyright__ = "Copyright (c) 2020, Israel Cabeza de Vaca Lopez"
__credits__ = ["Israel Cabeza de Vaca Lopez"]
__license__ = "MIT License"
__version__ = "2.1"
__maintainer__ = "Israel Cabeza de Vaca Lopez"
__email__ = "israel.cabezadevaca@icm.uu.se"
__status__ = "Production"

successFinalNote = '''            
    
                    LigParGen finished succesfully!
    
    Please do not forget to cite the following references: 

        1. LigParGen web server: an automatic OPLS-AA parameter generator for organic ligands  
           Leela S. Dodda  Israel Cabeza de Vaca  Julian Tirado-Rives William L. Jorgensen 
           Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W331–W336

        2. 1.14*CM1A-LBCC: Localized Bond-Charge Corrected CM1A Charges for Condensed-Phase Simulations
    '''


@dataclasses.dataclass
class LigParGen():
    """
    Class LigParGen (OPLS CM1A/CM1A-LBCC parameter Generator) to alchemical templates generation.
    Author: Israel Cabeza de Vaca Lopez
    Emails: israel.cabezadevaca@icm.uu.se // israel.cabezadevaca@yale.edu

    Parameters/ Atributes:
    
        ifile : str 
        
            Input file name in PDB, MOL, MOL2, .... or any BABEL accepted format 
        
        ifileB : str
            
            Input file name molecule B in PDB, MOL, MOL2, .... or any BABEL accepted format
        
        smile : str
        
            Input SMILE
        
        smileB : str
        
            Input SMILE molecule B

        molname : str
        
            Name for the output files 

        resname : str
        
            Residue name for the output files 

        numberOfOptimizations : int
        
            Number of optimizations perfomed to molecule A
        
        numberOfOptimizationsB : int
        
            Number of optimizations perfomed to molecule B

        charge : int
        
            Molecule charge - if not provided is autodetected
        
        chargeB : int
        
            Molecule B charge - if not provided is autodetected

        chargeModel : str
        
            Charge model applied to molecule A determine the charges (CM1A, CM1A-LBCC)
        
        chargeModelB : str
        
            Charge model applied to molecule B determine the charges (CM1A, CM1A-LBCC)

        workdir : str
        
            Working directory path 

        debug : bool
        
            if true, more output is provided and temporal files are not destroyed


    """

    ifile:      str = None
    ifileB:     str = None
    smile:      str = None
    smileB:     str = None

    molname:    str = None
    resname:    str = None

    numberOfOptimizations:  int = None
    numberOfOptimizationsB: int = None

    charge:     int = None
    chargeB:    int = None

    chargeAlgorithm:     str = None
    chargeAlgorithmB:    str = None

    workdir: str = None

    debug:   bool = False

    def __post_init__(self):

        self.workdir = utilities.generateFolder(self.workdir)

        molnameA = utilities.guessMoleculeName(self.ifile, self.molname, 'A')

        self.generateLogger(molnameA)

        self.generateSummaryInputs()

        self.alchemicalTransformation = self.checkArguments()

        self.checkRequiredExecutables()
        

        if self.ifile != None:  self.ifile = os.path.abspath(self.ifile)
        if self.ifileB != None: self.ifileB = os.path.abspath(self.ifileB)

        moleculeRDkit, newIndexToOriginalIndex, atomsNameOriginal, residueNameOriginal = utilities.generateRDkitMolecule(self.ifile, self.smile, self.workdir, molnameA, self.debug)

        self.charge = utilities.guessFormalMoleculeCharge(moleculeRDkit, self.charge)

        self.chargeAlgorithm = utilities.guessChargeAlgorithm(self.chargeAlgorithm, self.charge)

        self.resname = utilities.guessResname(residueNameOriginal, self.resname)
        
        moleculeA = Molecule.fromRDkitMolecule(moleculeRDkit)

        logger.info('--- SUMMARY MOLECULE STRUCTURE ---')
        moleculeA.report()

        zmatName = zmat.write(moleculeA, molnameA, self.workdir,  writeAtomParameters = True)

        outFile, pdbFile = boss.run(zmatName, self.chargeAlgorithm, self.numberOfOptimizations, self.charge, molnameA, 
            self.workdir, self.debug)

        moleculeA = Molecule.fromBOSS(zmatName, outFile, pdbFile, moleculeA.shiftX, moleculeA.shiftY, moleculeA.shiftZ)

        utilities.fixNonIntegerCharge(moleculeA)

        self.updateOriginalAtomIndexesAndSerials(moleculeA, newIndexToOriginalIndex, atomsNameOriginal)
        moleculeA.residueName = self.resname

        moleculeA_dual = moleculeA  # Fake possible dual topology for NAMD and others
        molnameAB_dual = molnameA

        if not self.debug:

            os.remove(outFile)
            os.remove(pdbFile)

        if self.alchemicalTransformation:
    
            molnameB = utilities.guessMoleculeName(self.ifileB, self.molname, 'B')

            moleculeRDkitB, newIndexToOriginalIndexB, atomsNameOriginalB, _ = utilities.generateRDkitMolecule(self.ifileB,self.smileB, self.workdir, molnameB)

            BtoAserialCorrespondency,umatchB = alchemify.alignMolecules(moleculeRDkit, moleculeRDkitB, self.workdir, self.debug)

            moleculeB = Molecule.fromRDkitMolecule(moleculeRDkitB)

            zmatNameB = zmat.write(moleculeB, molnameB, self.workdir, writeAtomParameters = True)

            outFileB, pdbFileB = boss.run(zmatNameB, self.chargeAlgorithmB, self.numberOfOptimizationsB, self.chargeB, molnameB, 
                                            self.workdir, self.debug)

            moleculeB = Molecule.fromBOSS(zmatNameB, outFileB, pdbFileB, moleculeB.shiftX, moleculeB.shiftY, moleculeB.shiftZ)

            utilities.fixNonIntegerCharge(moleculeB)

            self.updateOriginalAtomIndexesAndSerials(moleculeB, newIndexToOriginalIndexB, atomsNameOriginalB)

            moleculeAB_singleTopology = alchemify.generateMoleculeAB_singleTopology(moleculeA, moleculeB, BtoAserialCorrespondency, umatchB)

            moleculeAB_dualTopology = alchemify.generateMoleculeAB_dualTopology(moleculeA, moleculeB)
            
            molnameAB, molnameAB_dual = self.generateMoleculeABName(molnameA, molnameB)

            zmat.write(moleculeAB_singleTopology, molnameAB, self.workdir, writeAtomParameters = True)

            zmat.write(moleculeAB_dualTopology, molnameAB_dual, self.workdir, writeAtomParameters = True)

            moleculeA = moleculeAB_singleTopology
            moleculeA_dual = moleculeAB_dualTopology
            self.molname = molnameAB

            if not self.debug:
                
                os.remove(outFileB)
                os.remove(pdbFileB)
                if molnameAB != molnameB: os.remove(zmatNameB)
                if molnameAB != molnameA: os.remove(zmatName)

        else: self.molname = molnameA

        self.moleculeA = moleculeA
        self.moleculeA_dual = moleculeA_dual

    
    def generateMoleculeABName(self, molnameA, molnameB):
        """
        Generate a file name for the alchemical molecule

        Parameters
        ----------
        molnameA : str
            Molecule name
        molnameB : str
            Molecule name

        Returns
        -------
        str, str
            Molecule AB name for single and dual topology
        """

        if self.molname != None:

            molnameAB = self.molname
            molnameAB_dual = self.molname + '_dual'

        elif self.ifile == None and self.ifileB == None:

            molnameAB = 'mol_AB'
            molnameAB_dual = 'mol_AB_dual'

        else:

            if self.ifile == None: molnameA = molnameA
            if self.ifileB == None: molnameB = molnameB

            molnameAB = '_'.join([molnameA, molnameB, 'A_B'])
            molnameAB_dual = '_'.join([molnameA, molnameB, 'A_B_dual'])

        return molnameAB, molnameAB_dual


    def updateOriginalAtomIndexesAndSerials(self, molecule, newIndexToOriginalIndex, atomsNameOriginal):
        """
        Updates the original atom indexes to keep the original atom order

        Parameters
        ----------
        molecule : Molecule object
            Molecule object
        newIndexToOriginalIndex : dict
            Dictionary of new atom indexes and original indexes
        atomsNameOriginal : dict
            Dictionary of atom indexes to original atom names
        """

        for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]:

            atom.indexOriginal = newIndexToOriginalIndex[atom.index-molecule.numberOfStructuralDummyAtoms]
            atom.serialOriginal = atom.indexOriginal + 1 + molecule.numberOfStructuralDummyAtoms

            if len(atomsNameOriginal)>0:

                if len(atomsNameOriginal) != len(set(atomsNameOriginal)):

                    logger.warning('Duplicate names in input file found! Using new generated names')

                else:

                    atom.nameOriginal = atomsNameOriginal[atom.indexOriginal]


    def generateSummaryInputs(self):
        """
        Print out the summary log information
        """

        if self.ifile != None: logger.info('INPUT FILE: '+ self.ifile)
        if self.smile != None: logger.info('INPUT SMILE: '+ self.smile)

        if self.ifileB != None: logger.info('INPUT FILE B: '+ self.ifileB)
        if self.smileB != None: logger.info('INPUT SMILE B: '+ self.smileB)

        if self.numberOfOptimizations != None: 
            logger.info('NUMBER OF OPTIMIZATIONS: '+ str(self.numberOfOptimizations))

        if self.numberOfOptimizationsB != None: 
            logger.info('NUMBER OF OPTIMIZATIONS in B: '+ str(self.numberOfOptimizationsB))


    def checkArguments(self):
        """
        Function to check if arguments are right

        Returns
        -------
        bool
            True if an alchemical transformation will be performed
        """

        if self.ifile != None and not os.path.isfile(self.ifile): 
            logger.error('Input file: '+self.ifile+' does NOT exist')
            exit()


        if self.ifileB != None and not os.path.isfile(self.ifileB): 
            logger.error('Input file: '+self.ifileB+' does NOT exist')
            exit()

        if self.smile==None and self.ifile==None:

            logger.error("       ERROR!    Missing input molecule.")
            exit()

        if self.smile!=None and self.ifile != None:
            logger.error("       ERROR!    SMILE and INPUT file provided. Choose one!!!")
            exit()

        if self.smileB!=None and self.ifileB != None:
            logger.error("       ERROR!    SMILE molecule B and INPUT molecule B file provided. Choose one!!!")
            exit()


        alchemicalTransformation = False

        if self.smileB!=None or self.ifileB!=None: alchemicalTransformation = True

        if alchemicalTransformation:
            
            if self.chargeB==None: self.chargeB = self.charge
            if self.chargeAlgorithmB==None: self.chargeAlgorithmB = self.chargeAlgorithm
            if self.numberOfOptimizationsB==None: self.numberOfOptimizationsB = self.numberOfOptimizations

        return alchemicalTransformation



    def checkRequiredExecutables(self):
        """
        Check is required executables are installed and accesible from the system
        """

        babel = shutil.which("obabel")

        if babel==None:

            print("\n       ERROR!    Open Babel is not installed on your work station.\n")
            exit()

        csh = shutil.which("csh")
        if csh==None: 
            print("\n       ERROR!    csh is not installed on your work station.\n")
            exit()

        try:  os.environ["BOSSdir"]
        except KeyError: 

            print('''\n         ERROR!    BOSS is not found on your work station. $BOSSdir path is not set.
    BOSS software can be downloaded for free from  ---> http://zarbi.chem.yale.edu/software.html 

        DO NOT forget to export the BOSSdir path (export BOSSdir=/path/boss/intalled)
            ''')
            exit()

        if importlib.util.find_spec('rdkit') == None:

            print("\n       ERROR!    RDkit is not installed and it is needed for alchemical transformations.\n")
            exit()


    def writeAllOuputs(self):
        """
        Write molecule in the different output files format
        """

        logger.info('Generating Q inputs')
        Q.write(self.moleculeA, self.molname, self.workdir)

        if self.moleculeA_dual != self.moleculeA: Q.write(self.moleculeA_dual, self.molname + '_dual', self.workdir)

        logger.info('Generating Gromacs inputs')
        gromacs.write(self.moleculeA, self.molname, self.workdir)

        logger.info('Generating CHARMM/NAMD inputs')
        charmm.write(self.moleculeA_dual, self.molname, self.workdir)

        logger.info('Generating TINKER inputs')
        tinker.write(self.moleculeA, self.molname, self.workdir)

        if self.moleculeA_dual != self.moleculeA: tinker.write(self.moleculeA_dual, self.molname + '_dual', self.workdir)

        if not self.moleculeA.alchemicalTransformation:

            logger.info('Generating openMM inputs')
            openmm.write(self.moleculeA, self.molname, self.workdir)

            logger.info('Generating PQR inputs')
            pqr.write(self.moleculeA, self.molname, self.workdir)

            logger.info('Generating DESMOND inputs')
            desmond.write(self.moleculeA, self.molname, self.workdir)

            logger.info('Generating LAMMPS inputs')
            lammps.write(self.moleculeA, self.molname, self.workdir)

            logger.info('Generating XPLOR inputs')
            xplor.write(self.moleculeA, self.molname, self.workdir)


    def generateLogger(self, molname):
        """
        Generate a Logger from the logging module to print the log outputs

        Parameters
        ----------
        molname : str
            Molecule name

        Returns
        -------
        logging
            Logger generated to print the logs
        """


        FORMATTER = logging.Formatter("%(asctime)s — %(levelname)s — %(message)s", 
            datefmt='%Y-%m-%d %H:%M:%S')

        LOG_FILE = os.path.join(self.workdir, molname+'.log')
        file_handler = logging.FileHandler(LOG_FILE,'w')
        file_handler.setFormatter(FORMATTER)

        logger = logging.getLogger()
        logger.addHandler(file_handler)

        logger.setLevel(logging.INFO)


def main():

    args    = commandline.commandline_options()

    moleculeA = LigParGen(ifile=args.ifile, ifileB= args.ifileB, smile= args.smile, smileB= args.smileB,
        charge=args.charge, chargeB=args.chargeB, chargeAlgorithm=args.cgen, chargeAlgorithmB=args.cgenB,
        numberOfOptimizations= args.opt, numberOfOptimizationsB=args.optB, resname=args.resname, 
        molname=args.molname, workdir= args.path, debug= args.debug)

    moleculeA.writeAllOuputs()

    logger.info(successFinalNote)

    print('LigParGen finished succesfully!')


if __name__ == "__main__":

    main()

