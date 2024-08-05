"""

Molecule Class definition

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

import re
import sys
from ..tools import utilities
from ..topology import Atom, Bond, Angle, Torsion
from ..tools import geometry
from rdkit import Chem
import rdkit.Geometry.rdGeometry as rdGeometry 

import logging
logger = logging.getLogger(__name__)


def _parametersSanityCheck(bossData: list):
    """Sanity check to report BOSS potential issues detected during OPLS parameterization.

    Args:
        bossData (list): BOSS out log file containing the parameterization data
    """

    for line in bossData:

        if 'CHECK' in line: logger.warning(f'Detected potential problem with force parameterization: {line}')



def _generateNewAtom(moleculeRDkit, atom, posParentParentParent, posParentParent, posParent, posAtom, parent,parentParent,parentParentParent, dummyAtomsShift):
    """
    Generate molecule atom from RDkit molecule

    Parameters
    ----------
    moleculeRDkit : RDkit molecule
        RDkit molecule
    atom : int
        Atom index
    posParentParentParent : [float, float, float]
        Position grandGrandparent atom
    posParentParent : [float, float, float]
        Position grandparent atom
    posParent : [float, float, float]
        Position parent atom
    posAtom : [float, float, float]
        Position atom
    parent : int
        Parent Atom
    parentParent : int
        Grand parent atom
    parentParentParent : int
        Grand grand parent atom
    dummyAtomsShift : int
        Serial shift produced by structural dummy atoms

    Returns
    -------
    Atom : Atom class
        Atom object
    """
    
    distance = rdGeometry.Point3D.Distance(posParent, posAtom)
    angle = geometry.getAngle(posParentParent, posParent, posAtom)

    if angle>180.0: 
        logger.warning(f'Detected angle bigger than 180.0 for atom {atom}')
        angle = 180.0

    if posParentParentParent!=None: dihedral = geometry.getDihedral(posParentParentParent, posParentParent, posParent, posAtom)
    else: dihedral = posAtom.z

    element = moleculeRDkit.GetAtomWithIdx(atom).GetSymbol()
    atomicNumber = moleculeRDkit.GetAtomWithIdx(atom).GetAtomicNum()
    atomicMass = Chem.GetPeriodicTable().GetAtomicWeight(atomicNumber)

    return Atom.fromRDkit(atom, atom+dummyAtomsShift,parent,parentParent,parentParentParent,distance,\
        angle, dihedral, element, atomicNumber, atomicMass)


class Molecule(object):
    """
    docstring
    """
    def __init__(self, atoms, bondsVariable, bondsAdditional, anglesVariable, anglesAdditional, torsionsVariable, 
        torsionsAdditional, serial2IndexAtoms, lst14pairs, numberOfStructuralDummyAtoms, shifX, shifY, shifZ):

        self.atoms = atoms

        self.bondsVariable = bondsVariable
        self.bondsAdditional = bondsAdditional

        self.anglesVariable = anglesVariable
        self.anglesAdditional = anglesAdditional

        self.torsionsVariable = torsionsVariable
        self.torsionsAdditional = torsionsAdditional

        self.serial2IndexAtoms = serial2IndexAtoms

        self.lst14pairs = lst14pairs

        self.numberOfStructuralDummyAtoms = int(numberOfStructuralDummyAtoms)

        self.geometryVariations = []

        self.alchemicalTransformation = False
        self.dualTopology = False

        self.excludedList = []

        self.excludedListDict = {}

        self.shiftX = shifX # To restore initial coordinates
        self.shiftY = shifY # To restore initial coordinates
        self.shiftZ = shifZ # To restore initial coordinates

        self.residueName = 'MOL'


    @classmethod
    def fromBOSS(cls, zmatFile, outfile, pdbfile, shifX=0.0, shifY=0.0, shiftZ=0.0):
        """
        Alternative constructor of the molecule class

        Parameters
        ----------
        zmatFile : str
            Zmat file of the molecule
        outfile : str
            BOSS output file name
        pdbfile : str
            BOSS PDB file name
        shifX : float, optional
            X translation of the molecule, by default 0.0
        shifY : float, optional
            Y translation of the molecule, by default 0.0
        shiftZ : float, optional
            Z translation of the molecule, by default 0.0

        Returns
        -------
        Molecule
            Molecule object
        """

        logger.info('Parsing BOSS molecule outputs')

        with open(zmatFile) as f: zmatData = f.read()
        with open(outfile) as f: outFileData = f.read()

        serial2IndexAtoms = {}

        atoms, numberOfStructuralDummyAtoms = cls._getAtoms(cls, zmatData, pdbfile)

        for i,atom in enumerate(atoms): serial2IndexAtoms[atom.serial]=i

        bondsVariable, bondsAdditional = cls._getBonds(zmatData, outFileData, atoms, serial2IndexAtoms)

        cls._generateBondedAtoms(atoms, bondsVariable, bondsAdditional)

        anglesVariable, anglesAdditional = cls._getAngles(zmatData, outFileData, atoms, serial2IndexAtoms)
        torsionsVariable, torsionsAdditional = cls._getTorsions(atoms, zmatData, outFileData, serial2IndexAtoms)

        lst14pairs = cls._get14pairs(torsionsVariable, torsionsAdditional, anglesVariable, anglesAdditional, bondsVariable, bondsAdditional)

        return cls(atoms, bondsVariable, bondsAdditional, anglesVariable, anglesAdditional, torsionsVariable, 
            torsionsAdditional, serial2IndexAtoms, lst14pairs, numberOfStructuralDummyAtoms, shifX, shifY, shiftZ)


    @classmethod
    def fromRDkitMolecule(cls, moleculeRDkit):
        """
        Alternative constructor to generate the molecule object

        Parameters
        ----------
        moleculeRDkit : RDkit molecule
            RDkit molecule

        Returns
        -------
        moleculeObject
            A molecule object with the atributes of the class
        """

        dummyAtomsShift = 4 # dummy atoms (+3) and index to serial (+1)

        logger.info('Parsing Bonds, Angles and dihedrals')

        bondsVariable, bondsAdditional, anglesVariable, anglesAdditional, torsionsVariable, torsionsAdditional = \
            geometry.generateBondAnglesAndDihedrals(moleculeRDkit, dummyAtomsShift)

        logger.info('Parsing Atoms')

        cls._checkIfDiscontinuitiesInBondedLst(bondsVariable)
        cls._checkIfDiscontinuitiesInBondedLst(anglesVariable)
        cls._checkIfDiscontinuitiesInBondedLst(torsionsVariable)

        shiftcoordinates = utilities.translateToceroZcoord(moleculeRDkit)  # BOSS imposes Z coordinate of the first atom in cartesians to be cero (internal coordinates)

        atoms = cls.generateAtomsFromRDkit(moleculeRDkit, bondsVariable, anglesVariable, torsionsVariable, dummyAtomsShift)

        return cls(atoms, bondsVariable, bondsAdditional, anglesVariable, anglesAdditional, torsionsVariable, 
            torsionsAdditional, None, None, 0, shiftcoordinates[0], shiftcoordinates[1], shiftcoordinates[2])


    def generateAtomsFromRDkit(moleculeRDkit, variableBonds, variableAngles, variableDihedrals, dummyAtomsShift):
        """
        Method to generate Atoms from an RDkit molecule

        Parameters
        ----------
        moleculeRDkit : RDkit molecule
            An rdkit molecle object
        variableBonds : list
            Bond list of the molecule containing Bond objects
        variableAngles : list
            Angle list of the molecule containing Angle objects
        variableDihedrals : list
            Torsion list of the molecule containing Torsion objects
        dummyAtomsShift : int
            Serial shift produced by the number of structural dummy atoms

        Returns
        -------
        Atoms list
            List of Atom objects
        """
    
        conf = moleculeRDkit.GetConformer()

        atoms = []

        dummy1 = Atom.fromRDkit(0, 1,0,0,0,0.0,0.0,0.0,'DUM',-1,0.0)
        dummy2 = Atom.fromRDkit(0, 2,1,0,0,1.0,0.0,0.0,'DUM',-1,0.0)
        dummy3 = Atom.fromRDkit(0, 3,2,1,0,1.0,90.0,0.0,'DUM',-1,0.0)

        atoms.append(dummy1)
        atoms.append(dummy2)
        atoms.append(dummy3)

        positionDummy1 = rdGeometry.Point3D(0.0,0.0,0.0)
        positionDummy2 = rdGeometry.Point3D(1.0,0.0,0.0)
        positionDummy3 = rdGeometry.Point3D(1.0,1.0,0.0)

        # Add dummy atoms and build internal coordinates

        variableFirstAtomSet_bond = {bond.atomA.serial: bond for bond in variableBonds}
        variableFirstAtomSet_angle = {angle.atomA.serial: angle for angle in variableAngles}
        variableFirstAtomSet_dihedral = {dihedral.atomA.serial: dihedral for dihedral in variableDihedrals}


        for rdkitAtom in moleculeRDkit.GetAtoms():

            atomAserial = rdkitAtom.GetIdx() + dummyAtomsShift
            atomAindex = rdkitAtom.GetIdx()
            positionAtomA = conf.GetAtomPosition(atomAindex)
            

            if atomAserial in variableFirstAtomSet_dihedral:

                atomB = variableFirstAtomSet_dihedral[atomAserial].atomB
                atomC = variableFirstAtomSet_dihedral[atomAserial].atomC
                atomD = variableFirstAtomSet_dihedral[atomAserial].atomD

                positionAtomB = conf.GetAtomPosition(atomB.index)
                positionAtomC = conf.GetAtomPosition(atomC.index)
                positionAtomD = conf.GetAtomPosition(atomD.index)

                atomBserial = atomB.serial
                atomCserial = atomC.serial
                atomDserial = atomD.serial

            elif atomAserial in variableFirstAtomSet_angle:

                atomB = variableFirstAtomSet_angle[atomAserial].atomB
                atomC = variableFirstAtomSet_angle[atomAserial].atomC
                atomD = dummyAtomsShift-1

                positionAtomB = conf.GetAtomPosition(atomB.index)
                positionAtomC = conf.GetAtomPosition(atomC.index)
                positionAtomD = positionDummy3

                atomBserial = atomB.serial
                atomCserial = atomC.serial
                atomDserial = 3


            elif atomAserial in variableFirstAtomSet_bond:

                atomB = variableFirstAtomSet_bond[atomAserial].atomB
                atomC = dummyAtomsShift - 1
                atomD = dummyAtomsShift - 2

                positionAtomB = conf.GetAtomPosition(atomB.index)
                positionAtomC = positionDummy3
                positionAtomD = positionDummy2

                atomBserial = atomB.serial
                atomCserial = 3
                atomDserial = 2

            else:

                atomB = dummyAtomsShift - 1
                atomC = dummyAtomsShift - 2
                atomD = dummyAtomsShift - 3

                positionAtomB = positionDummy3
                positionAtomC = positionDummy2
                positionAtomD = positionDummy1

                atomBserial = 3
                atomCserial = 2
                atomDserial = 1
                

            atoms.append(_generateNewAtom(moleculeRDkit, atomAindex, positionAtomD, positionAtomC, positionAtomB, positionAtomA,
                                            atomBserial, atomCserial, atomDserial, dummyAtomsShift))


        return atoms

    def _get14pairs(torsionsVariable, torsionsAdditional, anglesVariable, anglesAdditional, bondsVariable, bondsAdditional):
        """
        Generate 14 pair list from the molecule torsions

        Parameters
        ----------
        torsionsVariable : List
            Variable torsion List
        torsionsAdditional : List
            Additional torsion List

        Returns
        -------
        pairLst : list 
            14 pair list
        """

        bondedUpToThree = [[bond.atomA.serial, bond.atomB.serial] for bond in bondsVariable] 
        bondedUpToThree += [[bond.atomA.serial, bond.atomB.serial] for bond in bondsAdditional]
        bondedUpToThree += [[angle.atomA.serial, angle.atomC.serial] for angle in anglesVariable]
        bondedUpToThree += [[angle.atomA.serial, angle.atomC.serial] for angle in anglesAdditional]
        
        bondedUpToThree_sorted = set(tuple(sorted(sublist)) for sublist in bondedUpToThree)

        total = torsionsVariable + torsionsAdditional

        pairSet = set()

        for torsion in total:
            
            if not torsion.improper and torsion.typeInitial!=0:

                pair14 = tuple(sorted([torsion.atomA,torsion.atomD], key= lambda x: x.serial))
                pair14serial = tuple([pair14[0].serial, pair14[1].serial])

                if pair14serial not in bondedUpToThree_sorted:
                
                    pairSet.add(pair14)

        pairLst = [list(item) for item in sorted(pairSet, key= lambda x: (x[0].serial, x[1].serial))]

        return pairLst 


    def _generateAtoms(atomsLines, atomsParameters, atomsLinesPDB):
        """
        Generate the atoms list from the zmat lines

        Parameters
        ----------
        atomsLines : List
            Atom lines in BOSS zmat
        atomsParameters : List
            Atom parameters in the ZMAT (VdW, Charge,..)
        atomsLinesPDB : List
            Atom lines in PDB

        Returns
        -------
        atoms : List
            Atoms list of atoms objects
        numberOfStructuralDummyAtoms : int
            Number of structural dummy atoms

        """
        
        atoms = []
        atomType2atomicNumber = {}
        atomType2atomTypeOPLS = {}
        atomType2charge = {}

        atomNumber2Mass = {}
        atomNumber2Element = {}

        atomType2sigma = {}
        atomType2epsilon = {}

        for atomParameter in atomsParameters:

            tmp = atomParameter.split()

            atomType2atomicNumber[tmp[0]]=tmp[1]
            atomType2atomTypeOPLS[tmp[0]]=tmp[2]
            atomType2charge[tmp[0]]=tmp[3]
            atomType2sigma[tmp[0]]=tmp[4]
            atomType2epsilon[tmp[0]]=tmp[5]

            atomNumber2Mass[tmp[0]]=Chem.GetPeriodicTable().GetAtomicWeight(int(tmp[1]))
            atomNumber2Element[tmp[0]]=Chem.GetPeriodicTable().GetElementSymbol(int(tmp[1]))
        
        numberOfStructuralDummyAtoms = len([line for line in atomsLines if line.split()[2]=='-1'])

        itera = 0

        for i,line in enumerate(atomsLines):

            atomTypeNumber = line.split()[2]

            if atomTypeNumber=='-1': 
                
                x = 0.0
                if i==1 or i==2: x = 1.0
                y = 0.0
                if i==2: y = 1.0
                z = 0.0

                A = Atom.fromZmat(line, 0, 'DU', 
                    0.0, 0.0, 0.0,x,y,z, 
                    0.0, 'DU' )

            else:                

                x = atomsLinesPDB[itera][28:38]
                y = atomsLinesPDB[itera][38:46]
                z = atomsLinesPDB[itera][46:]

                A = Atom.fromZmat(line, atomType2atomicNumber[atomTypeNumber], atomType2atomTypeOPLS[atomTypeNumber], 
                    atomType2charge[atomTypeNumber], atomType2sigma[atomTypeNumber], atomType2epsilon[atomTypeNumber],x,y,z, 
                    atomNumber2Mass[atomTypeNumber], atomNumber2Element[atomTypeNumber] )

                itera+=1

            atoms.append(A)

        return atoms, numberOfStructuralDummyAtoms

    def _checkIfDiscontinuitiesInBondedLst(bondedVariableLst):
        """Check if first atoms in each bonded term are consecutive 

        Parameters
        ----------
        bondedLst : List
            Bonded Term index list 

        """
        # verify that there are NOT discontinuities in variable list

        if len(bondedVariableLst)==0: return

        variableLst = [bondedTerm.atomA.serial for bondedTerm in bondedVariableLst]

        if variableLst != list(range(min(variableLst), max(variableLst)+1)):
            
            logger.error('Discontinuities in Bonded terms!!!')
            logger.error('Variable List: '+' '.join([str(item) for item in variableLst]))
            exit() 

    def _getAtoms(cls, zmatData, pdbfile):
        """
        Generate atoms from ZMAT data and PDB file

        Parameters
        ----------
        zmatData : lst
            zmat molecule data
        pdbfile : str
            PDB file name

        Returns
        -------
        atoms : list
            Atom list in atom object
        numberOfStructuralDummyAtoms : int
            Number of structural dummy atoms (default is three)
        """
        
        atomsLines = re.search(r'BOSS(.*?)Geometry', zmatData, re.DOTALL).group().splitlines()[1:-1]

        atomsParameters = re.search(r'Final Non-B(.*?)$', zmatData, re.DOTALL).group().splitlines()[1:]
        atomsParameters = [line for line in atomsParameters if not len(line.split())==0]
        
        with open(pdbfile) as f: pdbData = f.read()

        atomsLinesPDB = re.search(r'(ATOM|HETATM)(.*?)$', pdbData, re.DOTALL).group().splitlines()
        atomsLinesPDB = [line for line in atomsLinesPDB if 'ATOM' in line or 'HETATM' in line]

        if len(atomsParameters)==0 or len(atomsLinesPDB) ==0:
            logger.error('ERROR: Final Non-Bonded Parameters NOT found in zmat')
            sys.exit() 

        atoms, numberOfStructuralDummyAtoms = cls._generateAtoms(atomsLines, atomsParameters, atomsLinesPDB)

        return atoms, numberOfStructuralDummyAtoms


    def _getBonds(zmatData, outfileData, atoms, serial2IndexAtoms):
        """
        Generate bonds from zmatdata and BOSS out file

        Parameters
        ----------
        zmatData : str
            zmat file data
        outfileData : str
            BOSS out file data
        atoms : lst
            List of atoms
        serial2IndexAtoms : Dict
            Dictionary between serial and index atoms

        Returns
        -------
        bondsVariable : list
            Variable bonds list
        bondsAdditional : list
            Additional bonds list
        """

        additionalBondInfo = re.search(r'Additional Bonds follow(.*?)Harmonic Constraints follow', 
                        zmatData, re.DOTALL).group().splitlines()[1:-1]
        additionalBondInfo = ['_'.join(ele.split()) for ele in additionalBondInfo ]

        try:

            bondsFile = re.search(r'Atom1   Atom2      R0           K0         R1(.*?)Angle Bending Parameters', 
                outfileData, re.DOTALL).group().splitlines()[1:-1]
            
            _parametersSanityCheck(bondsFile)
            
            bondsFile = [bond[:40].split() for bond in bondsFile if not len(bond)==0 and 
                not 'Missing' in bond and 
                not 'Synonym' in bond and
                not 'CHECK' in bond]
            
            

        except:

            bondsFile = []

        bondsVariable = []
        bondsAdditional = []

        for bondF in bondsFile:

            atomA = atoms[serial2IndexAtoms[int(bondF[0])]]
            atomB = atoms[serial2IndexAtoms[int(bondF[1])]]

            A = Bond(atomA, atomB, *bondF[2:])


            if not '_'.join([str(A.atomA.serial),str(A.atomB.serial)]) in additionalBondInfo: bondsVariable.append(A)
            else: bondsAdditional.append(A)

        return bondsVariable, bondsAdditional

    def _getAngles(zmatData, outfileData, atoms, serial2IndexAtoms):
        """
        Generate angles from zmatdata and BOSS out file

        Parameters
        ----------
        zmatData : str
            zmat file data
        outfileData : str
            BOSS out file data
        atoms : lst
            List of atoms
        serial2IndexAtoms : Dict
            Dictionary between serial and index atoms

        Returns
        -------
        anglesVariable : list
            Variable angles list
        anglesAdditional : list
            Additional angles list
        """

    
        additionalAnglesInfo = re.search(r'Additional Bond Angles follow(.*?)Variable Dihedrals follow', zmatData, re.DOTALL).group().splitlines()[1:-1]
        additionalAnglesInfo = ['_'.join(ele.split()) for ele in additionalAnglesInfo ]


        try:
            
            anglesFile = re.search(r'Atom1 Atom2 Atom3    A0       K0       A1(.*?)(Quantum|                    Dipole Moment|SOLUTE SEEMS)', 
                outfileData, re.DOTALL).group().splitlines()[1:-1]
            
            _parametersSanityCheck(anglesFile)

            anglesFile = [angle[:39].split() for angle in anglesFile if not len(angle)==0 and
                not 'Missing' in angle and 
                not 'Synonym' in angle and
                not 'CHECK' in angle]

        except:

            anglesFile = []

        anglesVariable = []
        anglesAdditional = []

        for angleF in anglesFile:

            atomA = atoms[serial2IndexAtoms[int(angleF[0])]]
            atomB = atoms[serial2IndexAtoms[int(angleF[1])]]
            atomC = atoms[serial2IndexAtoms[int(angleF[2])]]

            if not '_'.join([str(atomA.serial),str(atomB.serial),str(atomC.serial)]) in additionalAnglesInfo: 

                A = Angle(atomC, atomB, atomA, *angleF[3:])
                anglesVariable.append(A)

            else: 

                A = Angle(atomA, atomB, atomC, *angleF[3:])
                anglesAdditional.append(A)

        return anglesVariable, anglesAdditional

    def _getTorsions(atoms, zmatData, outfileData, serial2IndexAtoms):
        """
        Generate torsions from zmatdata and BOSS out file

        Parameters
        ----------
        zmatData : str
            zmat file data
        outfileData : str
            BOSS out file data
        atoms : lst
            List of atoms
        serial2IndexAtoms : Dict
            Dictionary between serial and index atoms

        Returns
        -------
        TorsionsVariable : list
            Variable torsions list
        TorsionsAdditional : list
            Additional torsions list
        """

        zmatTorsionInfoVariables = re.search(r'Variable Dihedrals follow(.*?)Additional Dihedrals follow', zmatData, re.DOTALL).group().splitlines()[1:-1]
        zmatTorsionInfoVariables = [ele.split() for ele in zmatTorsionInfoVariables ]

        zmatTorsionInfoAditional = re.search(r'Additional Dihedrals follow(.*?)Domain Definitions follow', zmatData, re.DOTALL).group().splitlines()[1:-1]
        zmatTorsionInfoAditional = [ele.split() for ele in zmatTorsionInfoAditional ]

        try:

            TorsionsFile = re.search(r'  Angle    Atom  I Type F    V1        V2        V3        V4 (.*?)                Bond Stretching Parameters', outfileData, re.DOTALL).group().splitlines()[1:-1]
            TorsionsFile = [angle[:100].split() for angle in TorsionsFile if not len(angle)==0 and not 'Using' in angle]

        except:

            TorsionsFile = []

        TorsionsVariable = []
        TorsionsAdditional = []

        for i,torsionF in enumerate(TorsionsFile):

            if i < len(zmatTorsionInfoVariables):

                atom = atoms[serial2IndexAtoms[int(torsionF[1])]]
                
                torsionF.insert(2,atom.parentParentParent)
                torsionF.insert(2,atom.parentParent)
                torsionF.insert(2,atom.parent)

                atomA = atoms[serial2IndexAtoms[int(torsionF[1])]]
                atomB = atoms[serial2IndexAtoms[int(torsionF[2])]]
                atomC = atoms[serial2IndexAtoms[int(torsionF[3])]]
                atomD = atoms[serial2IndexAtoms[int(torsionF[4])]]

                A = Torsion(atomA, atomB, atomC, atomD, *torsionF[5:])

                TorsionsVariable.append(A)

            else:

                improper = False

                additionalsIndex = i - len(zmatTorsionInfoVariables)
                
                torsionF.insert(1,zmatTorsionInfoAditional[additionalsIndex][2])
                torsionF.insert(1,zmatTorsionInfoAditional[additionalsIndex][1])
                torsionF.insert(1,zmatTorsionInfoAditional[additionalsIndex][0])

                atomA = atoms[serial2IndexAtoms[int(torsionF[1])]]
                atomB = atoms[serial2IndexAtoms[int(torsionF[2])]]
                atomC = atoms[serial2IndexAtoms[int(torsionF[3])]]
                atomD = atoms[serial2IndexAtoms[int(torsionF[4])]]

                # TODO: CHECK IF OTHER OPTIONS

                if atomA.serial in atomB.bondedAtoms and atomC.serial in atomB.bondedAtoms and atomD.serial in atomB.bondedAtoms: 

                    if (atomB.element == 'C' or atomB.element == 'N') and len(atomB.bondedAtoms) == 3:

                        improper = True 
                    

                A = Torsion(atomA, atomB, atomC, atomD, *torsionF[5:], improper)
                TorsionsAdditional.append(A)

        return TorsionsVariable, TorsionsAdditional


    def _generateBondedAtoms(atoms, bondsVariable, bondsAdditional):
        """
        Generate bonded atoms property in atom class

        Parameters
        ----------
        atoms : list
            Atom list
        bondsVariable : list
            Variable bonds list
        bondsAdditional : list
            Additional bonds list
        """

        bonds = bondsVariable + bondsAdditional

        for atom in atoms:

            bonded = []

            for bond in bonds:
                if atom.serial==bond.atomA.serial: bonded.append(bond.atomB.serial)
                if atom.serial==bond.atomB.serial: bonded.append(bond.atomA.serial)

            atom.bondedAtoms = bonded

    # def isImproper():

    #     pass
        

    def report(self):
        """
        Print the molecule information in the logger
        """

        logger.info('Number of Atoms: '+ str(len(self.atoms)-self.numberOfStructuralDummyAtoms))

        logger.info('Number of Variable Bonds: '+ str(len(self.bondsVariable)))
        logger.info('Number of Additional Bonds: '+ str(len(self.bondsAdditional)))

        logger.info('Number of Variable Angles: '+ str(len(self.anglesVariable)))
        logger.info('Number of Additional Angles: '+ str(len(self.anglesAdditional)))

        logger.info('Number of Variable Torsions: '+ str(len(self.torsionsVariable)))
        logger.info('Number of Additional Torsions: '+ str(len(self.torsionsAdditional)))

        logger.info('Number of excluded Atoms: '+ str(len(self.excludedList)))


        
      
 