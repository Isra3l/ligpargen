"""

Module to generate hybrid topologies (single and dual topology)

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

import sys
import os
import copy
import string
import collections
from ..tools import geometry
import math

from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Geometry.rdGeometry as rdGeometry

import logging
logger = logging.getLogger(__name__)



def alignMolecules(molA, molB, workdir, debug):
    """
    Align molecule B to molecule A and get atom correspondency between molecules

    Parameters
    ----------
    molA : Molecule Object
        Molecule A
    molB : Molecule Object
        Molecule B
    workdir : str
        Working folder
    logger : logging
        Print log outputs
    debug : bool
        if True, intermediate files are kept

    Returns
    -------
    Dict, list
        BtoA is the index correspondecy between atoms from A and B
        umatchB are molecule B atoms indexes not in A
    """

    from rdkit.Chem.rdFMCS import FindMCS, AtomCompare, BondCompare
      
    _fmcs_params = dict(maximizeBonds=False, threshold=1.0, timeout=60,
                        verbose=False, matchValences=True,
                        ringMatchesRingOnly=True, completeRingsOnly=True,
                        atomCompare=AtomCompare.CompareAny,
                        bondCompare=BondCompare.CompareAny)

    try:
        mcs = FindMCS([molA, molB], **_fmcs_params)
    except ValueError:
        print('\n Max Common Substructure calculation \n failed for this molecule!! \n Please be judicious ')
        logger.error('Max Common Substructure calculation failed for this molecule!! --- Please be judicious')
        sys.exit()

    
    core = Chem.MolFromSmarts(mcs.smartsString)
    matchA = molA.GetSubstructMatch(core)
    matchB = molB.GetSubstructMatch(core)

    if debug:

        Chem.MolToPDBFile(molA, os.path.join(workdir,'molecule_A.pdb'))
        Chem.MolToPDBFile(molB, os.path.join(workdir,'molecule_B.pdb'))

    AllChem.AlignMol(molB, molA, atomMap=list(zip(matchB, matchA)))

    if debug: 
        
        Chem.MolToPDBFile(molB, os.path.join(workdir,'molecule_B_superposedToA.pdb'))

    
    ## Realign the molecule based on MCSS so as to have A and B
    ## States at the same place in space  
    BtoA = {j: i for i, j in zip(matchA, matchB)}
    umatchA = [i for i in range(molA.GetNumAtoms()) if i not in matchA]
    umatchB = [i for i in range(molB.GetNumAtoms()) if i not in matchB]
    score = 2 * core.GetNumHeavyAtoms() / (molA.GetNumHeavyAtoms() + molB.GetNumHeavyAtoms())

    score = 100* score

    logger.info('Overlap between heavy atoms of molecule A and B: %3.2f' % (score) + '%  100*(2*coreHeavyAtoms/(heavyAtomsA+heavyAtomsB))')
    logger.info('Number of heavy atoms in common substruture: %d' % (core.GetNumHeavyAtoms()))
    logger.info('Number of match atoms between these two molecules: %d' % (len(matchA)))
    logger.info('Number of unmatch atoms in molecule A: %d' % (len(umatchA)))
    logger.info('Number of unmatch atoms in molecule B: %d' % (len(umatchB)))

    logger.info('Matched atoms positions in molecule A:  %s' % ('  '.join([str(atomNumber+1) for atomNumber in matchA])))
    logger.info('Matched atoms positions in molecule B:  %s' % ('  '.join([str(atomNumber+1) for atomNumber in matchB])))

    logger.info('Unmatched atoms positions in molecule A:  %s' % ('  '.join([str(atomNumber+1) for atomNumber in umatchA])))
    logger.info('Unmatched atoms positions in molecule B:  %s' % ('  '.join([str(atomNumber+1) for atomNumber in umatchB])))

    return BtoA, umatchB

def generateMoleculeAB_singleTopology(moleculeA, moleculeB, BtoAIndexCorrespondency, umatchB):
    """
    Generate a single topology molecule from molecule A to molecule B

    Parameters
    ----------
    moleculeA : Molecule object
        Molecule A
    moleculeB : Molecule object
        Molecule B
    BtoAIndexCorrespondency : dict
        Molecule B atoms correspondency to molecule A index
    umatchB : list
        Index list of atoms in molecule B without correspondency with A
    logger : logging
        logger to print

    Returns
    -------
    Molecule Object
        Molecule AB with single topology
    """
    logger.info('Generating single topology AB molecule')

    moleculeAB = copy.deepcopy(moleculeA)

    moleculeAB.alchemicalTransformation = True

    serialShift = moleculeA.numberOfStructuralDummyAtoms + 1
    AtoBserialCorrespondency = {v+serialShift: k+serialShift for k, v in BtoAIndexCorrespondency.items()}
    BtoAserialCorrespondency = {k+serialShift: v+serialShift for k, v in BtoAIndexCorrespondency.items()}

    # add dummies

    dummies = {0:0, 1:1,2:2}
    AtoBserialCorrespondency.update(dummies)
    BtoAserialCorrespondency.update(dummies)

    umatchBSerial = [item+serialShift for item in umatchB]

    generateABAtoms(moleculeAB, moleculeB, BtoAIndexCorrespondency, BtoAserialCorrespondency, umatchBSerial)

    generateABbondingTerms(moleculeAB.atoms, moleculeAB.bondsVariable, moleculeB.bondsVariable, AtoBserialCorrespondency, BtoAserialCorrespondency, umatchBSerial)
    generateABbondingTerms(moleculeAB.atoms, moleculeAB.bondsAdditional, moleculeB.bondsAdditional, AtoBserialCorrespondency, BtoAserialCorrespondency, umatchBSerial)

    generateABbondingTerms(moleculeAB.atoms, moleculeAB.anglesVariable, moleculeB.anglesVariable, AtoBserialCorrespondency, BtoAserialCorrespondency, umatchBSerial)
    generateABbondingTerms(moleculeAB.atoms, moleculeAB.anglesAdditional, moleculeB.anglesAdditional, AtoBserialCorrespondency, BtoAserialCorrespondency, umatchBSerial)

    generateABbondingTerms(moleculeAB.atoms, moleculeAB.torsionsVariable, moleculeB.torsionsVariable, AtoBserialCorrespondency, BtoAserialCorrespondency, umatchBSerial)
    generateABbondingTerms(moleculeAB.atoms, moleculeAB.torsionsAdditional, moleculeB.torsionsAdditional, AtoBserialCorrespondency, BtoAserialCorrespondency, umatchBSerial)

    generateGeometryVariations(moleculeAB)
    
    generateExcludedAtomLst(moleculeAB)

    updateCartesianCoordinates(moleculeAB)

    return moleculeAB



def generateMoleculeAB_dualTopology(moleculeA, moleculeB):
    """
    Generate a dual topology molecule from molecule A to molecule B

    Parameters
    ----------
    moleculeA : Molecule object
        Molecule A
    moleculeB : Molecule object
        Molecule B
    BtoAIndexCorrespondency : dict
        Molecule B atoms correspondency to molecule A index
    umatchB : list
        Index list of atoms in molecule B without correspondency with A
    logger : logging
        logger to print

    Returns
    -------
    Molecule Object
        Molecule AB with single topology
    """

    logger.info('Generating dual topology AB molecule')

    moleculeAB = copy.deepcopy(moleculeA)

    moleculeAB.alchemicalTransformation = True

    generateABAtomsDual(moleculeAB, moleculeB)

    for bondTermB in moleculeB.bondsVariable: moleculeAB.bondsVariable.append(bondTermB)
    for bondTermB in moleculeB.bondsAdditional: moleculeAB.bondsAdditional.append(bondTermB)

    for bondTermB in moleculeB.anglesVariable: moleculeAB.anglesVariable.append(bondTermB)
    for bondTermB in moleculeB.anglesAdditional: moleculeAB.anglesAdditional.append(bondTermB)

    for bondTermB in moleculeB.torsionsVariable: moleculeAB.torsionsVariable.append(bondTermB)
    for bondTermB in moleculeB.torsionsAdditional: moleculeAB.torsionsAdditional.append(bondTermB)

    generateExcludedAtomLstDual(moleculeAB, moleculeA)

    return moleculeAB



def updateCartesianCoordinates(moleculeAB):
    """
    Update cartesian coordinates of a molecule

    Parameters
    ----------
    moleculeAB : Molecule Object
        Molecule
    """

    XYZ = [rdGeometry.Point3D(atom.x, atom.y, atom.z) for atom in moleculeAB.atoms]
    neighbors = [[atom.parent, atom.parentParent, atom.parentParentParent]for atom in moleculeAB.atoms]

    for atom in moleculeAB.atoms[moleculeAB.numberOfStructuralDummyAtoms+1:]:
    
        posAtom_r = atom.r
        posAtom_angle = math.radians(atom.angle)
        posAtom_dihedral = math.radians(atom.dihedral)

        atomParent = moleculeAB.atoms[atom.parent-1]
        posParent = [atomParent.x, atomParent.y, atomParent.z]

        atomParentParent = moleculeAB.atoms[atom.parentParent-1]
        posParentParent = [atomParentParent.x, atomParentParent.y, atomParentParent.z]

        atomParentParentParent = moleculeAB.atoms[atom.parentParentParent-1]
        posParentParentParent = [atomParentParentParent.x, atomParentParentParent.y, atomParentParentParent.z]

        x, y, z = geometry.getCartesiansFromInternals(posAtom_r, posAtom_angle, posAtom_dihedral, posParent, posParentParent, posParentParentParent)

        atom.x = x
        atom.y = y
        atom.z = z


def getNeighborAtomsUpTo13(dummyAtomA, bonds, angles):
    """
    Get atom neighbors up to 1-3 connections

    Parameters
    ----------
    dummyAtomA : Atom Object
        Atom to determine the neighbors up to 1-3
    bonds : List of Bond Objects
        Bonds in molecule 
    angles : List of Angle Objects
        Angles in molecule

    Returns
    -------
    list
        Neighbors atom serials
    """

    neighbors = [[bond.atomA.serial, bond.atomB.serial] for bond in bonds if dummyAtomA.serial==bond.atomA.serial or dummyAtomA.serial==bond.atomB.serial]


    neighbors += [[angle.atomA.serial, angle.atomB.serial, angle.atomC.serial] for angle in angles if dummyAtomA.serial==angle.atomA.serial 
        or dummyAtomA.serial==angle.atomB.serial or dummyAtomA.serial==angle.atomC.serial]

    neighbors = set([atomSerial for subLst in neighbors for atomSerial in subLst if atomSerial!=dummyAtomA.serial])

    return neighbors

def generateExcludedAtomLstDual(moleculeAB, moleculeA):
    """
    Generate the excluded atom list in dual topology

    Parameters
    ----------
    moleculeAB : Molecule object
        Molecule dual AB
    moleculeA : Molecule object
        Molecule A (reference)
    """

    firstAtomMoleculeA = moleculeA.atoms[moleculeA.numberOfStructuralDummyAtoms].serial
    lastAtomMoleculeA = len(moleculeA.atoms)

    firstAtomMoleculeB = len(moleculeA.atoms)+1
    lastAtomMoleculeB = len(moleculeAB.atoms)

    moleculeAB.excludedList = [[firstAtomMoleculeA, lastAtomMoleculeA, firstAtomMoleculeB, lastAtomMoleculeB]]

    excludeList = []

    for atomAserial in range(firstAtomMoleculeA - moleculeA.numberOfStructuralDummyAtoms, lastAtomMoleculeA - moleculeA.numberOfStructuralDummyAtoms + 1):
        for atomBserial in range(firstAtomMoleculeB - moleculeA.numberOfStructuralDummyAtoms, lastAtomMoleculeB -moleculeA.numberOfStructuralDummyAtoms + 1):
            excludeList.append([atomAserial, atomBserial])

    excludeListDict_1 = collections.defaultdict(list)
    for atomA,atomB in excludeList: excludeListDict_1[atomA].append(atomB)

    moleculeAB.excludedListDict = excludeListDict_1


def generateExcludedAtomLst(moleculeAB):
    """
    Generate the excluded atom list in single topology

    Parameters
    ----------
    moleculeAB : Molecule object
        Molecule dual AB
    """


    excludeList = []

    dummyAtomsInStateA = [atom for atom in moleculeAB.atoms if atom.typeA == 100]

    if len(dummyAtomsInStateA)!=0:

        totalBonds = moleculeAB.bondsVariable + moleculeAB.bondsAdditional
        totalAngles = moleculeAB.anglesVariable + moleculeAB.anglesAdditional

        dummyAtomsInStateB = [atom for atom in moleculeAB.atoms if atom.typeB == 100]

        for dummyAtomA in dummyAtomsInStateA:

            atomsUpTo13Neighbors = getNeighborAtomsUpTo13(dummyAtomA, totalBonds, totalAngles)

            for dummyAtomB in dummyAtomsInStateB:

                if dummyAtomB.serial not in atomsUpTo13Neighbors: 
                    
                    excludeList.append([dummyAtomA.serial, dummyAtomB.serial])

    # create Ranges

    excludeListDict_1 = collections.defaultdict(list)
    for atomA,atomB in excludeList: excludeListDict_1[atomA].append(atomB)

    moleculeAB.excludedListDict = excludeListDict_1

    newExcludedLst_1 = getExcludedAtomsRanges(excludeListDict_1)

    excludeListDict_2 = collections.defaultdict(list)
    for atomA,atomB,atomC,atomD in newExcludedLst_1: excludeListDict_2[(atomC,atomD)].append(atomA)

    newExcludedLst_2 = getExcludedAtomsRanges(excludeListDict_2, False)
    
    moleculeAB.excludedList = newExcludedLst_2


def getExcludedAtomsRanges(excludeListDict, first = True):
    """
    Generate exclude range atom list

    Parameters
    ----------
    excludeListDict : collections.defaultdict
        Exclude atom list
    first : bool, optional
        [description], by default True

    Returns
    -------
    lst
        Range excluded atom list 
    """

    newExcludedLst = []

    for atomA, atoms in excludeListDict.items():

        nums = sorted(set(atoms))
        gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
        edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
        rangeAtoms = list(zip(edges, edges))

        for rangeAtom in rangeAtoms:

            if first: newExcludedLst.append([atomA,atomA,rangeAtom[0],rangeAtom[1]])
            else: newExcludedLst.append([rangeAtom[0],rangeAtom[1],atomA[0],atomA[1]])

    return newExcludedLst

def generateGeometryVariations(moleculeAB):
    """
    Generate geometry variations for Zmat

    Parameters
    ----------
    moleculeAB : Molecule object
        Molecule 
    """

    for bond in moleculeAB.bondsVariable:

        if bond.R0 != bond.R0_B:

            moleculeAB.geometryVariations.append('%4d%4d%12.6f\n' % (bond.atomA.serial,1,bond.R0_B)) 

    for atom in moleculeAB.atoms: 

        if atom.typeA==100: atom.r = 0.3


def generateABbondingTerms(atomsAB, bondingList_AB, bondingList_B, AtoBserialCorrespondency, BtoAserialCorrespondency, umatchBSerial):
    """
    Generate bonding (bonds, angles and torsions) list in single topology

    Parameters
    ----------
    atomsAB : Atom lst
        Atom list of molecule AB
    bondingList_AB : lst
        Bonded list to build
    bondingList_B : lst
        Bonded list in molecule B
    AtoBserialCorrespondency : dict
        Atom A serials correspondency to atom B serials
    BtoAserialCorrespondency : dict
        Atom B serials correspondency to atom A serials
    umatchBSerial : lst
        Atom serials of molecule B not present in molecule A
    """

    for bondTermA in bondingList_AB:

        atomsBondTermA = bondTermA.getAtoms()

        alchemicalTransformation = False

        alchemicalTransformation = len([ True for atom in atomsBondTermA if atom.atomTypeOPLS != atom.atomTypeOPLS_B])!=0


        if alchemicalTransformation:

            try:

                bondinB = [AtoBserialCorrespondency[atom.serial] for atom in atomsBondTermA]

                bondinBreversed = bondinB[::-1]

                for bondB in bondingList_B:

                    bondBserial = [atom.serial for atom in bondB.getAtoms()]

                    if bondinB == bondBserial or bondinBreversed == bondBserial:

                        if bondTermA.__class__.__name__=='Bond': 

                            bondTermA.R0_B = bondB.R0
                            bondTermA.K0_B = bondB.K0

                        if bondTermA.__class__.__name__=='Angle': 

                            bondTermA.angle0_B = bondB.angle0
                            bondTermA.K0_B = bondB.K0

                        if bondTermA.__class__.__name__=='Torsion': 

                            bondTermA.V1_B = bondB.V1
                            bondTermA.V2_B = bondB.V2
                            bondTermA.V3_B = bondB.V3
                            bondTermA.V4_B = bondB.V4

                            bondTermA.typeFinal = bondB.typeInitial


            except:

                if bondTermA.__class__.__name__=='Bond': bondTermA.R0_B = 0.3


    for bondB in bondingList_B:

        umatchAtomInBondB = [True for atomBSerial in umatchBSerial if atomBSerial in [atom.serial for atom in bondB.getAtoms()]]
    
        if len(umatchAtomInBondB)!=0:
        
            newBond = copy.deepcopy(bondB)

            if bondB.__class__.__name__=='Bond': 

                newBond.atomA = atomsAB[BtoAserialCorrespondency[newBond.atomA.serial]-1]
                newBond.atomB = atomsAB[BtoAserialCorrespondency[newBond.atomB.serial]-1]

                newBond.R0 = 0.3

            if bondB.__class__.__name__=='Angle': 

                newBond.atomA = atomsAB[BtoAserialCorrespondency[newBond.atomA.serial]-1]
                newBond.atomB = atomsAB[BtoAserialCorrespondency[newBond.atomB.serial]-1]
                newBond.atomC = atomsAB[BtoAserialCorrespondency[newBond.atomC.serial]-1]


            if bondB.__class__.__name__=='Torsion': 

                newBond.atomA = atomsAB[BtoAserialCorrespondency[newBond.atomA.serial]-1]
                newBond.atomB = atomsAB[BtoAserialCorrespondency[newBond.atomB.serial]-1]
                newBond.atomC = atomsAB[BtoAserialCorrespondency[newBond.atomC.serial]-1]
                newBond.atomD = atomsAB[BtoAserialCorrespondency[newBond.atomD.serial]-1]


                newBond.typeFinal = bondB.typeInitial


            bondingList_AB.append(newBond)



def generateABAtomsDual(moleculeAB, moleculeB):
    """
    Generate atoms of the hybrid molecule in dual topology

    Parameters
    ----------
    moleculeAB : Molecule object
        Molecule AB
    moleculeB : Molecule object
        Molecule B
    """

    atomsAB = moleculeAB.atoms

    lastAtomSerialInMoleculeAWithoutDummies = len(atomsAB) - moleculeAB.numberOfStructuralDummyAtoms

    for atomA in atomsAB:

        if atomA.typeA ==-1: continue # Skip structural dummies

        atomA.charge_B = 0.0
        atomA.epsilon_B = 0.0
        atomA.sigma_B = 0.0

        atomA.typeB = 9000+atomA.typeB

        atomA.type_q_B =  'D100'


    nameSerialIterator = 0

    for atomB in moleculeB.atoms[moleculeAB.numberOfStructuralDummyAtoms:]:
    
        newSerial = atomB.serial + lastAtomSerialInMoleculeAWithoutDummies

        atomB.serial = newSerial
        atomB.serialOriginal = newSerial

        if atomB.parent > 3: atomB.parent += lastAtomSerialInMoleculeAWithoutDummies
        if atomB.parentParent > 3: atomB.parentParent += lastAtomSerialInMoleculeAWithoutDummies
        if atomB.parentParentParent > 3: atomB.parentParentParent += lastAtomSerialInMoleculeAWithoutDummies 

        atomB.charge = 0.0
        atomB.epsilon = 0.0
        atomB.sigma = 0.0

        atomB.type_q =  'D100'

        if len(atomB.element) == 1: 

            atomB.name = atomB.element + str(newSerial)
            atomB.nameOriginal = atomB.element + str(newSerial)
        
        else: 
        
            atomB.name = atomB.element + str(nameSerialIterator)
            atomB.nameOriginal = atomB.element + str(nameSerialIterator)
            nameSerialIterator +=1

        atomB.type_charmm = atomB.element + str(atomB.typeB+100)
        atomB.type_q_B = atomB.element + str(atomB.typeB+100)
        atomB.typeA = atomB.typeB + lastAtomSerialInMoleculeAWithoutDummies
        atomB.typeB = 9000+atomB.typeB + lastAtomSerialInMoleculeAWithoutDummies

        atomB.dual_just_stateB = True

        atomsAB.append(atomB)


def generateABAtoms(moleculeAB, moleculeB, BtoAIndexCorrespondency, BtoAserialCorrespondency, umatchBSerial):
    """
    Generate atoms of the hybrid molecule in single topology

    Parameters
    ----------
    moleculeAB : Molecule object
        Molecule AB
    moleculeB : Molecule object
        Molecule B
    BtoAIndexCorrespondency : dict
        Atom B index correspondency to atom A index
    BtoAserialCorrespondency : dict
        Atom B serials correspondency to atom A serials
    umatchBSerial : list   
        Atom serials of molecule B not present in molecule A
    """

    atomsAB = moleculeAB.atoms

    numberOfDummyAtoms = moleculeAB.numberOfStructuralDummyAtoms
    AtoBindexCorrespondency = {v+numberOfDummyAtoms: k+numberOfDummyAtoms for k, v in BtoAIndexCorrespondency.items()}

    for indexA,atomA in enumerate(atomsAB):

        if atomA.typeA ==-1: continue # Skip structural dummies

        try:

            atomA.typeB = 9000 + atomA.typeB

            atomB = moleculeB.atoms[AtoBindexCorrespondency[indexA]]

            atomA.atomicNumber_B = atomB.atomicNumber
            atomA.atomTypeOPLS_B = atomB.atomTypeOPLS
            atomA.charge_B = atomB.charge
            atomA.sigma_B = atomB.sigma
            atomA.epsilon_B = atomB.epsilon

            atomA.mass_B = atomB.mass

            if len(atomB.element)==1: atomA.type_q_B = atomB.element + str(atomA.typeB-9000+100) # Q does NOT support 5>= length atomtypes
            else: atomA.type_q_B = atomB.element + str(atomA.typeB-9800)

            atomA.type_gmx_B = 'opls_'+str(atomA.typeB)

        except:

            # NOT matched atom in A molecule so ..... future dummy

            atomA.typeB = 100

            atomA.atomicNumber_B = 99
            atomA.atomTypeOPLS_B = 'DM'
            atomA.charge_B = 0.00000
            atomA.sigma_B = 0.00000
            atomA.epsilon_B = 0.0000

            atomA.type_q_B = 'D' + str(atomA.typeB)
            atomA.type_gmx_B = 'opls_DM'


    newSerialpadding = len(atomsAB)+1 # dummies
    newTypepadding = 800+len(atomsAB)-4 # dummies

    for i,indexB in enumerate(umatchBSerial):
        
        atomB =  copy.deepcopy(moleculeB.atoms[indexB-1])

        newSerial = newSerialpadding+i

        BtoAserialCorrespondency[atomB.serial] = newSerial

        atomB.name = 'Du'+string.ascii_uppercase[i]
        atomB.nameOriginal = 'Du'+string.ascii_uppercase[i]

        newTypepadding += 1
        atomB.typeB = 9000+newTypepadding 
        atomB.typeA = 100
        
        atomB.atomicNumber = 99
        atomB.charge = 0.00000
        atomB.sigma = 0.00000
        atomB.epsilon = 0.0000

        atomB.atomTypeOPLS = 'DM'
        atomB.type_q = 'D' + str(newTypepadding)
        atomB.type_gmx = 'opls_'+str(newTypepadding)
        atomB.type_gmx_B = 'opls_'+str(9000+newTypepadding)

        # atomB.mass = 0.0000
        
        atomB.parent = updateParentOfB_Molecule(atomB.parent,  BtoAserialCorrespondency)
        atomB.parentParent = updateParentOfB_Molecule(atomB.parentParent, BtoAserialCorrespondency)
        atomB.parentParentParent = updateParentOfB_Molecule(atomB.parentParentParent, BtoAserialCorrespondency)

        atomB.serial = newSerial

        # atomB.indexOriginal = newSerial 
        atomB.serialOriginal = newSerial
        
        moleculeAB.atoms.append(atomB)


def updateParentOfB_Molecule(serialB,BtoA):

    newSerial = 0

    try:

        newSerial = BtoA[serialB]

    except KeyError:

        print(serialB, BtoA)

        raise KeyError('No key bar in dict foo') from None

    return newSerial
