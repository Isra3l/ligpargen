"""

Module with geometry related functions (compute angles, distances, ...)
Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""


import math
import copy
from ..topology import Atom, Bond, Angle, Torsion

from rdkit.Geometry.rdGeometry import Point3D, ComputeSignedDihedralAngle

import logging
logger = logging.getLogger(__name__)


def generateAngles(molecule, bonds):
    """
    Generate angles form the bonds list

    Parameters
    ----------
    molecule : RDkit molecule
        RDkit molecule
    bonds : list
        Bonds index list

    Returns
    -------
    angles : list
        Index angles atoms list
    """

    angles = []

    for bond in bonds:

        atomAindex = bond[0]
        atomBindex = bond[1]

        atomB = molecule.GetAtomWithIdx(atomBindex)

        bondedToB = [x.GetIdx() for x in atomB.GetNeighbors() if not x.GetIdx() == atomAindex]

        angles.extend([bond+[atom] for atom in bondedToB])

    for bond in bonds:

        atomAindex = bond[0]
        atomBindex = bond[1]

        atomA = molecule.GetAtomWithIdx(atomAindex)

        bondedToA = [x.GetIdx() for x in atomA.GetNeighbors() if not x.GetIdx()== atomBindex]

        angles.extend([[atom]+bond for atom in bondedToA])

    return angles

def generateDihedrals(molecule, angles):
    """
    Generate dihedrals angles form the bonds list

    Parameters
    ----------
    molecule : RDkit molecule
        RDkit molecule
    bonds : list
        Bonds index list

    Returns
    -------
    dihedrals : list
        Index angles atoms list
    """

    dihedrals = []

    for angle in angles:

        atomAindex = angle[0]
        atomBindex = angle[1]
        atomC = molecule.GetAtomWithIdx(angle[2])
        bondedToC = [x.GetIdx() for x in atomC.GetNeighbors() if x.GetIdx() != atomBindex and x.GetIdx() != atomAindex]

        dihedrals.extend([angle+[atom] for atom in bondedToC])

    for angle in angles:
    
        atomCindex = angle[2]
        atomBindex = angle[1]
        atomA = molecule.GetAtomWithIdx(angle[0])
        bondedToA = [x.GetIdx() for x in atomA.GetNeighbors() if x.GetIdx() != atomBindex and x.GetIdx() != atomCindex]

        dihedrals.extend([[atom]+angle for atom in bondedToA])

    return dihedrals

def isVariableCandidate(atomsList):
    """
    Check is atoms list from an atom term is variable or not. 
    It checks if the atoms are in decrease order.

    Parameters
    ----------
    atomsList : list
        Bonded term list of atoms

    Returns
    -------
    variable : bool
        True if are a variable candidate
    """

    for ele in atomsList[1:]:
        if atomsList[0]<ele:
            return False
    
    return True

def removeDuplicatesUsingFirstAtom(variablesExtra):
    """
    Remove duplicates from the list

    Parameters
    ----------
    variablesExtra : list
        List of data

    Returns
    -------
    variablesExtra : list
        List without duplicates
    """

    tmp = set()

    tmpLst = []

    for item in variablesExtra:
        if item[0] in tmp: continue
        
        tmpLst.append(item)
        tmp.add(item[0])

    return tmpLst



def splitIntoVariableAndAdditionalFromDihedrals(bondedList, variableDihedrals):
    """Function to split between variable and additional the bonds and angles using the dihedrals

    Args:
        bondedList (lst): List of bond/angles atom indexes 
        variableDihedrals (lst): List of variable dihedrals 

    Returns:
        variables, additional: Two list with bonded terms splitted
    """

    numberOfAtomsInBondedLst = len(bondedList[0])

    variables = [dihedral[0:numberOfAtomsInBondedLst] for dihedral in variableDihedrals]

    if len(variables)>0: 
        
        smallestAtomIndex = variables[0][0]

        variablesExtra = list(dict.fromkeys([tuple(bondedTerm) for bondedTerm in bondedList 
            if isVariableCandidate(bondedTerm) and bondedTerm[0] < smallestAtomIndex]))

    else: 
        
        variablesExtra = list(dict.fromkeys([tuple(bondedTerm) for bondedTerm in bondedList 
            if isVariableCandidate(bondedTerm)]))

    variablesExtra = [list(bondedTerm) for bondedTerm in variablesExtra]

    # remove first atom duplicates

    variablesExtra = removeDuplicatesUsingFirstAtom(variablesExtra)

    variables.extend(variablesExtra)

    variables = sorted(variables, key=lambda x: x[0])

    additional = [bondedTerm for bondedTerm in bondedList if not bondedTerm in variables]

    return variables, additional


def splitIntoVariableAndAdditional(bondedList):
    """
    Split bonded terms list into variable and additional angles

    Parameters
    ----------
    bondedList : list
        Bonded terms lists

    Returns
    -------
    variables : list
        Variables bonded terms list
    additionals : list
        Additional bonded terms list
    """

    variableAtomsSelected = set()

    variables = []
    additionals = []

    for atomsList in bondedList:

        firstAtom = atomsList[0]

        if isVariableCandidate(atomsList) and not firstAtom in variableAtomsSelected: 
            variables.append(atomsList)
            variableAtomsSelected.add(atomsList[0]) 
        else:
            additionals.append(atomsList)

    return variables, additionals

def removeDuplicatesInAdditional(variable, additional):
    """Remove duplicated in additional bonded terms checking the variables

    Args:
        variable (lst): Variable bonded terms
        additional (lst): Additional bonded terms

    Returns:
        lst : Aditional bonded term list without duplicates
    """
    
    newAdditionalLst = []

    totalBondedTerms = copy.deepcopy(variable)

    for bondedTermA in additional:

        found = False

        for bondedTermIntotal in totalBondedTerms:

            bondedTermInTotalBreversed = bondedTermIntotal[::-1]

            if bondedTermA==bondedTermIntotal or bondedTermA==bondedTermInTotalBreversed: 
                found = True
                break

        if not found:
            newAdditionalLst.append(bondedTermA)
            totalBondedTerms.append(bondedTermA)

    return newAdditionalLst

def generateImproperDihedrals(molecule):

    impropers = []

    for atom in molecule.GetAtoms():

        atomNeighbors = atom.GetNeighbors()

        if len(atomNeighbors)==3:

            atomIndex = atom.GetIdx()
            atomNeighborsIndex = [atomA.GetIdx() for atomA in atomNeighbors]

            imp = [atomNeighborsIndex[2],atomNeighborsIndex[1], atomIndex, atomNeighborsIndex[0]]

            impropers.append(imp)

    return impropers

def generateBondAnglesAndDihedrals(molecule, dummyAtomsShift):
    """
    Generate Bond, Angles and Torsion List from RDkit molecules

    Parameters
    ----------
    molecule : RDkit molecule
        An RDkit molecule
    dummyAtomsShift : int
        Atoms shift produced by the addition of structural dummy atoms

    Returns
    -------
    variableBonds : List
        Bond list using atom index
    additionalBonds : List
        Bond list using atom index
    variableAngles : List
        Angle list using atom index
    additionalAngles : List
        Angle list using atom index
    variableDihedrals : List
        Angle list using atom index
    additionalDihedrals : List
        Angle list using atom index
    
    """

    # Get all bonds in the molecule

    bonds = sorted([sorted([x.GetBeginAtomIdx(), x.GetEndAtomIdx()],reverse=True) for x in molecule.GetBonds()], key=lambda x: x[0])

    angles = generateAngles(molecule, bonds)

    dihedrals = generateDihedrals(molecule, angles)

    variableDihedrals, additionalDihedrals = splitIntoVariableAndAdditional(dihedrals)
    variableAngles, additionalAngles = splitIntoVariableAndAdditionalFromDihedrals(angles, variableDihedrals)
    variableBonds, additionalBonds = splitIntoVariableAndAdditionalFromDihedrals(bonds, variableDihedrals)
    

    improperDihedrals = generateImproperDihedrals(molecule)  # Impropers are always considered additional

    additionalAngles = removeDuplicatesInAdditional(variableAngles, additionalAngles)

    additionalDihedrals = removeDuplicatesInAdditional(variableDihedrals, additionalDihedrals)


    def generateFakeAtom(index):
    
        serial = index + dummyAtomsShift 

        return Atom(index, serial,  0.0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 'MOL', 
        1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, '0')


    variableBonds = [Bond(generateFakeAtom(bond[0]),generateFakeAtom(bond[1]),0.0,0.0) for bond in variableBonds]
    additionalBonds = [Bond(generateFakeAtom(bond[0]),generateFakeAtom(bond[1]),0.0,0.0) for bond in additionalBonds]

    variableAngles = [Angle(generateFakeAtom(angle[0]),generateFakeAtom(angle[1]),generateFakeAtom(angle[2]),0.0,0.0) 
        for angle in variableAngles]
    additionalAngles = [Angle(generateFakeAtom(angle[0]),generateFakeAtom(angle[1]),generateFakeAtom(angle[2]),0.0,0.0) 
        for angle in additionalAngles]

    variableDihedrals = [Torsion(generateFakeAtom(dihedral[0]),generateFakeAtom(dihedral[1]),generateFakeAtom(dihedral[2]),
    generateFakeAtom(dihedral[3]),-1.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) for dihedral in variableDihedrals]

    additionalDihedrals = [Torsion(generateFakeAtom(dihedral[0]),generateFakeAtom(dihedral[1]),generateFakeAtom(dihedral[2]),
    generateFakeAtom(dihedral[3]),-1.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) for dihedral in additionalDihedrals]

    improperDihedrals = [Torsion(generateFakeAtom(dihedral[0]),generateFakeAtom(dihedral[1]),generateFakeAtom(dihedral[2]),
    generateFakeAtom(dihedral[3]),-1.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, True) for dihedral in improperDihedrals]

    additionalDihedrals.extend(improperDihedrals)


    return variableBonds, additionalBonds, variableAngles, additionalAngles, variableDihedrals, additionalDihedrals


def getAngle(a, b, c):
    """
    Compute angle between a, b, c positions

    Parameters
    ----------
    a : Point3D
        Position A
    b : Point3D
        Position B
    c : Point3D
        Position C

    Returns
    -------
    angle : float
        Angle between a, b, c in degrees
    """

    ab = b - a
    cb = b - c

    return math.degrees(Point3D.AngleTo(ab, cb))

def getDihedral(a, b, c, d):
    """
    Compute angle between a, b, c, d positions

    Parameters
    ----------
    a : Point3D
        Position A
    b : Point3D
        Position B
    c : Point3D
        Position C
    d : Point3D
        Position D

    Returns
    -------
    angle : float
        Angle between a, b, c in degrees
    """

    return math.degrees(ComputeSignedDihedralAngle(a,b,c,d))


def getCartesiansFromInternals(r0, angle, dihedral, posParent, posParentParent, posParentParentParent):
    """
    Generate cartesian coordinates of an atom from internal coordinates

    Parameters
    ----------
    r0 : [type]
        Distance atom to parent
    angle : float
        Angle between atom and parent and grand parent
    dihedral : float
        Dihedral angle between atom and parent, grand and grand grand parents
    posParent : [float,float,float]
        Position parent atom
    posParentParent : [float,float,float]
        Position grand parent atom
    posParentParentParent : [float,float,float]
        Position grand grand parent atom

    Returns
    -------
    position : [float,float,float]
        Position atom
    """

    CX = posParent[0]
    CY = posParent[1]
    CZ = posParent[2]

    BX = posParentParent[0]
    BY = posParentParent[1]
    BZ = posParentParent[2]

    AX = posParentParentParent[0]
    AY = posParentParentParent[1]
    AZ = posParentParentParent[2]

    st = math.sin(angle)
    xpd = r0*st*math.cos(dihedral)
    ypd = r0*st*math.sin(dihedral)
    zpd = -r0*math.cos(angle)

    vca = [0.0,0.0,0.0]
    vcb = [0.0,0.0,0.0]

    vca[0] = AX-CX
    vca[1] = AY-CY
    vca[2] = AZ-CZ

    vcb[0] = BX-CX
    vcb[1] = BY-CY
    vcb[2] = BZ-CZ

    YP = crossProduct3Normalized(vca,vcb)
    XP = crossProduct3Normalized(vcb,YP)
    ZP = crossProduct3Normalized(XP,YP)

    DX = XP[0]*xpd+YP[0]*ypd+ZP[0]*zpd+CX
    DY = XP[1]*xpd+YP[1]*ypd+ZP[1]*zpd+CY
    DZ = XP[2]*xpd+YP[2]*ypd+ZP[2]*zpd+CZ

    return DX,DY,DZ

def crossProduct3Normalized(a, b):
    
    ab = crossProduct3(a, b)

    c = 1.0/math.sqrt(module2(ab))

    return [ab[0]*c, ab[1]*c, ab[2]*c]

def crossProduct3(a, b):
    
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def module2(vect):
    
    return vect[0]**2+vect[1]**2+vect[2]**2


def getDistance(posA, posB):
    """
    Compute distance between two cartesian positions

    Parameters
    ----------
    posA : list[float, float, float]
        Cartesian position of atom A
    posB : list[float, float, float]
        Cartesian position of atom B

    Returns
    -------
    distance : float
        Distance between points
    """
    
    dx = posA[0]-posB[0]
    dy = posA[1]-posB[1]
    dz = posA[2]-posB[2]

    return math.sqrt(dx**2 + dy**2 + dz**2)