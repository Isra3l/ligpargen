"""

Module with utilities
Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

from ..tools import geometry
import os
import shutil
import subprocess
import contextlib
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
import re

import logging
logger = logging.getLogger(__name__)


def translateToceroZcoord(moleculeRDkit):
    """
    Translate the molecule to put the first atom in the origin of the coordinates

    Parameters
    ----------
    moleculeRDkit : RDkit molecule
        An RDkit molecule

    Returns
    -------
    List
        List with the shift value applied to X, Y, Z
    """
    
    from rdkit.Chem import rdMolTransforms

    conf = moleculeRDkit.GetConformer()

    # avoid first atom overlap with dummy 3


    if abs(conf.GetAtomPosition(0).x-1.0)<1e-3 and abs(conf.GetAtomPosition(0).y-1.0)<1e-3 and abs(conf.GetAtomPosition(0).z-0.0)<1e-3:
            

        shiftX = conf.GetAtomPosition(0).x - 1.0
        shiftY = conf.GetAtomPosition(0).y - 1.0
        shiftZ = conf.GetAtomPosition(0).z

        translationMatrix  = np.array(  [[1, 0, 0, -shiftX],
                                        [0, 1, 0, -shiftY],
                                        [0, 0, 1, -shiftZ],
                                        [0, 0, 0, 1]], dtype=np.double)

        rdMolTransforms.TransformConformer(conf, translationMatrix)

    else: 

        shiftX = 0.0
        shiftY = 0.0
        shiftZ = 0.0

    return [shiftX, shiftY, shiftZ]


def generateFolder(folderName):
    """
    Generate folder to place files and return full path

    Parameters
    ----------
    folderName : str
        Folder name to be created

    Returns
    -------
    [str]
        Full path of the new created folder

    """

    newfolderName = os.path.abspath(folderName)

    if os.getcwd()!=newfolderName: 

        if os.path.exists(newfolderName): shutil.rmtree(newfolderName, ignore_errors=True)
        
        os.makedirs(newfolderName)

    return newfolderName


def generateRDkitMolecule(ifile, smile, workdir, molname, debug = False):
    """
    Generate the RDkit molecule including reordering atoms in the molecule to avoid BOSS
    conflicts

    Parameters
    ----------
    ifile : str
        Input file name
    smile : str
        Input SMILE string
    workdir : str
        Working directory
    molname : str
        Molecule name
    debug : bool, optional
        Debug output in logger, by default False

    Returns
    -------
    RDkit molecule
        An RDkit molecule object
    """

    if smile: logger.info('Generating molecule from input SMILES: '+ smile)
    else: logger.info('Parsing molecule from input: '+ ifile)

    sfile   = os.path.join(workdir, molname+"-debug.pdb")

    atomsNameOriginal, residueNameOriginal = [], []

    if smile: 

        molecule = Chem.MolFromSmiles(smile)
        molecule = Chem.AddHs(molecule)
        AllChem.EmbedMolecule(molecule,randomSeed=0xf00d)

    else:

        if ifile.lower().endswith('.pdb'): 
            
            molecule = Chem.MolFromPDBFile(ifile,removeHs=False)

            atomsNameOriginal, residueNameOriginal = getAtomsNameAndMolFromPDB(ifile)

        elif ifile.lower().endswith('.mol2'): 
            
            molecule = Chem.MolFromMol2File(ifile,removeHs=False)

            atomsNameOriginal, residueNameOriginal = getAtomsNameAndMolFromMOL2(ifile)

        elif ifile.lower().endswith('.mol'): molecule = Chem.MolFromMolFile(ifile,removeHs=False)
        else:

            babel = shutil.which("obabel")

            if babel != None:

                subprocess.run([babel, ifile, "-opdb", "-O", sfile], stdin =subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                molecule = Chem.MolFromPDBFile(sfile,removeHs=False)

            else:

                logger.error('obale NOT found in the system')
                logger.error('Input file format CAN NOT be process: '+ ifile)
                logger.error('Please install open BABEL in your system or provide an input in PDB or MOL file')
                exit()
    
    atomsIndexLstWithRightOrder = buildProperOrderAtomIndexLst(molecule, molname, workdir)

    molecule, newIndexToOriginalIndex = generateNewMoleculeWithProperOrder(molecule, atomsIndexLstWithRightOrder)

    if smile: newIndexToOriginalIndex = {i:i for i in range(len(atomsIndexLstWithRightOrder))}

    if debug: 
    
        Chem.MolToPDBFile(molecule, sfile)

    if isAtomOrderWrong(molecule):

        logger.error('ERROR in atom order!!!')
        exit()

    return molecule, newIndexToOriginalIndex, atomsNameOriginal, residueNameOriginal


def getAtomsNameAndMolFromPDB(ifile):

    logger.info('Retreiving atom names form PDB:' + ifile)

    atomsData = [[line[12:16], line[17:20]] for line in open(ifile) if line.startswith('ATOM') or line.startswith('HETATM')]

    atomsName = [item[0] for item in atomsData]
    resName = atomsData[0][1]

    return atomsName, resName


def guessFormalMoleculeCharge(moleculeRDkit, userCharge):
    """
    Guess formal charge of the input molecules or keep original user defined charge

    Parameters
    ----------
    moleculeRDkit : rdkit molecule object
        RDkit molecule object 
    userCharge : int
        Charge defined by the user

    Returns
    -------
    int
        Final charge of the molecule
    """

    autoCharge = Chem.rdmolops.GetFormalCharge(moleculeRDkit)

    if userCharge == None:

        userCharge = autoCharge

        logger.info('AUTO DETECTED CHARGE: ' + str(autoCharge))

    else:

        if autoCharge != userCharge: 
            
            logger.warning('AUTO DETECTED CHARGE: ' + str(autoCharge) + '  USER INPUT CHARGE: ' + str(userCharge) + '  (NOT MATCHING CHARGES)')
    
        logger.info('MOLECULE CHARGE USED: ' + str(userCharge))

    return userCharge

def guessResname(residueNameOriginal, userResname):
    """
    Define the residue Name if user is not defining it

    Parameters
    ----------
    residueNameOriginal : str
        Residue name found in input file (if there is one)
    userResname : str
        Residue name defined by user

    Returns
    -------
    str
        Residue name
    """

    if userResname == None:

        if len(residueNameOriginal) != 0: userResname = residueNameOriginal
        else: userResname = 'MOL'

        logger.info('AUTO ASSIGNED RESIDUE NAME: '+ userResname)

    else: logger.info('RESIDUE NAME USED: '+ userResname)

    return userResname[:3]


def guessMoleculeName(ifileName, userMolName, label):
    """
    Guess the molecule from the inputs

    Parameters
    ----------
    ifileName : str
        Input file name
    userMolName : str
        Input molecule name
    label : str
        label to be added to the molecule name (such as A or B to have mol_A or mol_B)

    Returns
    -------
    str
        Final generated molecule name
    """

    if userMolName == None:

        if ifileName == None: userMolName = 'mol_' + label
        else: userMolName = os.path.splitext(os.path.basename(ifileName))[0]

        logger.info('AUTO ASSIGNED MOLECULE NAME: '+ userMolName)

    else: 
        
        if label != 'A': userMolName = userMolName + '_B'
        logger.info('MOLECULE NAME USED: '+ userMolName)

    return userMolName



def guessChargeAlgorithm(userChargeAlgorithm, charge):
    """
    Define the charge method applied to generate the molecular charges. CM1A-LBCC just valid
    for neutral molecules.

    Parameters
    ----------
    userChargeAlgorithm : str
        Charge method defined by the user (CM1A or CM1A-LBCC)
    charge : int
        Formal charge of the input molecule

    Returns
    -------
    str
        Charge method (CM1A or CM1A-LBCC)
    """

    if userChargeAlgorithm == None:

        if charge == 0: userChargeAlgorithm = 'CM1A-LBCC'
        else: userChargeAlgorithm = 'CM1A'

        logger.info('AUTO DETECTED CHARGE ALGORITHM: ' + userChargeAlgorithm)

    else: 
        
        if charge != 0 and userChargeAlgorithm == 'CM1A-LBCC': 

            userChargeAlgorithm = 'CM1A'

            logger.warning('INCOMPATIBLE CHARGE ALGORTIHM CM1A-LBCC with charged molecules')
            logger.warning('AUTO ASSIGNED CM1A charge model')

    logger.info('CHARGE ALGORITHM USED: ' + userChargeAlgorithm)

    return userChargeAlgorithm




def getAtomsNameAndMolFromMOL2(ifile):

    logger.info('Retreiving atom names form MOL:' + ifile)

    with open(ifile) as f: data = f.read()

    atomsDat = re.search(r'@<TRIPOS>ATOM(.*?)@<TRIPOS>BOND', data, re.DOTALL).group().splitlines()[1:-1]

    atomsData = [[line[8:12], line[59:62]] for line in atomsDat]

    atomsName = [item[0] for item in atomsData]
    resName = atomsData[0][1]

    return atomsName, resName


def getClosestAtomToTheCenterOfMassIndex(molecule):
    """
    Get closest atom to the center of mass

    Parameters
    ----------
    molecule : RDkit molecule
        RDkit molecule

    Returns
    -------
    atomClosestToCenterOfMass : int
        Index of the atom closest to the center
    """

    conf = molecule.GetConformer(0)

    comPoint = Chem.rdMolTransforms.ComputeCentroid(conf)

    distances = []

    for atom in molecule.GetAtoms():

        if atom.GetSymbol() == 'H': continue  # Skip hydrogens

        atom_position = conf.GetAtomPosition(atom.GetIdx())

        distances.append([geometry.getDistance([atom_position.x, atom_position.y, atom_position.z], comPoint), atom])

    atomClosestToCenterOfMass = sorted(distances, key=lambda x: x[0])[0]

    logger.info('Center of mass position: %6.3f %6.3f %6.3f' % (comPoint.x, comPoint.y, comPoint.z))

    logger.info('Closest atom to the center of mass: ' + atomClosestToCenterOfMass[1].GetSymbol()+ '  ' +str(atomClosestToCenterOfMass[1].GetIdx())
                    + ' with distance : %6.3f ' % atomClosestToCenterOfMass[0] + ' (A) ')

    return atomClosestToCenterOfMass[1].GetIdx() 



def generateNewMoleculeWithProperOrder(molecule, atomsIndexLstWithRightOrder):
    """
    Regenerate the same molecule but using the right atoms order

    Parameters
    ----------
    molecule : RDkit molecule
        RDkit molecule with the original atom order
    atomsIndexLstWithRightOrder : List
        Index atoms list with the right atom order

    Returns
    -------
    newMolecule : RDkit molecule
        New RDkit molecule with right atom order for zmat
    newIndexToOriginalIndex : Dict
        Dictionary with original atom order to new right atom order correspondency
    """

    conf = molecule.GetConformer(0)

    newAtomLabels = {}

    newMoleculeBlock = ''

    newIndexToOriginalIndex = {}

    for i, item in enumerate(atomsIndexLstWithRightOrder): 

        newIndexToOriginalIndex[i] = item

        atom = molecule.GetAtomWithIdx(item)
        atom_position = conf.GetAtomPosition(item)

        try:

            residueInfo = atom.GetPDBResidueInfo()

            resname = residueInfo.GetResidueName()
            atomName = residueInfo.GetName() 

        except:

            element = atom.GetSymbol()

            if not element in newAtomLabels.keys(): newAtomLabels[element] = 1
            else: newAtomLabels[element] += 1

            atomName = element + str(newAtomLabels[element])

            resname = 'MOL'

        line = 'ATOM%7d%5s%4s%6d%12.3f%8.3f%8.3f  1.00  0.00%12s  \n' % (i, atomName, resname, 1,atom_position.x, 
                    atom_position.y, atom_position.z, atom.GetSymbol())
        
        newMoleculeBlock += line


    newMolecule = Chem.MolFromPDBBlock(newMoleculeBlock, removeHs=False)

    return newMolecule, newIndexToOriginalIndex


def isAtomOrderWrong(moleculeA):
    """
    Check if atoms in the molecule are in a proper order. 
    So, each atom position can be generated from previous atoms.

    Parameters
    ----------
    moleculeA : RDkit molecule
        An RDkit molecule object

    Returns
    -------
    Bool
        True if there is a wrong order.
    """


    bonds = sorted([sorted([x.GetBeginAtomIdx(), x.GetEndAtomIdx()],reverse=True) for x in moleculeA.GetBonds()],key=lambda x: x[0])

    bondsAtomA = list(dict.fromkeys([bond[0] for bond in bonds]))

    missingAtomsInBondLst = set(range(min(bondsAtomA), max(bondsAtomA)+1)) - set(bondsAtomA)

    if len(missingAtomsInBondLst) >0: return True
    else: return False


def buildProperOrderAtomIndexLst(molecule, molname, workdir):
    """
    Generate atoms index list using the closest atom the centre of mass as the first one and considering
    that atoms have to be generated from the previous atoms positions (Zmat condition)

    Parameters
    ----------
    molecule : RDkit molecule
        Rdkit molecule object
    molname : str
        Molecule name
    workdir : str
        Working folder

    Returns
    -------
    List
        Atoms index List in the right order
    """

    heavyAtomsIndexLst = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetSymbol() != 'H']

    firstAtomIndex = getClosestAtomToTheCenterOfMassIndex(molecule)

    atomsIndexLstWithRightOrder = [firstAtomIndex]

    addAtomsToTheList(molecule, atomsIndexLstWithRightOrder, heavyAtomsIndexLst)

    hydrogenAtomsIndexLst = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetSymbol() == 'H']

    addAtomsToTheList(molecule, atomsIndexLstWithRightOrder, hydrogenAtomsIndexLst)

    return atomsIndexLstWithRightOrder

def addAtomsToTheList(molecule, atomsIndexLst, heavyAtomsIndexLst):
    """
    Add atoms to the list using recursivity

    Parameters
    ----------
    molecule : RDkit molecule
        RDkit molecule
    atomsIndexLst : list
        Atoms index list
    heavyAtomsIndexLst : list
        heavy atoms index
    """

    missingAtomsToAdd = []

    for atomIndex in heavyAtomsIndexLst:

        atom = molecule.GetAtomWithIdx(atomIndex)

        neighborsIndex = [x.GetIdx() for x in atom.GetNeighbors()]

        if atomIndex not in atomsIndexLst:

            if any([True for neighIndex in neighborsIndex if neighIndex in atomsIndexLst]):

                atomsIndexLst.append(atom.GetIdx())

            else:

                missingAtomsToAdd.append(atom.GetIdx())

    if len(missingAtomsToAdd)==0: return
    else: addAtomsToTheList(molecule, atomsIndexLst, heavyAtomsIndexLst)


        
def generateLogger(filename, workdir):


    FORMATTER = logging.Formatter("%(asctime)s — %(levelname)s — %(message)s", 
        datefmt='%Y-%m-%d %H:%M:%S')

    # %(funcName)s 

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.ERROR)

    LOG_FILE = os.path.join(workdir, filename+'.log')
    file_handler = logging.FileHandler(LOG_FILE,'w')
    file_handler.setFormatter(FORMATTER)


    logger = logging.getLogger(__name__)
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    logger.setLevel(logging.INFO)

    return logger


@contextlib.contextmanager
def changedir(dirname):
    """
    Context manager to change temporally the working folder

    Parameters
    ----------
    dirname : str
        Folder to work
    """

    currentdir = os.getcwd()
  
    try:
    
        os.chdir(dirname)
        yield

    finally:
        os.chdir(currentdir) 