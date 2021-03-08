"""

Module with functions to generate PQR input (PQR)

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

import os

kcalToKj = 4.184

def writePQR(molecule, pqrFile):
    """Generate PQR file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    pqrFile : str
        PQR file name
    """

    sigmaFactor = 2.0**(1.0/6.0)  # 2^(1/6)

    with open(pqrFile, 'w') as ofile:

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for i, atom in enumerate(atomsToWrite, start =1):

            ofile.write('ATOM  %5d %4s %3s  %4d    %8.3f%8.3f%8.3f%8.4f%7.4f\n' % (i, atom.nameOriginal, molecule.residueName, 1,atom.x + molecule.shiftX, 
                atom.y + molecule.shiftY, atom.z + molecule.shiftZ, atom.charge, (atom.sigma/2.0)*sigmaFactor))


def getFileNames(molname, workdir):
    """Return output file names

    Parameters
    ----------
    molname : str
        Molecule name
    workdir : str
        Working folder path

    Returns
    -------
    pqrFile : str
        PQR file name
    """

    return os.path.join(workdir,molname+'.pqr')


def write(molecule, molName, workdir):

    pqrFile = getFileNames(molName, workdir)

    writePQR(molecule, pqrFile)

