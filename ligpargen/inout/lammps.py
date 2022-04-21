"""

Module with functions to generate LAMMPS software inputs (LMP)

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

import os

headerLMP = 'LAMMPS data file generated with LigParGen (israel.cabezadevaca@yale.edu)\n\n'

def writeLMP(molecule, lmpFile):
    """Generate XML file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    lmpFile : str
        LMP file name
    """

    with open(lmpFile,'w') as ofile:

        ofile.write(headerLMP)

        bonds = molecule.bondsVariable + molecule.bondsAdditional
        angles = molecule.anglesVariable + molecule.anglesAdditional
        torsions = molecule.torsionsVariable + molecule.torsionsAdditional
        propers = [torsion for torsion in torsions if not torsion.improper]
        impropers = [torsion for torsion in torsions if torsion.improper]

        xlist = [atom.x for atom in molecule.atoms]
        ylist = [atom.y for atom in molecule.atoms]
        zlist = [atom.z for atom in molecule.atoms]

        max_mol_size = 50


        ofile.write('%8d atoms\n' % len(molecule.atoms[molecule.numberOfStructuralDummyAtoms:]))
        ofile.write('%8d bonds\n' % len(bonds))
        ofile.write('%8d angles\n' % len(angles))
        ofile.write('%8d dihedrals\n' % len(propers))
        ofile.write('%8d impropers\n\n' % len(impropers))

        ofile.write('%8d atom types\n' % len(molecule.atoms))
        ofile.write('%8d bond types\n' % len(bonds))
        ofile.write('%8d angle types\n' % len(angles))
        ofile.write('%8d dihedral types\n' % len(propers))
        ofile.write('%8d improper types\n\n' % len(impropers))

        ofile.write('%12.6f %12.6f xlo xhi\n' % (min(xlist), min(xlist) + max_mol_size))
        ofile.write('%12.6f %12.6f ylo yhi\n' % (min(ylist), min(ylist) + max_mol_size))
        ofile.write('%12.6f %12.6f zlo zhi\n\n' % (min(zlist), min(zlist) + max_mol_size))


        ofile.write('\nMasses\n\n')

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for i, atom in enumerate(atomsToWrite, start =1): 
            ofile.write('%8d %10.3f\n' %( i, atom.mass))

        ofile.write('\nPair Coeffs\n\n')

        for i, atom in enumerate(atomsToWrite, start =1): 
            ofile.write('%8d%11.3f%11.7f\n' %( i, atom.epsilon, atom.sigma))

        ofile.write('\nBond Coeffs\n\n')

        for i, bond in enumerate(bonds, start=1): 
            ofile.write('%8d%11.4f%11.4f\n' %( i, bond.K0, bond.R0))

        ofile.write('\nAngle Coeffs\n\n')

        for i, angle in enumerate(angles, start=1): 
            ofile.write('%8d%11.3f%11.3f\n' %( i, angle.K0, angle.angle0))

        ofile.write('\nDihedral Coeffs\n\n')

        for i, torsion in enumerate(propers, start=1): 
            ofile.write('%8d%11.3f%11.3f%11.3f%11.3f\n' %( i, torsion.V1, torsion.V2, torsion.V3, torsion.V4))

        ofile.write('\nImproper Coeffs\n\n')

        for i, torsion in enumerate(impropers, start=1): 
            ofile.write('%8d%11.3f%8d%8d\n' %( i, torsion.V2/2.0, -1, 2.0)) #TODO: check if V2 or V2/2

        ofile.write('\nAtoms\n\n')

        for i, atom in enumerate(atomsToWrite, start =1): 
            ofile.write('%6d %6d %6d %10.6f %8.5f %8.5f %8.5f\n' %( i, 1, i, atom.charge, atom.x, atom.y, atom.z))

        shift = molecule.numberOfStructuralDummyAtoms

        ofile.write('\nBonds\n\n')

        for i, bond in enumerate(bonds, start=1): 
            ofile.write('%6d %6d %6d %6d\n' %( i, i, bond.atomA.serialOriginal - shift, bond.atomB.serialOriginal - shift))

        ofile.write('\nAngles\n\n')

        for i, angle in enumerate(angles, start=1): 
            ofile.write('%6d %6d %6d %6d %6d\n' %( i, i, angle.atomA.serialOriginal - shift, angle.atomB.serialOriginal - shift, angle.atomC.serialOriginal - shift))

        ofile.write('\nDihedrals\n\n')

        for i, torsion in enumerate(propers, start=1): 
            ofile.write('%6d %6d %6d %6d %6d %6d\n' %( i, i, torsion.atomA.serialOriginal - shift, torsion.atomB.serialOriginal - shift, \
                torsion.atomC.serialOriginal - shift, torsion.atomD.serialOriginal - shift))

        ofile.write('\nImpropers\n\n')

        for i, torsion in enumerate(impropers, start=1): 
            ofile.write('%6d %6d %6d %6d %6d %6d\n' %( i, i, torsion.atomA.serialOriginal - shift, torsion.atomB.serialOriginal - shift, \
                torsion.atomC.serialOriginal - shift, torsion.atomD.serialOriginal - shift))

        ofile.write('\n\n')

# def writePDB(molecule, pdbFile):
#     """Generate PDB file

#     Parameters
#     ----------
#     molecule : molecule class
#         Molecule class
#     pdbFile : str
#         PDB file name
#     """
    
#     with open(pdbFile, 'w') as ofile:

#         atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

#         for i, atom in enumerate(atomsToWrite, start =1):

#             ofile.write('ATOM%7d%5s%4s%6d%12.3f%8.3f%8.3f  1.00  0.00%12s  \n' % (i, atom.nameOriginal, atom.resname, 1,atom.x + molecule.shiftX, 
#                 atom.y + molecule.shiftY, atom.z + molecule.shiftZ, atom.element))


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
    lmpFile : str
        PDB file name
    """

    lmpFile = os.path.join(workdir,molname+'.lammps.lmp')

    return lmpFile


def write(molecule, molName, workdir):

    lmpFile = getFileNames(molName, workdir)

    writeLMP(molecule, lmpFile)


