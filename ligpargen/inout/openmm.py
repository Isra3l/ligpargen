"""

Module with functions to generate OpenMM software inputs (PDB, XML)

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

import os
import math

headerXML = '''<!--       Generated with LigParGen   --> 
<!--        William L. Jorgensen Lab  -->
<!--     Author: israel.cabezadevaca@yale.edu  -->
<!--    OPLS Force Field with CM1A derived Atomic Charges  -->
'''

kcalToKj = 4.184


def printBondedAtoms(bond, shiftdummySerial):
    """Generate bond atoms line

    Parameters
    ----------
    bond : Bond Object
        Bond Object
    shiftdummySerial : int
        Shift to remove structural dummy atoms and get index (should be 4)

    Returns
    -------
    bondLine : str
        Bond line data
    """

    return '<Bond from=\"%d\" to=\"%d\"/>\n' % (bond.atomA.serialOriginal- shiftdummySerial, bond.atomB.serialOriginal - shiftdummySerial)


def printBond(bond):
    """Generate bond line

    Parameters
    ----------
    bond : Bond Object
        Bond Object

    Returns
    -------
    bondLine : str
        Bond line data
    """

    R0 = 0.1*bond.R0
    K0 = bond.K0*836.80

    return '<Bond class1=\"%s\" class2=\"%s\" length=\"%0.6f\" k=\"%0.6f\"/>\n' % (bond.atomA.type_q, bond.atomB.type_q, R0, K0)



def printAngle(angle):
    """Generate angle line

    Parameters
    ----------
    angle : Angle Object
        Angle Object

    Returns
    -------
    angleLine : str
        Angle line data
    """
    
    k0 = angle.K0*8.3680
    angle0 = math.radians(angle.angle0)

    return '<Angle class1=\"%s\" class2=\"%s\" class3=\"%s\" angle=\"%0.6f\" k=\"%0.6f\"/>\n' % \
        (angle.atomA.type_q, angle.atomB.type_q, angle.atomC.type_q, angle0, k0)


def printDihedral(dihedral):
    """Generate dihedral line

    Parameters
    ----------
    dihedral : dihedral Object
        dihedral Object

    Returns
    -------
    dihedralLine : str
        dihedral line data
    """

    V1 = dihedral.V1*0.5*kcalToKj
    V2 = dihedral.V2*0.5*kcalToKj
    V3 = dihedral.V3*0.5*kcalToKj
    V4 = dihedral.V4*0.5*kcalToKj

    label = '<Proper class1=\"%s\" class2=\"%s\" class3=\"%s\" class4=\"%s\" ' % (dihedral.atomA.type_q, dihedral.atomB.type_q, dihedral.atomC.type_q, dihedral.atomD.type_q)
    label += 'k1=\"%0.6f\" k2=\"%0.6f\" k3=\"%0.6f\" k4=\"%0.6f\" ' % (V1, V2, V3, V4)
    label += 'periodicity1=\"1\" periodicity2=\"2\" periodicity3=\"3\" periodicity4=\"4\" '
    label += 'phase1=\"0.00\" phase2=\"3.141592653589793\" phase3=\"0.00\" phase4=\"3.141592653589793\"/>\n'

    return label


def printImproperDihedral(dihedral):
    """Generate improper dihedral line

    Parameters
    ----------
    dihedral : dihedral Object
        dihedral Object

    Returns
    -------
    dihedralLine : str
        Improper dihedral line data
    """

    V1 = dihedral.V1*0.5*kcalToKj
    V2 = dihedral.V2*0.5*kcalToKj
    V3 = dihedral.V3*0.5*kcalToKj
    V4 = dihedral.V4*0.5*kcalToKj

    label = '<Improper class1=\"%s\" class2=\"%s\" class3=\"%s\" class4=\"%s\" ' % (dihedral.atomA.type_q, dihedral.atomB.type_q, dihedral.atomC.type_q, dihedral.atomD.type_q)
    label += 'k1=\"%0.6f\" k2=\"%0.6f\" k3=\"%0.6f\" k4=\"%0.6f\" ' % (V1, V2, V3, V4)
    label += 'periodicity1=\"1\" periodicity2=\"2\" periodicity3=\"3\" periodicity4=\"4\" '
    label += 'phase1=\"0.00\" phase2=\"3.141592653589793\" phase3=\"0.00\" phase4=\"3.141592653589793\"/>\n'

    return label


def writeXML(molecule, xmlFile):
    """Generate XML file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    xmlFile : str
        XML file name
    """

    with open(xmlFile,'w') as ofile:

        ofile.write(headerXML)

        ofile.write('<ForceField>\n')

        ofile.write('<AtomTypes>\n')

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for i, atom in enumerate(atomsToWrite, start =1):

            ofile.write('<Type name=\"%s\" class=\"%s\" element=\"%s\" mass=\"%0.3f\"/>\n' %(atom.type_gmx, atom.type_q, atom.element, atom.mass))

        ofile.write('</AtomTypes>\n')

        ofile.write("<Residues>\n")
        ofile.write("<Residue name=\"%s\">\n" % molecule.residueName)

        for i, atom in enumerate(atomsToWrite, start =1):

            ofile.write("<Atom name=\"%s\" type=\"%s\"/>\n" % (atom.nameOriginal.replace(" ",""), atom.type_gmx))

        shift = molecule.numberOfStructuralDummyAtoms + 1
        
        for bond in molecule.bondsVariable: ofile.write(printBondedAtoms(bond, shift))
        for bond in molecule.bondsAdditional: ofile.write(printBondedAtoms(bond, shift))

        ofile.write('</Residue>\n')
        ofile.write('</Residues>\n')
        
        ofile.write('<HarmonicBondForce>\n')

        for bond in molecule.bondsVariable: ofile.write(printBond(bond))
        for bond in molecule.bondsAdditional: ofile.write(printBond(bond))

        ofile.write('</HarmonicBondForce>\n')

        ofile.write('<HarmonicAngleForce>\n')

        for angle in molecule.anglesVariable: ofile.write(printAngle(angle))
        for angle in molecule.anglesAdditional: ofile.write(printAngle(angle))

        ofile.write('</HarmonicAngleForce>\n')

        ofile.write('<PeriodicTorsionForce>\n')

        for dihedral in molecule.torsionsVariable: 
            if dihedral.improper==False: ofile.write(printDihedral(dihedral))
        for dihedral in molecule.torsionsAdditional: 
            if dihedral.improper==False: ofile.write(printDihedral(dihedral))


        for dihedral in molecule.torsionsVariable: 
            if dihedral.improper==True: ofile.write(printImproperDihedral(dihedral))
        for dihedral in molecule.torsionsAdditional: 
            if dihedral.improper==True: ofile.write(printImproperDihedral(dihedral))

        ofile.write('</PeriodicTorsionForce>\n')

        ofile.write("<NonbondedForce coulomb14scale=\"0.5\" lj14scale=\"0.5\">\n")

        for i, atom in enumerate(atomsToWrite, start =1):
    
            ofile.write('<Atom type=\"%s\" charge=\"%0.6f\" sigma=\"%0.6f\" epsilon=\"%0.6f\"/>\n' % \
                (atom.type_gmx, atom.charge, atom.sigma*0.1, atom.epsilon*kcalToKj))

        ofile.write("</NonbondedForce>\n")
        ofile.write("</ForceField>\n")


def writePDB(molecule, pdbFile):
    """Generate PDB file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    pdbFile : str
        PDB file name
    """
    
    with open(pdbFile, 'w') as ofile:

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for i, atom in enumerate(atomsToWrite, start =1):

            ofile.write('ATOM%7d%5s%4s%6d%12.3f%8.3f%8.3f  1.00  0.00%12s  \n' % (i, atom.nameOriginal, molecule.residueName, 1,atom.x + molecule.shiftX, 
                atom.y + molecule.shiftY, atom.z + molecule.shiftZ, atom.element))


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
    pdbFile : str
        PDB file name
    xmlFile : str
        XML file name
    """

    pdbFile = os.path.join(workdir,molname+'.openmm.pdb')
    xmlFile = os.path.join(workdir,molname+'.openmm.xml')

    return pdbFile, xmlFile


def write(molecule, molName, workdir):

    pdbFile, xmlFile = getFileNames(molName, workdir)

    writeXML(molecule, xmlFile)

    writePDB(molecule, pdbFile)

